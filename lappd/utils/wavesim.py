from multiprocessing import Pool
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from lappd.lappdevent import LAPPDEvent
from lappd.matching import plot_hitmap, strip_positions, x_to_t
from lappd.utils import lappdcfg as cfg
from lappd.utils.lognormal import LogNormal
from lappd.utils.pdf import ChargeSharePDF
from lappd.utils.pyhist import pyhist
from lappd.utils.roothist import roothist
from scipy.stats import exponnorm

root.gSystem.Load("libRATEvent.so")
strip_positions = strip_positions
nsamples = cfg.NSAMPLES - cfg.NREMOVEDSAMPLES
xs = np.arange(1, nsamples+1) * cfg.NS_PER_SAMPLE
medianx = np.max(xs)
NOISE_RMS = 0.63  # calculated from rms of 10,000 empty pulses in dark rate data
TIME_OFFSET = 10.0  # ns, arbitrary offset for when first photon is registered

# chargex = np.arange(0, 2.0, 0.05)
# chargepdf = exponnorm.pdf(chargex, 0.0001, loc=0.3, scale=0.30)
# plt.plot(chargex, chargepdf)
# plt.xlim(0, 2.0)
# plt.show()
# breakpoint()

# Look at getting the MCPMT + MCPhoton info as well
# FrontEndTime is already smeared by time resolution
# Hit time is time when photon strikes cathode
# FrontEndTime is time after the dynode chain?
# Need to match up MCPMTs with PMTs, GetPMTCount doesn't always give the same
# number for each? -> True for LSD2, not for LSD1


class Hit():

    def __init__(self, pmt, mcpmt) -> None:
        self.pmt = pmt
        self.mcpmt = mcpmt


class Event():

    def __init__(self, tree: root.TTree, pmtinfo, evnumber: int) -> None:
        self.tree = tree
        self.pmtinfo = pmtinfo
        self.evnumber = evnumber
        self.tree.GetEntry(self.evnumber)
        self.ds = tree.ds
        self.evcount = self.ds.GetEVCount()
        self.ev = self.ds.GetEV(0) if self.evcount else None
        self.trigger_time = self.ev.GetCalibratedTriggerTime() if self.ev else None
        self.delta_t = self.ev.GetDeltaT() if self.ev else None
        self.total_charge = self.ev.GetTotalCharge() if self.ev else None
        self.pmtcount = self.ev.GetPMTCount() if self.ev else None
        self.pmts = [self.ev.GetPMT(ipmt) for ipmt in range(
            self.pmtcount)] if self.ev else None
        self.pmt_times = [pmt.GetTime()
                          for pmt in self.pmts] if self.pmts else None
        self.pmt_charges = [pmt.GetCharge()
                            for pmt in self.pmts] if self.pmts else None
        self.pmt_ids = [pmt.GetID()
                        for pmt in self.pmts] if self.pmts else None
        self.pmt_types = [self.get_pmt_type(pmt)
                          for pmt in self.pmts] if self.pmts else None
        self.mc = tree.ds.GetMC() if self.ev else None
        self.mcpmtcount = self.mc.GetMCPMTCount() if self.mc else None
        self.mcpmts = [self.mc.GetMCPMT(ipmt)
                       for ipmt in range(self.mcpmtcount)] if self.mc else None
        self.mcpmt_types = [
            pmt.GetType() for pmt in self.mcpmts] if self.mc else None
        self.mcpmt_ids = [
            pmt.GetID() for pmt in self.mcpmts] if self.mc else None

    def get_pmt_pos(self, pmt):
        return tuple(self.pmtinfo.GetPosition(pmt.GetID()))

    def get_pmt_type(self, pmt):
        return self.pmtinfo.GetType(pmt.GetID())

    def get_mcphoton_times(self, mcpmt):
        return [mcpmt.GetMCPhoton(i).GetFrontEndTime() for i in range(mcpmt.GetMCPhotonCount())]

    def get_mcphoton_pos(self, mcpmt):
        return [tuple(mcpmt.GetMCPhoton(i).GetPosition()) for i in range(mcpmt.GetMCPhotonCount())]


def gen_event(filepath: str):
    f = root.TFile(filepath, "READ")
    t = f.Get("T")
    rt = f.Get("runT")
    assert rt.GetEntries() == 1
    rt.GetEntry(0)
    pmtinfo = rt.run.GetPMTInfo()
    n = 0
    while n < t.GetEntries():
        #print(f"Yielding entry {n} out of {t.GetEntries()}")
        yield Event(t, pmtinfo, n)
        n += 1


def pc_to_mvns(pc, termination_ohms=50):
    # pc = charge / 1e3 * 1e-9 / 50 / 1e-12
    coulombs = pc * 1e-12
    vs = coulombs * termination_ohms
    mvs = vs * 1000
    mvns = mvs / 1e-9
    return mvns


def find_nearest_MChit(reco_hit, mchits):
    closest = 999999999.9
    index = 999999
    for i, mchit in enumerate(mchits):
        mcx = mchit[0]
        mcy = mchit[1]
        #mcz = mchits[2]
        mct = mchit[2]
        deltaR = np.sqrt((mcx - reco_hit.recox)**2 + (mcy - reco_hit.recoy)**2 + (
            (mct * 180) + ((reco_hit.x * (cfg.NS_PER_SAMPLE / cfg.INTERPFACTOR)) * 180)))
        if deltaR < closest:
            closest = deltaR
            index = i
    if index > 999:
        breakpoint()
    return index


def finder(a, b):
    # Taken from:
    # https://stackoverflow.com/questions/44526121/finding-closest-values-in-two-numpy-arrays
    dup = np.searchsorted(a, b)
    uni = np.unique(dup)
    uni = uni[uni < a.shape[0]]
    ret_b = np.zeros(uni.shape[0])
    for idx, val in enumerate(uni):
        bw = np.argmin(np.abs(a[val]-b[dup == val]))
        tt = dup == val
        ret_b[idx] = np.where(tt == True)[0][bw]
    return np.column_stack((uni, ret_b))


def gen_noise(mean, std, shape):
    samples = np.random.normal(mean, std, size=shape)
    return samples


def add_noise(matrix, rms, baseline=0):
    noise = gen_noise(baseline, rms, matrix.shape)
    return matrix + noise


def gen_arb_wave(position, charge, time):
    ldistance = position[0] - -100
    rdistance = 100 - position[0]
    ltime = ldistance / (cfg.STRIPVELOCITY * cfg.SOL)
    rtime = rdistance / (cfg.STRIPVELOCITY * cfg.SOL)
    leftmatrix = np.zeros((28, nsamples))
    rightmatrix = np.zeros((28, nsamples))
    # xs = np.arange(cfg.NS_PER_SAMPLE, nsamples *
    #                cfg.NS_PER_SAMPLE, cfg.NS_PER_SAMPLE)
    ln = LogNormal()
    # Create the waveform readout for this mcphoton
    # Assume the charge share distribution is gaussian. Variable is distance from strip to hit
    # The cloud radius is ~4mm on striplines, assume this is 2sigma value so 2mm is the PDF sigma
    chargepdf = ChargeSharePDF(cfg.CLOUDRADIUS, 0, 200, 1000)
    #print(f"x: {position[0]}, y: {position[1]}, z: {position[2]}")
    charges = []
    for strip in strip_positions:
        distance = abs(strip - position[1])
        charge_deposited = chargepdf.sample(distance)
        charges.append(charge_deposited)
    scaledcharges = [i*(charge/sum(charges)) for i in charges]
    for i, strip in enumerate(strip_positions):
        Q = pc_to_mvns(scaledcharges[i])
        # waveform parameters
        # use 100 ns as default event trigger time, this is arbitrary unless
        # this is integrated into wider simulation
        lcentre = time + ltime
        rcentre = time + rtime
        # the pulse width (but in log form because of the parameterisation of
        # lognormal function used) - this gives ~1.5 ns wide pulses
        sigma = 0.003
        baseline = 0
        # lwave = ln.npcall(xs, [Q, lcentre, sigma, baseline])
        # rwave = ln.npcall(xs, [Q, rcentre, sigma, baseline])
        lwave = ln.npcall2(xs, Q, medianx-lcentre, medianx, sigma)
        rwave = ln.npcall2(xs, Q, medianx-rcentre, medianx, sigma)
        if np.isnan(lwave).any():
            breakpoint()
        # if np.min(lwave) < -7:
        #     plt.plot(lwave)
        #     plt.show()
        #     breakpoint()
        leftmatrix[i] = lwave
        rightmatrix[i] = rwave
    # print(f"Total modelled charge is: {sum(scaledcharges)}")
    # print(f"Real total charge is: {charge}")
    # plt.plot(strip_positions, scaledcharges, "x")
    # plt.show()
    # plt.imshow(leftmatrix, aspect="auto")
    # plt.show()
    return leftmatrix, rightmatrix


filepath = sys.argv[1]
try:
    n_samples = int(sys.argv[2])
except IndexError:
    n_samples = 0
try:
    mp = True if sys.argv[3].lower() == "mp" else False
except IndexError:
    mp = False
if mp:
    print("Running wavesim in multiprocessing mode")
else:
    print("Running wavesim in default mode")
rng = np.random.default_rng(seed=1)

pmt_times = []
lappd_times = []
trigger_times = []
mc_pos = []
reco_hits = []
x_res = []
y_res = []
t_res = []
pos_diffs = []
scint_times = []
ch_times = []
rem_times = []
n_hits = []
n_photons = []
effs = []
n_scint = 0
n_rem = 0
n_ch = 0


def mp_wavesim(event_no):
    x_res = []
    y_res = []
    t_res = []
    trigger_times = []
    f = root.TFile(filepath, "READ")
    t = f.Get("T")
    rt = f.Get("runT")
    assert rt.GetEntries() == 1
    rt.GetEntry(0)
    pmtinfo = rt.run.GetPMTInfo()
    event = Event(t, pmtinfo, event_no)
    if event.pmts:
        for mcpmt in event.mcpmts:
            if mcpmt.GetType() == 2:
                leftmatrix = np.zeros((28, cfg.NSAMPLES-cfg.NREMOVEDSAMPLES))
                rightmatrix = np.zeros((28, cfg.NSAMPLES-cfg.NREMOVEDSAMPLES))
                leftmatrix = add_noise(leftmatrix, NOISE_RMS)
                rightmatrix = add_noise(rightmatrix, NOISE_RMS)
                photon_count = mcpmt.GetMCPhotonCount()
                if n_samples:
                    if n_samples > photon_count:
                        continue
                    samples = rng.choice(
                        photon_count, size=n_samples, replace=False)
                else:
                    samples = range(photon_count)
                print(f"Photon count: {photon_count}")
                mc_hits = []
                for mcpid in samples:
                    mcphoton = mcpmt.GetMCPhoton(int(mcpid))
                    time = mcphoton.GetFrontEndTime() + TIME_OFFSET
                    if time > 200:
                        continue
                    position = np.asarray(mcphoton.GetPosition())
                    charge = mcphoton.GetCharge()
                    mc_hits.append((position[0], position[1], time))
                    left, right = gen_arb_wave(position, charge, time)
                    leftmatrix += left
                    rightmatrix += right
                    lappd_times.append(time)
                levent = LAPPDEvent.build_from_matrix(
                    leftmatrix, rightmatrix, event_no=event.evnumber)
                if np.max(levent.leftmatrix) < 4.0 or np.max(levent.rightmatrix) < 4:
                    continue
                levent.reconstruct(
                    plot=False, template=cfg.TEMPLATE)  # + gen_noise(0, NOISE_RMS/4, cfg.TEMPLATE.shape))
                print(f"Number of hits reconstructed: {len(levent.hits)}")
                print(levent.hits)
                n_photons.append(len(samples))
                n_hits.append(len(levent.hits))
                effs.append(len(levent.hits) / len(samples))
                reco_matches = []
                for recohit in levent.hits:
                    index = find_nearest_MChit(recohit, mc_hits)
                    reco_matches.append(index)
                    x_res.append(recohit.recox - mc_hits[index][0])
                    y_res.append(recohit.recoy - mc_hits[index][1])
                    t_res.append(
                        (recohit.recot) - mc_hits[index][2])
                this_trigger_times = [x_to_t(hit.x) - np.random.normal(TIME_OFFSET, 0.27)
                                      for hit in levent.hits]
                for trig_time in this_trigger_times:
                    trigger_times.append(trig_time)
    return x_res, y_res, t_res, trigger_times


if mp:
    with Pool(16) as pool:
        success = pool.map(mp_wavesim, range(1000))
        for result in success:
            x_res += result[0]
            y_res += result[1]
            t_res += result[2]
            trigger_times += result[3]
else:
    for event in gen_event(filepath):
        print(event.evnumber)
        # if event.evnumber > 1000:
        #     break
        if event.pmts:
            for mcpmt in event.mcpmts:
                if mcpmt.GetType() == 2:
                    leftmatrix = np.zeros(
                        (28, cfg.NSAMPLES-cfg.NREMOVEDSAMPLES))
                    rightmatrix = np.zeros(
                        (28, cfg.NSAMPLES-cfg.NREMOVEDSAMPLES))
                    leftmatrix = add_noise(leftmatrix, NOISE_RMS)
                    rightmatrix = add_noise(rightmatrix, NOISE_RMS)
                    photon_count = mcpmt.GetMCPhotonCount()
                    if n_samples:
                        if n_samples > photon_count:
                            continue
                        samples = rng.choice(
                            photon_count, size=n_samples, replace=False)
                    else:
                        samples = range(photon_count)
                    print(f"Photon count: {photon_count}")
                    mc_hits = []
                    for mcpid in samples:
                        mcphoton = mcpmt.GetMCPhoton(int(mcpid))
                        time = mcphoton.GetFrontEndTime() + TIME_OFFSET
                        if time > 200:
                            continue
                        position = np.asarray(mcphoton.GetPosition())
                        charge = mcphoton.GetCharge()
                        # print(f"Creator process: {mcphoton.GetProcess()}")
                        if mcphoton.GetProcess() == "Cerenkov":
                            n_ch += 1
                        elif mcphoton.GetProcess() == "Scintillation":
                            n_scint += 1
                        elif mcphoton.GetProcess() == "Reemission":
                            n_rem += 1
                        else:
                            breakpoint()
                        mc_hits.append((position[0], position[1], time))
                        left, right = gen_arb_wave(position, charge, time)
                        leftmatrix += left
                        rightmatrix += right
                        lappd_times.append(time)
                    levent = LAPPDEvent.build_from_matrix(
                        leftmatrix, rightmatrix, event_no=event.evnumber)
                    if np.max(levent.leftmatrix) < 4.0 or np.max(levent.rightmatrix) < 4:
                        continue
                    levent.reconstruct(
                        plot=False, template=cfg.TEMPLATE)  # + gen_noise(0, NOISE_RMS/4, cfg.TEMPLATE.shape))
                    print(f"Number of hits reconstructed: {len(levent.hits)}")
                    print(levent.hits)
                    n_photons.append(len(samples))
                    n_hits.append(len(levent.hits))
                    effs.append(len(levent.hits) / len(samples))
                    reco_matches = []
                    for recohit in levent.hits:
                        index = find_nearest_MChit(recohit, mc_hits)
                        reco_matches.append(index)
                        x_res.append(recohit.recox - mc_hits[index][0])
                        y_res.append(recohit.recoy - mc_hits[index][1])
                        t_res.append(
                            (recohit.recot) - mc_hits[index][2])
                    trigger_times += [x_to_t(hit.x) - np.random.normal(TIME_OFFSET, 0.27)
                                      for hit in levent.hits]

pyhist(trigger_times, fit=False, binwidth=0.2)
xposhist = root.TH1D("xpos", "xpos", 300, -100, 100)
yposhist = root.TH1D("ypos", "ypos", 300, -100, 100)
thist = root.TH1D("t", "t", 6000, -100, 100)
for xval, yval, tval in zip(x_res, y_res, t_res):
    xposhist.Fill(xval)
    yposhist.Fill(yval)
    thist.Fill(tval)
xposhist.Fit("gaus")
yposhist.Fit("gaus")
thist.Fit("gaus")
times = root.TH1D("times", "times", 300, -10, 50)
for value in trigger_times:
    times.Fill(value)
times.Draw()
breakpoint()
print("EFFS: ", effs)
breakpoint()
pmtth1d = root.TH1D("test", "test", 40, -10, 10)
for value in pmt_times:
    pmtth1d.Fill(value)
pmtth1d.Draw()
breakpoint()
lth1d = root.TH1D("test2", "test2", 40, -2, 2)
for value in lappd_times:
    lth1d.Fill(value)
lth1d.Draw()
breakpoint()
