import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from lappd.matching import strip_positions
from lappd.utils import lappdcfg as cfg
from lappd.utils.lognormal import LogNormal
from lappd.utils.pdf import ChargeSharePDF
from lappd.lappdevent import LAPPDEvent

root.gSystem.Load("libRATEvent.so")
strip_positions = strip_positions - 100

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


def pe_to_mvns(pe, termination_ohms=50):
    gain = 1e6
    coulombs = pe * 1.6e-19 * gain
    vs = coulombs * termination_ohms
    mvs = vs * 1000
    mvns = mvs / 1e-9
    return mvns


def gen_wave(mcphoton):
    position = mcphoton.GetPosition()
    charge = mcphoton.GetCharge()
    time = mcphoton.GetFrontEndTime()  # or GetHitTime() ?
    ldistance = position.x() - -100
    rdistance = 100 - position.x()
    ltime = ldistance / (cfg.STRIPVELOCITY * cfg.SOL)
    rtime = rdistance / (cfg.STRIPVELOCITY * cfg.SOL)
    nsamples = cfg.NSAMPLES - cfg.NREMOVEDSAMPLES
    leftmatrix = np.zeros((28, nsamples))
    rightmatrix = np.zeros((28, nsamples))
    xs = np.arange(cfg.NS_PER_SAMPLE, nsamples *
                   cfg.NS_PER_SAMPLE, cfg.NS_PER_SAMPLE)
    ln = LogNormal()
    # Create the waveform readout for this mcphoton
    # Assume the charge distribution is gaussian. Variable is distance from strip to hit
    chargepdf = ChargeSharePDF(3, 0, 200, 1000)
    print(f"x: {position.x()}, y: {position.y()}, z: {position.z()}")
    charges = []
    for strip in strip_positions:
        distance = abs(strip - position.y())
        charge_deposited = chargepdf.sample(distance)
        charges.append(charge_deposited)
    scaledcharges = [i*(charge/sum(charges)) for i in charges]
    for i, strip in enumerate(strip_positions):
        Q = pe_to_mvns(scaledcharges[i])
        # waveform parameters
        lcentre = 100 + time + ltime
        rcentre = 100 + time + rtime
        sigma = 0.005
        baseline = 0
        lwave = ln.npcall(xs, [Q, lcentre, sigma, baseline])
        rwave = ln.npcall(xs, [Q, rcentre, sigma, baseline])
        leftmatrix[i] = lwave
        rightmatrix[i] = rwave
    print(f"Total modelled charge is: {sum(scaledcharges)}")
    print(f"Real total charge is: {charge}")
    # plt.plot(strip_positions, scaledcharges, "x")
    # plt.show()
    # plt.imshow(leftmatrix, aspect="auto")
    # plt.show()
    return leftmatrix, rightmatrix


filepath = sys.argv[1]
# f = root.TFile(filepath, "READ")
# t = f.Get("T")
# rt = f.Get("runT")
# assert rt.GetEntries() == 1
# rt.GetEntry(0)
# PMTINFO = rt.run.GetPMTInfo()
pmt_times = []
lappd_times = []
for event in gen_event(filepath):
    print(f"PMT IDs: {event.pmt_ids}. MCPMT IDs: {event.mcpmt_ids}")
    print(f"PMT types: {event.pmt_types}")
    if event.pmts:
        for mcpmt in event.mcpmts:
            if mcpmt.GetType() == 2:
                print(f"Photon count: {mcpmt.GetMCPhotonCount()}")
                leftmatrix = np.zeros((28, 1014))
                rightmatrix = np.zeros((28, 1014))
                for mcpid in range(mcpmt.GetMCPhotonCount()):
                    mcphoton = mcpmt.GetMCPhoton(mcpid)
                    time = mcphoton.GetFrontEndTime()
                    left, right = gen_wave(mcphoton)
                    leftmatrix += left
                    rightmatrix += right
                    lappd_times.append(time)
                levent = LAPPDEvent.build_from_matrix(leftmatrix, rightmatrix)
                levent.plot_both()
                breakpoint()
pmtth1d = root.TH1D("test", "test", 40, -10, 10)
#delt = list(np.array(delt) + 5e-9)
for value in pmt_times:
    pmtth1d.Fill(value)
pmtth1d.Draw()
breakpoint()
lth1d = root.TH1D("test2", "test2", 40, -2, 2)
#delt = list(np.array(delt) + 5e-9)
for value in lappd_times:
    lth1d.Fill(value)
lth1d.Draw()
breakpoint()
