import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root

from lappd.lappdevent import LAPPDEvent
from lappd.strip import StripEvent
from lappd.utils import gimmedatwave as gdw
from lappd.utils import lappdcfg as cfg
from lappd.utils import sigutils as su
from lappd.utils.roothist import roothist
from lappd.utils.cxxbindings import caenreader

base_dir = sys.argv[1]
stripnumber = int(sys.argv[2])


def peak_to_gain(peak_val, termination_ohms=50):
    """Assume histogram is in mV*ns"""
    elementary_charge = 1.6e-19
    coulombs = (peak_val / 1e3 * 1e-9) / termination_ohms
    gain = coulombs / elementary_charge
    return gain


def pc_to_gain(pc):
    elementary_charge = 1.6e-19
    gain = (pc * 1e-12) / elementary_charge
    return gain


pulse_charges = []
other_charges = []
charges = []
heights = []
id_pulse_heights = []
ldelta_t = []
rdelta_t = []


sample_start = int(80 / cfg.NS_PER_SAMPLE)
sample_end = int(100 / cfg.NS_PER_SAMPLE)
# int_start = int(54 / cfg.NSAMPLES)
# int_end = int(60 / cfg.NSAMPLES)

sample_times = np.arange(sample_start*cfg.NS_PER_SAMPLE,
                         sample_end*cfg.NS_PER_SAMPLE, cfg.NS_PER_SAMPLE)

debug = False


total_pulses = 0
total_skipped = 0
times = []
trig_times = []

trig_filename = gdw.get_filename(0, base_dir, prefix="TR_0_")

for sevent in StripEvent.itr_file_raw(stripnumber, base_dir):
    trig_caen_wave = caenreader.readOne(trig_filename, sevent.event_no)
    caenreader.preprocessWave(trig_caen_wave)
    trig_wave = np.asarray(trig_caen_wave)
    if np.max(trig_wave) - np.min(trig_wave) < 200:
        continue
    trig_time = np.argmax(trig_wave) * cfg.NS_PER_SAMPLE
    lwave = sevent.leftwaveform
    rwave = sevent.rightwaveform
    # if (np.max(lwave) > 10) or (np.max(rwave) > 10):
    #     print("skipping:", sevent.event_no)
    #     continue
    # plt.plot(lwave)
    # plt.plot(rwave)
    # plt.show()
    lwindow = lwave[sample_start:sample_end]
    rwindow = rwave[sample_start:sample_end]
    lpeak_time = (np.argmax(lwave < -8))
    rpeak_time = np.argmax(rwave < -8)
    peak_time = (lpeak_time+rpeak_time) / 2 * 0.2
    if peak_time > 1:
        if lpeak_time < 1:
            peak_time = rpeak_time
        if rpeak_time < 1:
            peak_time = lpeak_time
        times.append(peak_time)
    inwindowpeak = int((np.argmin(lwindow) + np.argmin(rwindow))/2)
    start_int = inwindowpeak - 10 if inwindowpeak > 10 else 0
    end_int = inwindowpeak + \
        20 if (inwindowpeak + 20) < (sample_end -
                                     sample_start) else (sample_end-sample_start) if start_int != 0 else 30
    if end_int == (sample_end-sample_start):
        start_int = end_int - 30
    totalcharge = 0
    leftcharge = np.trapz(
        lwindow[start_int:end_int], dx=cfg.NS_PER_SAMPLE) * -1
    rightcharge = np.trapz(
        rwindow[start_int: end_int], dx=cfg.NS_PER_SAMPLE) * -1
    totalcharge = leftcharge + rightcharge
    lheight = -1*np.min(lwindow)
    rheight = -1*np.min(rwindow)
    # if abs(lheight - rheight) > 3:
    #     print("Height mismatch: ", sevent.event_no)
    #     continue
    height = (lheight+rheight) / 2
    heights.append(height)
    picocoulombs = totalcharge/1e3*1e-9/50/1.6e-19
    pulse_charges.append(picocoulombs)
    if picocoulombs > 5e7:
        debug = True
    if debug:
        print("Height: ", height)
        print(f"Gain: {picocoulombs:2E}")
        print("start int: ", start_int)
        print("end int: ", end_int)
        plt.plot(lwave)
        plt.plot(rwave)
        plt.axvline(lpeak_time, c='blue', linestyle='dashed')
        plt.axvline(rpeak_time, c='orange', linestyle='dashed')
        plt.axvline(sample_start+start_int, c='purple')
        plt.axvline(sample_start+end_int, c='purple')
        plt.show()
    debug = False
    if np.min(lwindow) < -8 or np.min(rwindow) < -8:
        total_pulses += 1
        # plt.plot(lx, ly)
        # plt.show()
        trig_time = np.argmin(abs(
            trig_wave - ((np.max(trig_wave) + np.min(trig_wave))*0.5))) * cfg.NS_PER_SAMPLE
        lx, ly = su.cubic_spline(sample_times, lwindow, ratio=10)
        rx, ry = su.cubic_spline(sample_times, rwindow, ratio=10)
        tx, ty = su.cubic_spline(cfg.TIMES[int(trig_time/cfg.NS_PER_SAMPLE)-30:int(trig_time/cfg.NS_PER_SAMPLE)+30],
                                 trig_wave[int(trig_time/cfg.NS_PER_SAMPLE)-30:int(trig_time/cfg.NS_PER_SAMPLE)+30], ratio=10)
        # The -ty+ty[0] inverts the waveform for the CFD and performs rudimentary baseline subtraction
        fine_trig_time = su.cfd(-ty+ty[0], 0.5, times=tx)[0]
        lcrossing_time = su.cfd(ly, 0.8, times=lx)[0]
        rcrossing_time = su.cfd(ry, 0.8, times=rx)[0]
        ldelta_t.append(lcrossing_time - fine_trig_time)
        rdelta_t.append(rcrossing_time - fine_trig_time)
        trig_times.append(trig_time)
        # if lcrossing_time - fine_trig_time < 45.5:
        #     plt.plot(lwave)
        #     plt.plot(trig_wave)
        #     plt.show()


print(f"Number of waveforms above threshold: {total_pulses}")
final_vals = []
for lvalue, rvalue in zip(ldelta_t, rdelta_t):
    final_vals.append((lvalue+rvalue)/2.0)
np.save(
    f'./data/tts/{base_dir.split("/")[-2]}_tts.npy', np.asarray(final_vals))
sys.exit()

th1 = root.TH1F("hist", "hist", 300, 38, 44)
for lvalue, rvalue in zip(np.asarray(ldelta_t), np.asarray(rdelta_t)):
    th1.Fill((lvalue + rvalue) / 2.0)
th1.SetTitle("Average TTS;ns;Count")
th1.Draw()
breakpoint()
del th1

tts_values = np.asarray(ldelta_t)
th1 = root.TH1F("hist", "hist", 300, 38, 44)
for value in tts_values:
    th1.Fill(value)
th1.SetTitle("Left TTS;ns;Count")
th1.Draw()
breakpoint()
del th1

tts_values = np.asarray(rdelta_t)
th1 = root.TH1F("hist", "hist", 300, 38, 44)
for value in tts_values:
    th1.Fill(value)
th1.SetTitle("Right TTS;ns;Count")
th1.Draw()
breakpoint()
del th1

time_values = np.asarray(times)
th1 = root.TH1F("hist", "hist", 400, 0, 200)
for value in time_values:
    th1.Fill(value)
th1.SetTitle("Pulse timing;ns;Count")
th1.Draw()
print("Press c to continue")
breakpoint()
del th1

trig_time_values = np.asarray(trig_times)
th1 = root.TH1F("hist", "hist", 400, 0, 200)
for value in trig_time_values:
    th1.Fill(value)
th1.SetTitle("Trigger timing;ns;Count")
th1.Draw()
print("Press c to continue")
breakpoint()
del th1

heights = np.asarray(heights)
th1 = root.TH1F("hist", "hist", 300, 0, 50)
for value in heights:
    th1.Fill(value)
th1.SetTitle("Height;mV;Count")
th1.Draw()
breakpoint()
del th1

pulse_charges = np.asarray(pulse_charges)
th1 = root.TH1F("hist", "hist", 250, -5e6, 2.5e7)
for value in pulse_charges:
    th1.Fill(value)
th1.SetTitle("Gain;Gain;Count")
print(f"################ Total pulses: {total_pulses} ####################")
print(f"Total skipped: {total_skipped}")
th1.Draw()
breakpoint()
savefilename = input("enter file name")
np.save(f"./data/gain/{savefilename}_gain.npy", pulse_charges)
np.save(f"./data/gain/{savefilename}_heights.npy", heights)
print(np.mean(heights[heights > 4]))
print(len(heights[heights > 4]))
