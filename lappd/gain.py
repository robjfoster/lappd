import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root

from lappd.lappdevent import LAPPDEvent
from lappd.strip import StripEvent
from lappd.utils import lappdcfg as cfg
from lappd.utils import sigutils as su
from lappd.utils.roothist import roothist

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


sample_start = int(80 / cfg.NS_PER_SAMPLE)
sample_end = int(100 / cfg.NS_PER_SAMPLE)
# int_start = int(54 / cfg.NSAMPLES)
# int_end = int(60 / cfg.NSAMPLES)

debug = False


total_pulses = 0
total_skipped = 0
times = []
for sevent in StripEvent.itr_file_raw(stripnumber, base_dir):
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
    if np.min(lwindow) < -8 or np.min(rwindow) < -8:
        total_pulses += 1


print(f"Number of waveforms above threshold: {total_pulses}")
time_values = np.asarray(times)
th1 = root.TH1F("hist", "hist", 400, 0, 200)
for value in time_values:
    th1.Fill(value)
th1.SetTitle("Timing;ns;Count")
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
th1 = root.TH1F("hist", "hist", 250, -5e6, 2e7)
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
