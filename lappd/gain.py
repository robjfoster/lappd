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
sample_start = int(75 / cfg.NS_PER_SAMPLE)
sample_end = int(100 / cfg.NS_PER_SAMPLE)
int_start = int(54 / cfg.NSAMPLES)
int_end = int(60 / cfg.NSAMPLES)


# for sevent in StripEvent.itr_file(stripnumber, base_dir):
#     if sevent.pulses:
#         for pulse in sevent.pulses:
#             charge = -1 * np.trapz(pulse.left.wave, dx=NS_PER_SAMPLE)
#             pulse_charges.append(charge / 1e3 * 1e-9 / 50 / 1e-12)
#     else:
#         charge = -1 * np.trapz(sevent.leftwaveform[int_start:int_end])
#         other_charges.append(charge / 1e3 * 1e-9 / 50 / 1e-12)

# for levent in LAPPDEvent.itr_all(base_dir):
#     sevent = levent.stripevents[stripnumber]
#     if sevent.pulses:
#         for pulse in sevent.pulses:
#             charge = -1 * np.trapz(pulse.left.wave, dx=cfg.NS_PER_SAMPLE)
#             pulse_charges.append(charge / 1e3 * 1e-9 / 50 / 1e-12)
#     else:
#         charge = -1 * np.trapz(sevent.leftwaveform[int_start:int_end])
#         other_charges.append(charge / 1e3 * 1e-9 / 50 / 1e-12)

# for levent in LAPPDEvent.itr_all_raw(base_dir):
#     sevent = levent.stripevents[stripnumber]
#     totalcharge = 0
#     totalcharge += -1 * \
#         np.trapz(
#             sevent.leftwaveform[sample_start:sample_end], dx=cfg.NS_PER_SAMPLE)
#     totalcharge += -1 * \
#         np.trapz(
#             sevent.rightwaveform[sample_start:sample_end], dx=cfg.NS_PER_SAMPLE)
#     pulse_charges.append(totalcharge / 1e3 * 1e-9 / 50 / 1e-12)


total_pulses = 0
total_skipped = 0
times = []
for sevent in StripEvent.itr_file_raw(stripnumber, base_dir):
    if su.threshold(sevent.leftwaveform, 3, polarity="positive") or su.threshold(sevent.rightwaveform, 3, polarity="positive"):
        total_skipped += 1
        continue
    if su.threshold(sevent.leftwaveform, -4) and su.threshold(sevent.rightwaveform, -4):
        total_pulses += 1
        # sevent.plot_raw()
        times.append(np.argmin(sevent.leftwaveform) * cfg.NS_PER_SAMPLE)
    totalcharge = 0
    totalcharge += np.trapz(
        sevent.leftwaveform[sample_start:sample_end], dx=cfg.NS_PER_SAMPLE)
    totalcharge += np.trapz(
        sevent.rightwaveform[sample_start:sample_end], dx=cfg.NS_PER_SAMPLE)
    totalcharge = -1 * totalcharge / 1e3 * 1e-9 / 50 / 1e-12
    # if totalcharge < -0.3:
    #     plt.plot(sevent.leftwaveform[sample_start:sample_end])
    #     plt.show()
    #     plt.plot(sevent.rightwaveform[sample_start:sample_end])
    #     plt.show()
    pulse_charges.append(totalcharge)

pulse_charges = np.asarray(pulse_charges)
th1 = root.TH1F("hist", "hist", 250, -0.5, 1.5)
for value in pulse_charges:
    th1.Fill(value)
th1.SetTitle("Gain;pC;Count")
# tf1 = double_gauss()
# tf1.SetParameters(1000, -0.2, 0.25, 1000, 0.5, 0.15)
# th1.Fit("doublegauss")
print(f"################ Total pulses: {total_pulses} ####################")
print(f"Total skipped: {total_skipped}")
th1.Draw()
breakpoint()
