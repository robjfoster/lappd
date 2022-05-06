import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root

from .lappdevent import LAPPDEvent
from .strip import StripEvent
from .utils.lappdcfg import config as lcfg
from .utils.lognormal import double_gauss

daqconfig = lcfg['DAQCONFIG']

NS_PER_SAMPLE = daqconfig.getfloat("nspersample")
NSAMPLES = daqconfig.getint("nsamples")

base_dir = sys.argv[1]
stripnumber = sys.argv[2]


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
sample_start = int(50 / NS_PER_SAMPLE)
sample_end = int(75 / NS_PER_SAMPLE)
int_start = int(192 / NS_PER_SAMPLE)
int_end = int(198 / NS_PER_SAMPLE)

# for sevent in StripEvent.itr_file(stripnumber, base_dir):
#     leftwave = sevent.leftwaveform[sample_start:sample_end]
#     charge = -1 * np.trapz(leftwave, dx=NS_PER_SAMPLE)
#     if np.min(leftwave) > -4.0 or np.max(leftwave) > 5.0:
#         if charge > 10:
#             continue
#     charges.append(charge)

for sevent in StripEvent.itr_file(stripnumber, base_dir):
    if sevent.pulses:
        for pulse in sevent.pulses:
            charge = -1 * np.trapz(pulse.left.wave, dx=NS_PER_SAMPLE)
            pulse_charges.append(charge / 1e3 * 1e-9 / 50 / 1e-12)
    else:
        charge = -1 * np.trapz(sevent.leftwaveform[int_start:int_end])
        other_charges.append(charge / 1e3 * 1e-9 / 50 / 1e-12)

th1 = root.TH1F("hist", "hist", 150, -1, 3)
for value in pulse_charges+other_charges:
    th1.Fill(value)
th1.SetTitle("Gain;pC;Count")
tf1 = double_gauss()
tf1.SetParameters(1000, -0.2, 0.25, 1000, 0.5, 0.15)
th1.Fit("doublegauss")
th1.Draw()
breakpoint()
