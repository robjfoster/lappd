import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root

from .lappdevent import LAPPDEvent
from .strip import StripEvent
from .utils import lappdcfg as cfg

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
sample_start = int(50 / cfg.NS_PER_SAMPLE)
sample_end = int(75 / cfg.NS_PER_SAMPLE)
int_start = int(54 / cfg.NSAMPLES)
int_end = int(60 / cfg.NSAMPLES)

for sevent in StripEvent.itr_file(stripnumber, base_dir):
    leftwave = sevent.leftwaveform[sample_start:sample_end]
    charge = -1 * np.trapz(leftwave, dx=cfg.NS_PER_SAMPLE)
    if np.min(leftwave) > -4.0 or np.max(leftwave) > 5.0:
        if charge > 10:
            continue
    charges.append(charge)

th1 = root.TH1F("hist", "hist", 50, -5, 50)
for value in charges:
    th1.Fill(value)
th1.SetTitle("Gain;pC;Count")
tf1 = double_gauss()
tf1.SetParameters(1000, -0.2, 0.25, 1000, 0.5, 0.15)
th1.Fit("doublegauss")
th1.Draw()
breakpoint()
