import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root

from .lappdevent import LAPPDEvent
from .strip import StripEvent
from .utils.lappdcfg import config as lcfg

daqconfig = lcfg['DAQCONFIG']

NS_PER_SAMPLE = daqconfig.getfloat("nspersample")
NSAMPLES = daqconfig.getint("nsamples")

base_dir = sys.argv[1]
stripnumber = sys.argv[2]

charges = []
sample_start = int(50 / NS_PER_SAMPLE)
sample_end = int(75 / NS_PER_SAMPLE)
int_start = int(54 / NSAMPLES)
int_end = int(60 / NSAMPLES)

for sevent in StripEvent.itr_file(stripnumber, base_dir):
    leftwave = sevent.leftwaveform[sample_start:sample_end]
    charge = -1 * np.trapz(leftwave, dx=NS_PER_SAMPLE)
    if np.min(leftwave) > -4.0 or np.max(leftwave) > 5.0:
        if charge > 10:
            continue
    charges.append(charge)

th1 = root.TH1F("hist", "hist", 50, -5, 50)
for value in charges:
    th1.Fill(value)
th1.Draw()
breakpoint()
