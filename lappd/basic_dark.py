import os
import pdb
import sys
import time
from contextlib import ExitStack

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root

import lappd.utils.gimmedatwave as gdw
import lappd.utils.sigutils as su
from lappd.utils.cxxbindings import caenreader

try:
    base_directory = sys.argv[1]
except IndexError:
    sys.exit("Specify the directory")


LOWTHRESHOLD = -5.0  # mV
HIGHTHRESHOLD = 3.0  # mV
RISETHRESHOLD = 2.0  # ns

sample_times = gdw.ns(0.2, 1014)

files = gdw.find_dats(base_directory)
for file in files:
    print(f"Reading {file}")
    n_waves = caenreader.getNumberEntries(file)
    print(f"Number of waves: {n_waves}")
    breakpoint()
    total_time = sample_times[-1] * n_waves * 1e-9
    rt_waves = []
    start = time.time()
    coarse_output = caenreader.coarseDarkSearch(
        file, LOWTHRESHOLD, HIGHTHRESHOLD)
    end = time.time()
    print(f"Processing time: {end - start}")
    for wave in coarse_output.waves:
        peaks = caenreader.findPeaks(wave.wave, -4, 5, 2)
        if len(peaks) > 1:
            pdb.set_trace()
        try:
            if su.rise_time(np.asarray(wave.wave), RISETHRESHOLD, fraction1=0.5, fraction2=1.0):
                breakpoint()
                rt_waves.append(wave)
        except IndexError:
            pdb.set_trace()
            pass
    dark_rate = len(rt_waves) / total_time
    dark_rate_noise_rej = len(
        rt_waves) / (sample_times[-1] * (n_waves - coarse_output.rejectedMax) * 1e-9)
    print(f"Number of waves passed threshold: {len(rt_waves)}")
    print(f"Dark rate = {dark_rate} Hz = {dark_rate / 13.5} Hz/cm2")
    print(
        f"Dark rate (noise rejected) = {dark_rate_noise_rej} Hz = {dark_rate_noise_rej / 13.5} Hz/cm2")
    # for caenwave in rt_waves:
    #    np.save("data/examplewaves/" + str(caenwave.eventNo) + ".npy", caenwave.wave)
    # pdb.set_trace()
    # with open("results_250.txt", "a") as f:
    #    f.write(str(dark_rate_noise_rej / 13.5))
    #    f.write("\n")
