import os
import pdb
import sys
from contextlib import ExitStack

import numpy as np
import matplotlib.pyplot as plt

from utils import gimmedatwave as gdw
from utils import sigutils as su

base_directory = sys.argv[1]

dat_files = gdw.find_dats(base_directory)


mv_threshold = 2.5
threshold = su.mv_to_adc(mv_threshold)
print(f"Threshold at -{mv_threshold} mV, = {threshold} ADC counts")
count = 0
plot = True

with ExitStack() as stack:
    files = [stack.enter_context(open(fname, "rb")) for fname in dat_files]
    for i, f in enumerate(files):
        f.seek(0)
        while f.tell() < os.stat(f.name).st_size:
            header, wave = gdw.read_one(f)
            wave = wave[:-9]
            samples = range(len(wave))
            reduced_wave = wave[np.abs(wave) < wave + np.std(wave)]
            mean = np.mean(reduced_wave)
            wave = wave - mean
            # and np.std(wave) < 2.8:#np.count_nonzero(wave > 0) < 450:
            if np.any(wave < -threshold):
                count += 1
                header2, wave2 = gdw.find_wave(files[i+1], f.tell()-header[0])
                wave2 = wave2[:-9]
                reduced_wave2 = wave2[np.abs(wave2) < wave2 + np.std(wave2)]
                mean2 = np.mean(reduced_wave2)
                wave2 = wave2 - mean2
                if plot:
                    fig, ax = plt.subplots(2, 1)
                    ax[0].plot(samples, wave)
                    ax[1].plot(samples, wave2)
                    ax[0].axhline(0, c='r', alpha=0.5)
                    ax[1].axhline(0, c='r', alpha=0.5)
                    ax[0].axhline(-threshold, c='purple', alpha=0.5)
                    ax[1].axhline(-threshold, c='purple', alpha=0.5)
                    fig.suptitle(f"File: {f.name.split('/')[-1]}")
                    ax[0].set_title(
                        f"Event byte start {f.tell() - header[0]}. File {f.name.split('/')[-1]}")
                    ax[1].set_title(
                        f"Event byte start {files[i+1].tell() - header[0]}. File {files[i+1].name.split('/')[-1]}")
                    plt.show()
    print(f"Total pulses: {count}")
