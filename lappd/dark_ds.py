import os
import pdb
import sys
from contextlib import ExitStack

import matplotlib.pyplot as plt
import numpy as np

from utils import gimmedatwave as gdw
from utils import sigutils as su

try:
    base_directory = sys.argv[1]
except IndexError:
    base_directory = "data/20211210-PrelimDarkNoise/"

CAEN_HDTYPE = gdw.def_dtype(np.int32, 6)
CAEN_RDTYPE = gdw.def_dtype(np.float32, 1024)
THRESHOLD = 3.0  # mV
adc_thresh = su.mv_to_adc(THRESHOLD)

dat_files = gdw.find_dats(base_directory)
# Pair together the files for each end of a strip
for pairs in gdw.chunks(dat_files, 2):
    # Open each in a stack for simultaneous analysis
    with ExitStack() as stack:
        files = [stack.enter_context(open(fname, "rb")) for fname in pairs]
        left = files[0]
        right = files[1]
        left.seek(0)
        right.seek(0)
        # analyse left first
        while left.tell() < os.stat(left.name).st_size:
            header, wave = gdw.read_one(left)
            # consistent issue with last few samples
            wave = wave[:-10]
            baseline = su.reduced_mean(wave)
            # subtract baseline
            wave = wave - baseline
            if su.threshold(wave, -adc_thresh):
                plt.plot(range(len(wave)), wave)
                plt.axhline(-adc_thresh, c='purple')
                plt.show()
                pdb.set_trace()
