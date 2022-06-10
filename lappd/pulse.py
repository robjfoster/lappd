from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from .utils import lappdcfg as cfg
from .utils import sigutils as su
from .utils.lognormal import LogNormal, root_ln


class Pulse():

    def __init__(self, rawpulse) -> None:
        self.rawpulse = rawpulse
        self.ns_per_sample = cfg.NS_PER_SAMPLE
        self.wave = np.asarray(self.rawpulse.wave)
        self.times = np.asarray(self.rawpulse.times)
        self.smoothedtimes, self.smoothedwave = self._interpolate()
        self.rawpeak = self.rawpulse.peakSample
        self.rawinterppeak = self.rawpeak * cfg.INTERPFACTOR
        self.peak = self._getpeak()
        if self.peak:
            self.peak_present = True
            self.peaktime = self.smoothedtimes[self.peak]
            self.height = self.smoothedwave[self.peak]
            self.cfpeak = su.cfd(self.smoothedwave, 0.2,
                                 times=self.smoothedtimes, userpeak=self.peak)
            if self.cfpeak is None:
                self.peak_present = False
        else:
            self.peak_present = False
            self.peaktime = None
            self.height = self.smoothedwave[self.rawinterppeak]
            self.cfpeak = None

    def _interpolate(self, interpfactor: int = cfg.INTERPFACTOR
                     ) -> Tuple[np.ndarray, np.ndarray]:
        x, y = su.cubic_spline(self.times, self.wave, interpfactor)
        return x, y

    def _getpeak(self) -> int:
        peaks, peakinfo = signal.find_peaks(
            self.smoothedwave*-1,
            height=cfg.MINHEIGHT,
            distance=3.0/cfg.NS_PER_SAMPLE * cfg.INTERPFACTOR,
            width=1.0/cfg.NS_PER_SAMPLE * cfg.INTERPFACTOR)
        # TODO: Handle finding multiple peaks in a pulse
        if len(peaks) != 1:
            # self.plot()
            # if len(peaks) > 1:
            #    breakpoint()
            return None
        return peaks[0]

    def fit(self):
        return root_ln(self.times, self.wave)

    def plot(self) -> None:
        plt.plot(self.smoothedtimes, self.smoothedwave)
        plt.plot(self.times, self.wave, "x")
        plt.plot(self.times[self.rawpeak],
                 self.wave[self.rawpeak], "x", c="r")
        plt.show()
