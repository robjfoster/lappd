from logging import root
from typing import Tuple
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from .utils import sigutils as su
from .utils.lappdcfg import config as lcfg
from .utils.lognormal import root_ln

peakparams = lcfg['PEAKPARAMS']
daqconfig = lcfg['DAQCONFIG']

NS_PER_SAMPLE = daqconfig.getfloat("nspersample")
MINHEIGHT = peakparams.getfloat("minheight")  # mV
MINDISTANCE = peakparams.getfloat("mindistance")
RELHEIGHT = peakparams.getfloat("relheight")
MINRELPULSEHEIGHT = peakparams.getfloat("minrelpulseheight")
MAXRELPULSEHEIGHT = peakparams.getfloat("maxrelpulseheight")
MAXOFFSET = peakparams.getfloat("maxoffset")
INTERPFACTOR = peakparams.getint("interpfactor")
LOOKBACK = peakparams.getfloat("slicelookback")
LOOKFORWARD = peakparams.getfloat("slicelookforward")


class Pulse():

    def __init__(self, rawpulse) -> None:
        self.rawpulse = rawpulse
        self.ns_per_sample = NS_PER_SAMPLE
        self.wave = np.asarray(self.rawpulse.wave)
        self.times = np.asarray(self.rawpulse.times)
        self.smoothedtimes, self.smoothedwave = self._interpolate()
        self.rawpeak = self.rawpulse.peakSample
        self.rawinterppeak = self.rawpeak * INTERPFACTOR
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

    def _interpolate(self, interpfactor: int = INTERPFACTOR
                     ) -> Tuple[np.ndarray, np.ndarray]:
        x, y = su.cubic_spline(self.times, self.wave, interpfactor)
        return x, y

    def _getpeak(self) -> int:
        peaks, peakinfo = signal.find_peaks(
            self.smoothedwave*-1,
            height=MINHEIGHT,
            distance=3.0/NS_PER_SAMPLE * INTERPFACTOR,
            width=1.0/NS_PER_SAMPLE * INTERPFACTOR)
        # TODO: Handle finding multiple peaks in a pulse
        if len(peaks) != 1:
            # self.plot()
            # if len(peaks) > 1:
            #    breakpoint()
            return None
        return peaks[0]

    def fit(self):
        root_ln(self.times, self.wave)

    def plot(self) -> None:
        plt.plot(self.smoothedtimes, self.smoothedwave)
        plt.plot(self.times, self.wave, "x")
        plt.plot(self.times[self.rawpeak],
                 self.wave[self.rawpeak], "x", c="r")
        plt.show()
