import pdb
import sys
from typing import Generator

import matplotlib.pyplot as plt
import numpy as np

import utils.gimmedatwave as gdw
from pulse import Pulse
from utils.cxxbindings import caenreader
from utils.lappdcfg import config as lcfg

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


class StripPulse():

    # Identified pulse on a single strip

    def __init__(self,
                 leftpulse: Pulse,
                 rightpulse: Pulse,
                 slice_point: int = None,
                 event_no: int = None
                 ) -> None:
        # --- Contains 2x Pulse
        self.left = Pulse(leftpulse)
        self.right = Pulse(rightpulse)
        # ---
        self.slice_point = slice_point
        self.event_no = event_no
        if (self.left.peak is not None) and (self.right.peak is not None):
            self.offset = self._get_offset()
            self.peaks_present = True
        else:
            self.peaks_present = False
            self.offset = None
        if self.offset is not None:
            if self.left.cfpeak is not None and self.right.cfpeak is not None:
                self.cfd_offset = self.left.cfpeak - self.right.cfpeak
            else:
                self.cfd_offset = None
        else:
            self.cfd_offset = None

    def _get_offset(self):
        offset = self.left.peaktime - self.right.peaktime
        return offset if offset < MINDISTANCE else None

    def plot(self):
        plt.plot(self.left.smoothedtimes, self.left.smoothedwave, alpha=0.75)
        plt.plot(self.right.smoothedtimes, self.right.smoothedwave, alpha=0.75)
        plt.show()

# /////////////////////////////////////////////////////////////////////////////


class StripEvent():

    def __init__(self, leftwaveform, rightwaveform, event_no=None, cfg=None, peaks=None) -> None:

        # Raw Event information for each strip

        self.rawleftwaveform = leftwaveform
        self.leftwaveform = np.asarray(leftwaveform.wave)
        self.rawrightwaveform = rightwaveform
        self.rightwaveform = np.asarray(rightwaveform.wave)
        self.event_no = event_no
        self.cfg = cfg  # Or just use global?
        # The peaks for each PAIR of Pulses
        self.leftpeaks, self.rightpeaks = self._find_peaks_custom()
        self.leftpeak_heights, self.rightpeak_heights = self._find_peak_heights()
        if peaks is None:
            self._peaks, self.coarse_offsets = self.coarse_correlate(MAXOFFSET)
            self.triggered = True
        else:
            self._peaks = peaks
            self.coarse_offsets = None
            self.triggered = False
        if self._peaks:
            self.passed = True
        else:
            self.passed = False
        self.pulses = []
        self.offsets = []
        for peak in self._peaks:
            # Should both waves be sliced at the same point? i.e. midpoint between two peaks
            if self.triggered:
                self.search_pulse(
                    int((peak[0] + peak[1]) / 2), LOOKBACK, LOOKFORWARD)
            else:
                self.add_pulse(
                    int((peak[0] + peak[1]) / 2), LOOKBACK, LOOKFORWARD)

    @classmethod
    def build(cls, leftfile: str, rightfile: str, event_no: int
              ) -> "StripEvent":
        leftwave = caenreader.readCAENWave(leftfile, event_no)
        rightwave = caenreader.readCAENWave(rightfile, event_no)
        caenreader.preprocessWave(leftwave)
        caenreader.preprocessWave(rightwave)
        return cls(leftwave, rightwave, event_no=event_no)

    @classmethod
    def itr_file(cls, stripnumber: int, dir: str
                 ) -> Generator["StripEvent", None, None]:
        """Yields all events in a file from a single strip. Finds files from dir."""
        leftchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"L"]
        leftfile = gdw.get_filename(leftchannel, base=dir)
        rightchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"R"]
        rightfile = gdw.get_filename(rightchannel, base=dir)
        leftentries = caenreader.getNumberEntries(leftfile)
        rightentries = caenreader.getNumberEntries(rightfile)
        print(f"Left channel file: {leftfile}")
        print(f"Right channel file: {rightfile}")
        if leftentries != rightentries:
            sys.exit("Each side of strip does not have same number of entries")
        print(f"{leftentries} entries for this strip.")
        for entry in range(leftentries):
            yield cls.build(leftfile, rightfile, entry)

    def _find_peaks_custom(self):
        # print("Finding left peaks")
        leftpeaks = caenreader.findPeaks(
            self.rawleftwaveform.wave, -MINHEIGHT, int(MINDISTANCE / NS_PER_SAMPLE), 5.0)
        # print("Finding right peaks")
        rightpeaks = caenreader.findPeaks(
            self.rawrightwaveform.wave, -MINHEIGHT, int(MINDISTANCE / NS_PER_SAMPLE), 5.0)
        return np.asarray(leftpeaks), np.asarray(rightpeaks)

    def _find_peak_heights(self):
        leftpeak_heights = self.leftwaveform[self.leftpeaks]
        rightpeak_heights = self.rightwaveform[self.rightpeaks]
        return leftpeak_heights, rightpeak_heights

    def add_pulse(self,
                  slice_point: int,
                  lookback: float = LOOKBACK,
                  lookforward: float = LOOKFORWARD
                  ) -> None:
        """Adds a pulse centered around slice_point to the StripEvent"""
        slicedleft = caenreader.sliceAroundPeak(
            self.rawleftwaveform, slice_point, lookback, lookforward)
        slicedright = caenreader.sliceAroundPeak(
            self.rawrightwaveform, slice_point, lookback, lookforward)
        pulsepair = StripPulse(slicedleft, slicedright,
                               slice_point=slice_point, event_no=self.event_no)
        self.pulses.append(pulsepair)

    def search_pulse(self,
                     slice_point: int,
                     lookback: float = LOOKBACK,
                     lookforward: float = LOOKFORWARD
                     ) -> None:
        """Searches for a pulse centered around slice_point to the StripEvent,
        appends it if it finds a peak but ignores if not."""
        slicedleft = caenreader.sliceAroundPeak(
            self.rawleftwaveform, slice_point, lookback, lookforward)
        slicedright = caenreader.sliceAroundPeak(
            self.rawrightwaveform, slice_point, lookback, lookforward)
        pulsepair = StripPulse(slicedleft, slicedright,
                               slice_point=slice_point, event_no=self.event_no)
        if pulsepair.peaks_present:
            self.pulses.append(pulsepair)

    def coarse_correlate(self,
                         max_offset,
                         minrelpulseheight=MINRELPULSEHEIGHT,
                         maxrelpulseheight=MAXRELPULSEHEIGHT,
                         ns_per_sample=NS_PER_SAMPLE):
        # An attempt at correlating peaks on two waveform produced from either side of a stripline
        # Looks for pulses within max_offset (ns) and chooses the one with the largest pulse height
        peak_pairs = []
        offsets = []
        # Loop through left peaks first
        for i, ipeak in enumerate(self.leftpeaks):
            potential_partners = []
            partner_offsets = []
            partner_heights = []
            for j, jpeak in enumerate(self.rightpeaks):
                offset = (ipeak - jpeak) * ns_per_sample
                # check that pulses are within max_offset of each other in time
                if abs(offset) < max_offset:
                    # check that the pulse height is within rel_threshold
                    try:
                        relative_height = self.leftpeak_heights[i] / \
                            self.rightpeak_heights[j]
                    except IndexError:
                        pdb.set_trace()
                    if minrelpulseheight < relative_height < maxrelpulseheight:
                        potential_partners.append(jpeak)
                        partner_offsets.append(offset)
                        partner_heights.append(relative_height)
                        # peak_pairs.append((ipeak, jpeak))
                        # offsets.append(offset)
            try:
                partner_peak = np.argmax(partner_heights)
                peak_pairs.append((ipeak, potential_partners[partner_peak]))
                offsets.append(partner_offsets[partner_peak])
            except ValueError:
                # No matching peaks
                return peak_pairs, offsets
        # save the relative heights as well
        return peak_pairs, offsets

    def plot_raw(self):
        plt.plot(self.rawleftwaveform.times, self.rawleftwaveform.wave)
        plt.plot(self.rawrightwaveform.times, self.rawrightwaveform.wave)
        plt.show()

    def plot_peak(self, peak, smoothed=True):
        if smoothed:
            plt.plot(self.pulses[peak][0].smoothedtimes,
                     self.pulses[peak][0].smoothedwave)
            plt.plot(self.pulses[peak][1].smoothedtimes,
                     self.pulses[peak][1].smoothedwave)
        else:
            plt.plot(self.pulses[peak][0].times, self.pulses[peak][0].wave)
            plt.plot(self.pulses[peak][1].times, self.pulses[peak][1].wave)
        plt.show()