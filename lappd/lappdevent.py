from abc import ABC
import pdb
import sys
from typing import Generator

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from mpl_toolkits import mplot3d
from scipy import signal

import utils.gimmedatwave as gdw
import utils.sigutils as su
from utils.lappdcfg import config as lcfg
from utils.cxxbindings import caenreader

# root.gSystem.Load(os.path.dirname(os.path.realpath(__file__))
#                  + "/utils/caenreader_cpp.so")
#caenreader = root.CAENReader

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


class BaseEvent(ABC):

    @property
    def parent(self) -> "BaseEvent":
        return self._parent

    @parent.setter
    def parent(self, parent: "BaseEvent") -> None:
        self._parent = parent

    def add(self, event: "BaseEvent") -> None:
        pass

    def remove(self, event: "BaseEvent") -> None:
        pass

# /////////////////////////////////////////////////////////////////////////////


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
        else:
            self.peak_present = False
            self.peaktime = None
            self.height = self.smoothedwave[self.rawinterppeak]
            self.cfpeak = None

    def _interpolate(self, interpfactor=INTERPFACTOR):
        x, y = su.cubic_spline(self.times, self.wave, interpfactor)
        return x, y

    def _getpeak(self):
        peaks, peakinfo = signal.find_peaks(
            self.smoothedwave*-1,
            height=MINHEIGHT,
            distance=3.0/NS_PER_SAMPLE * INTERPFACTOR,
            width=1.0/NS_PER_SAMPLE * INTERPFACTOR)
        # TODO: Handle finding multiple peaks in a pulse
        if len(peaks) != 1:
            # self.plot()
            if len(peaks) > 1:
                pdb.set_trace()
            return None
        return peaks[0]

    def plot(self):
        plt.plot(self.smoothedtimes, self.smoothedwave)
        plt.plot(self.times, self.wave, "x")
        plt.plot(self.times[self.rawpeak],
                 self.wave[self.rawpeak], "x", c="r")
        plt.show()

# /////////////////////////////////////////////////////////////////////////////


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

# /////////////////////////////////////////////////////////////////////////////


class LAPPDEvent():

    # Hold raw events for all strips

    def __init__(self, stripevents, event_no=None) -> None:
        self.stripevents = stripevents
        self.event_no = event_no

    @classmethod
    def build(cls, stripfiles, event_no):
        stripevents = {}
        for strip in stripfiles:
            leftfile = stripfiles[strip][0]
            rightfile = stripfiles[strip][1]
            stripevent = StripEvent.build(leftfile, rightfile, event_no)
            stripevents[strip] = stripevent
        return cls(stripevents, event_no=event_no)

    @classmethod
    def search_strip(cls, stripnumber, dir):
        stripfiles = cls.get_stripfiles(dir)
        for stripevent in StripEvent.itr_file(stripnumber, base_dir):
            if stripevent.pulses:
                event_no = stripevent.rawleftwaveform.eventNo
                levent = LAPPDEvent.build(stripfiles, event_no)
                lpulses = []
                for pulse in stripevent.pulses:
                    # breakpoint()
                    lpulses.append(LAPPDPulse(
                        levent, pulse.slice_point, event_no=event_no))
                yield levent, lpulses

    @staticmethod
    def get_stripfiles(dir):
        stripfiles = {}
        for stripnumber in range(14, -15, -1):
            if stripnumber == 0:
                continue
            print(f"Looking for strip {stripnumber}")
            leftchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"L"]
            rightchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"R"]
            try:
                leftfile = gdw.get_filename(leftchannel, base=dir)
                rightfile = gdw.get_filename(rightchannel, base=dir)
            except FileNotFoundError:
                print(f"Did not find files for strip {stripnumber}, skipping.")
                continue
            leftentries = caenreader.getNumberEntries(leftfile)
            rightentries = caenreader.getNumberEntries(rightfile)
            print(f"Left channel file: {leftfile}")
            print(f"Right channel file: {rightfile}")
            if leftentries != rightentries:
                sys.exit("Each side of strip does not have same number of entries")
            print(f"{leftentries} entries for this strip.")
            stripfiles[stripnumber] = (leftfile, rightfile)
        return stripfiles

    @classmethod
    def itr_num(cls, dir, n_events=1):
        stripfiles = {}
        for stripnumber in range(14, -15, -1):
            if stripnumber == 0:
                continue
            print(f"Looking for strip {stripnumber}")
            leftchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"L"]
            rightchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"R"]
            try:
                leftfile = gdw.get_filename(leftchannel, base=dir)
                rightfile = gdw.get_filename(rightchannel, base=dir)
            except FileNotFoundError:
                print(f"Did not find files for strip {stripnumber}, skipping.")
                continue
            leftentries = caenreader.getNumberEntries(leftfile)
            rightentries = caenreader.getNumberEntries(rightfile)
            print(f"Left channel file: {leftfile}")
            print(f"Right channel file: {rightfile}")
            if leftentries != rightentries:
                sys.exit("Each side of strip does not have same number of entries")
            print(f"{leftentries} entries for this strip.")
            stripfiles[stripnumber] = (leftfile, rightfile)
        for event_no in range(n_events):
            yield cls.build(stripfiles, event_no)

    def create_LAPPD_pulse(self, stripnumber, slice_point):
        # thisstripevent = self.stripevents[stripnumber]
        # for pulsepair in thisstripevent.pulses:
        #     for strip in self.stripevents.keys():
        #         if strip == stripnumber:
        pass

    def plot_samples(self):
        @np.vectorize
        def get_sample(strip, sample):
            return self.stripevents[strip].leftwaveform[sample] * -1
        strips = [strip for strip in self.stripevents.keys()]
        samples = [i for i in range(1014)]
        times = [i*NS_PER_SAMPLE for i in samples]
        x, y = np.meshgrid(strips, samples)
        z = get_sample(x, y)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.view_init(elev=25.0, azim=-135)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Stripnumber")
        ax.set_zlabel("Amplitude (mV)")
        #ax.plot_surface(x, y, z, cmap='binary', cstride=8, rstride=1014)
        for strip in strips:
            ax.plot3D(times, [strip for i in range(1014)],
                      get_sample(strip, [i for i in range(1014)]), linestyle="-", marker="o", markersize=1, linewidth=1)
        plt.show()

# /////////////////////////////////////////////////////////////////////////////


class LAPPDPulse():
    """A container for multiple StripPulses corresponding to an entire LAPPD pulse"""

    def __init__(self, lappdevent: LAPPDEvent, slice_point: int, event_no: int = None) -> None:
        self.lappdevent = lappdevent
        self.slice_point = slice_point
        self.event_no = event_no
        self.strippulses = {}
        for strip in self.lappdevent.stripevents:
            self._add_pulse(strip, slice_point)

    @classmethod
    def search(cls, lappdevent: LAPPDEvent, strip=13):
        # Build the whole LAPPDPulse, depending on the search.
        # search = stripnumber or "all" etc
        # Needs to return multiple LAPPDPulses if necessary
        searchevent: StripEvent = lappdevent.stripevents[strip]

    def _add_pulse(self,
                   strip: int,
                   slice_point: int,
                   lookback: float = LOOKBACK,
                   lookforward: float = LOOKFORWARD
                   ) -> None:
        """Adds a StripPulse centered around slice_point to the LAPPDPulse"""
        slicedleft = caenreader.sliceAroundPeak(
            self.lappdevent.stripevents[strip].rawleftwaveform, slice_point, lookback, lookforward)
        slicedright = caenreader.sliceAroundPeak(
            self.lappdevent.stripevents[strip].rawrightwaveform, slice_point, lookback, lookforward)
        strippulse = StripPulse(slicedleft, slicedright,
                                slice_point=slice_point, event_no=self.event_no)
        self.strippulses[strip] = strippulse

    def centroid(self):
        leftvals = []
        rightvals = []
        leftints = []
        rightints = []
        for strip in self.strippulses.keys():
            strippulse = self.strippulses[strip]
            # leftvals.append(
            #    strippulse.left.smoothedwave[strippulse.left.interppeak] * -1.0)
            # rightvals.append(
            #    strippulse.right.smoothedwave[strippulse.right.interppeak] * -1.0)
            leftvals.append(strippulse.left.height)
            rightvals.append(strippulse.right.height)
            leftints.append(caenreader.integratedCharge(
                strippulse.left.rawpulse.wave))
            rightints.append(caenreader.integratedCharge(
                strippulse.right.rawpulse.wave))
        # breakpoint()
        #plt.plot(self.strippulses.keys(), leftvals)
        #plt.plot(self.strippulses.keys(), rightvals)
        # plt.show()
        #plt.plot(self.strippulses.keys(), leftints)
        #plt.plot(self.strippulses.keys(), rightints)
        # plt.show()

    def plot(self):
        @np.vectorize
        def get_sample(strip, sample):
            return self.strippulses[strip].left.smoothedwave[sample] * -1
        strips = [strip for strip in self.strippulses.keys()]
        samples = [i for i in range(len(
            self.strippulses[strips[0]].left.smoothedtimes))]
        times = self.strippulses[strips[0]].left.smoothedtimes
        x, y = np.meshgrid(strips, samples)
        z = get_sample(x, y)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.view_init(elev=25.0, azim=-135)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Stripnumber")
        ax.set_zlabel("Amplitude (mV)")
        #ax.plot_surface(x, y, z, cmap='binary', cstride=8, rstride=1014)
        for strip in strips:
            ax.plot3D(times, [strip for i in range(len(samples))],
                      get_sample(strip, [i for i in range(len(samples))]),
                      linestyle="-",
                      marker="o",
                      markersize=1,
                      linewidth=1,
                      alpha=0.60)
        plt.show()

# /////////////////////////////////////////////////////////////////////////////


if __name__ == "__main__":
    try:
        base_dir = sys.argv[1]
    except IndexError:
        sys.exit("Specify a directory.")
    try:
        stripnumber = sys.argv[2]
    except IndexError:
        sys.exit("Specify a stripnumber")
    lheights = []
    rheights = []
    offsets = []
    for levent, lpulses in LAPPDEvent.search_strip(stripnumber, base_dir):
        lpulses[0].centroid()
