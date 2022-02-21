import os
import pdb
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from scipy import signal

from utils.lappdcfg import config as lcfg
import utils.gimmedatwave as gdw
import utils.sigutils as su

root.gSystem.Load(os.path.dirname(os.path.realpath(__file__))
                  + "/utils/caenreader_cpp.so")
caenreader = root.CAENReader


# Attempt at event structure
# Each LAPPDEvent contains 28 StripEvents
# StripEvent contains the raw waveforms for each side of that strip
# Also contains Pulses, sliced waveforms around each identified peak, these
# pulses *should* be correlated

peakparams = lcfg['PEAKPARAMS']
daqconfig = lcfg['DAQCONFIG']

NS_PER_SAMPLE = daqconfig.getfloat("nspersample")
MINHEIGHT = -peakparams.getfloat("minheight")  # mV
MINDISTANCE = peakparams.getfloat("mindistance")
RELHEIGHT = peakparams.getfloat("relheight")
MINRELPULSEHEIGHT = peakparams.getfloat("minrelpulseheight")
MAXRELPULSEHEIGHT = peakparams.getfloat("maxrelpulseheight")
MAXOFFSET = peakparams.getfloat("maxoffset")
INTERPFACTOR = peakparams.getint("interpfactor")
LOOKBACK = peakparams.getfloat("slicelookback")
LOOKFORWARD = peakparams.getfloat("slicelookforward")


class LAPPDEvent():

    def __init__(self) -> None:
        pass

# NEED CONSISTENT NAMING CONVENTION FOR WAVE, WAVEFORM etc


class StripEvent():

    def __init__(self, leftwaveform, rightwaveform, cfg=None) -> None:
        # Pulses should be a list of tuple pairs
        # or should __init__ find all pulses from the left and right waveforms?
        # Waveforms should be preprocessed before creating StripEvent
        self.rawleftwaveform = leftwaveform
        self.leftwaveform = np.asarray(leftwaveform.wave)
        self.rawrightwaveform = rightwaveform
        self.rightwaveform = np.asarray(rightwaveform.wave)
        self.cfg = cfg  # Or just use global?
        # The peaks for each PAIR of Pulses
        self.leftpeaks, self.rightpeaks = self._find_peaks_custom()
        self.leftpeak_heights, self.rightpeak_heights = self._find_peak_heights()
        self.peaks, self.offsets = self.coarse_correlate(MAXOFFSET)
        if self.peaks:
            self.passed = True
        else:
            self.passed = False
        self.pulses = []
        for peak in self.peaks:
            slicedleft = caenreader.sliceAroundPeak(
                self.rawleftwaveform, int(peak[0]), LOOKBACK, LOOKFORWARD)
            slicedright = caenreader.sliceAroundPeak(
                self.rawrightwaveform, int(peak[1]), LOOKBACK, LOOKFORWARD)
            leftpulse = Pulse(slicedleft)
            rightpulse = Pulse(slicedright)
            self.pulses.append((leftpulse, rightpulse))

    def _find_peaks_custom(self):
        leftpeaks = caenreader.findPeaks(
            self.rawleftwaveform.wave, MINHEIGHT, int(MINDISTANCE / NS_PER_SAMPLE), 3.0)
        rightpeaks = caenreader.findPeaks(
            self.rawrightwaveform.wave, MINHEIGHT, int(MINDISTANCE / NS_PER_SAMPLE), 4.0)
        return np.asarray(leftpeaks), np.asarray(rightpeaks)

    def _find_peak_heights(self):
        leftpeak_heights = self.leftwaveform[self.leftpeaks]
        rightpeak_heights = self.rightwaveform[self.rightpeaks]
        return leftpeak_heights, rightpeak_heights

    def coarse_correlate(self,
                         max_offset,
                         minrelpulseheight=MINRELPULSEHEIGHT,
                         maxrelpulseheight=MAXRELPULSEHEIGHT,
                         ns_per_sample=NS_PER_SAMPLE):
        # An attempt at correlating peaks on two waveform produced from either side of a stripline
        # Looks for pulses within max_offset (ns) and chooses the one with the largest pulse height
        # Probably some edge cases that have been missed, but scipy.find_peaks ensures that there
        # should not be multiple identified peaks within a given distance anyway
        peak_pairs = []
        offsets = []
        # Loop through the shorter of the two lists
        for i, ipeak in enumerate(self.leftpeaks):
            potential_partners = []
            partner_offsets = []
            partner_heights = []
            for j, jpeak in enumerate(self.rightpeaks):
                offset = (ipeak - jpeak) * ns_per_sample
                # if first_peaks is rightpeaks then ipeak is right end of strip
                # we want to do left - right
                #offset = offset * -1 if first_peaks is self.rightpeaks else offset
                # check that pulses are within max_offset of each other in time
                if abs(offset) < max_offset:
                    # check that the pulse height is within rel_threshold
                    try:
                        relative_height = self.leftpeak_heights[i] / \
                            self.rightpeak_heights[j]
                    except IndexError:
                        breakpoint()
                    if minrelpulseheight < relative_height < maxrelpulseheight:
                        potential_partners.append(jpeak)
                        partner_offsets.append(offset)
                        partner_heights.append(relative_height)
                        #peak_pairs.append((ipeak, jpeak))
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


class Pulse():
    # Probably shouldn't interact with Pulse directly, only through StripEvent

    def __init__(self, wave, ns_per_sample=NS_PER_SAMPLE) -> None:
        self._rawwave = wave.wave  # c++ vector<float> object
        self.wave = np.asarray(wave.wave)
        self._rawtimes = wave.times  # c++ vector<float> object
        self.times = np.asarray(wave.times)
        self.ns_per_sample = ns_per_sample
        self.smoothedtimes, self.smoothedwave = self._interpolate()
        #self.peak = np.argmin(self.smoothedwave) if not peak else peak
        #self.peak = int(len(self.smoothedwave) / 2.0)
        self.peak = self._getpeak()
        if self.peak:
            self.peaktime = self.smoothedtimes[self.peak]
            self.height = self.smoothedwave[self.peak]
        else:
            self.peaktime = None
            self.height = None

    def _interpolate(self, interpfactor=INTERPFACTOR):
        x, y = su.cubic_spline(self.times, self.wave, interpfactor)
        return x, y

    def _getpeak(self):
        # if len(self.smoothedwave) == (LOOKBACK+LOOKFORWARD)/NS_PER_SAMPLE*INTERPFACTOR:
        #     peak = int(len(self.smoothedwave) / 2)
        # else:
        #     peak = np.argmin(self.smoothedwave)
        # peak = signal.argrelmin(
        #     self.smoothedwave, order=int(LOOKBACK/NS_PER_SAMPLE*INTERPFACTOR))
        # try:
        #     return peak[0][0]
        # except IndexError:
        peaks, peakinfo = signal.find_peaks(
            self.smoothedwave*-1, height=4, distance=3.0/NS_PER_SAMPLE * INTERPFACTOR, width=1.0/NS_PER_SAMPLE * INTERPFACTOR)
        if len(peaks) != 1:
            self.plot()
            return None
        return peaks[0]

    def plot(self):
        plt.plot(self.smoothedtimes, self.smoothedwave)
        plt.plot(self.times, self.wave, "x")
        # plt.plot(self.smoothedtimes[self.peak],
        #          self.smoothedwave[self.peak], "x", c="r")
        plt.show()

# "An event" = two pulses with calculated offset and amplitude


def build_strip_event(stripnumber, event_no, dir):
    leftchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"L"]
    leftfile = gdw.get_filename(leftchannel, base=dir)
    leftwave = caenreader.readCAENWave(leftfile, event_no)
    rightchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"R"]
    rightfile = gdw.get_filename(rightchannel, base=dir)
    rightwave = caenreader.readCAENWave(rightfile, event_no)
    caenreader.preprocessWave(leftwave)
    caenreader.preprocessWave(rightwave)
    event = StripEvent(leftwave, rightwave)
    return event


def get_strip_events(stripnumber, dir):
    leftchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"L"]
    leftfile = gdw.get_filename(leftchannel, base=dir)
    rightchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"R"]
    rightfile = gdw.get_filename(rightchannel, base=dir)
    leftentries = caenreader.getNumberEntries(leftfile)
    rightentries = caenreader.getNumberEntries(rightfile)
    if leftentries != rightentries:
        sys.exit("Each side of strip does not have same number of entries")
    print(f"{leftentries} entries for this strip.")
    for entry in range(leftentries):
        yield build_strip_event(stripnumber, entry, dir)


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
    for event in get_strip_events(7, base_dir):
        if event.pulses:
            for pulse in event.pulses:
                if pulse[0].height is not None:
                    lheights.append(pulse[0].height)
    th1d = root.TH1D("test", "test", 30, 0, 50)
    for value in lheights:
        try:
            th1d.Fill(-value)
        except TypeError:
            breakpoint()
    th1d.Draw()
    breakpoint()
