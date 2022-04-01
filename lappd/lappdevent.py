import sys
from abc import ABC

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from mpl_toolkits import mplot3d

import utils.gimmedatwave as gdw
from strip import StripEvent, StripPulse
from utils.cxxbindings import caenreader
from utils.lappdcfg import config as lcfg

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
