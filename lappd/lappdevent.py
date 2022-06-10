import os
import sys
from abc import ABC
from typing import Dict, Generator, Tuple

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from mpl_toolkits import mplot3d

from .matching import detect_peaks, match_peaks
from .strip import StripEvent, StripPulse
from .utils import gimmedatwave as gdw
from .utils import lappdcfg as cfg
from .utils.cxxbindings import caenreader
from .utils.interpolation import interp_matrix
from .utils.lappdcfg import allstrips
from .utils.wiener import do_wiener


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

    def __init__(self, stripevents: Dict[int, StripEvent], event_no: int = None) -> None:
        self.stripevents = stripevents
        self.event_no = event_no
        self.leftmatrix, self.rightmatrix = self._build_matrix()
        self.rootfile = None
        # Available after running reconstruct()
        self.deconvolved = None
        self.interped = None
        self.peaks = None
        self.pairs = None

    @classmethod
    def build(cls, stripfiles: Dict[int, Tuple[str, str]], event_no: int) -> "LAPPDEvent":
        stripevents = {}
        for strip in stripfiles:
            leftfile = stripfiles[strip][0]
            rightfile = stripfiles[strip][1]
            stripevent = StripEvent.build(leftfile, rightfile, event_no)
            stripevents[strip] = stripevent
        return cls(stripevents, event_no=event_no)

    @classmethod
    def search_strip(cls, stripnumber: int, dir: str
                     ) -> Generator[Tuple["LAPPDEvent", "LAPPDPulse"], None, None]:
        stripfiles = cls.get_stripfiles(dir)
        for stripevent in StripEvent.itr_file(stripnumber, dir):
            if stripevent.pulses:
                event_no = stripevent.rawleftwaveform.eventNo
                levent = cls.build(stripfiles, event_no)
                lpulses = []
                for pulse in stripevent.pulses:
                    lpulses.append(LAPPDPulse(
                        levent, pulse.slice_point, event_no=event_no))
                yield levent, lpulses

    @classmethod
    def search_all(cls, dir: str):
        stripfiles = cls.get_stripfiles(dir)
        n_entries = caenreader.getNumberEntries(
            list(stripfiles.values())[0][0])
        print(f"Found {n_entries} events in this directory.")
        for event_no in range(n_entries):
            levent = cls.build(stripfiles, event_no)
            for stripevent in levent.stripevents.values():
                if stripevent.passed:
                    yield levent
                    break

    @classmethod
    def itr_all(cls, dir: str):
        stripfiles = cls.get_stripfiles(dir)
        n_entries = caenreader.getNumberEntries(
            list(stripfiles.values())[0][0])
        print(f"Found {n_entries} events in this directory.")
        for event_no in range(n_entries):
            levent = cls.build(stripfiles, event_no)
            yield levent

    @staticmethod
    def get_stripfiles(dir: str) -> Dict[int, Tuple[str, str]]:
        stripfiles = {}
        for stripnumber in range(14, -15, -1):
            if stripnumber == 0:
                continue
            print(f"Looking for strip {stripnumber}")
            leftchannel = cfg.config['STRIPTODAQ'][str(stripnumber)+"L"]
            rightchannel = cfg.config['STRIPTODAQ'][str(stripnumber)+"R"]
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
    def itr_num(cls, dir: str, n_events: int = 1) -> Generator["LAPPDEvent", None, None]:
        stripfiles = {}
        for stripnumber in range(14, -15, -1):
            if stripnumber == 0:
                continue
            print(f"Looking for strip {stripnumber}")
            leftchannel = cfg.config['STRIPTODAQ'][str(stripnumber)+"L"]
            rightchannel = cfg.config['STRIPTODAQ'][str(stripnumber)+"R"]
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

    def find_lappd_pulses(self, time_delta):
        pulse_centres = []
        pulse_positions = []
        for sevent in self.stripevents.values():
            for spulse in sevent.pulses:
                pulse_centres.append(spulse.slice_point)
                pulse_positions.append(spulse.get_transverse_position())
        time_centres = [i * cfg.NS_PER_SAMPLE for i in pulse_centres]
        breakpoint()

    def _build_matrix(self):
        nsamples = (cfg.NSAMPLES - cfg.NREMOVEDSAMPLES)
        leftmatrix = np.zeros((28, nsamples))
        rightmatrix = np.zeros((28, nsamples))
        for i, strip in enumerate(allstrips):
            try:
                stripevent = self.stripevents[strip]
                leftvalues = stripevent.leftwaveform
                rightvalues = stripevent.rightwaveform
            except KeyError:
                leftvalues = [0 for i in range(nsamples)]
                rightvalues = [0 for i in range(nsamples)]
                pass
            try:
                leftmatrix[i, :] = leftvalues
                rightmatrix[i, :] = rightvalues
            except:
                breakpoint()
        return leftmatrix * -1, rightmatrix * -1

    def reconstruct(self):
        left_deconvolved = do_wiener(self.leftmatrix, template=cfg.TEMPLATE)
        right_deconvolved = do_wiener(self.rightmatrix, template=cfg.TEMPLATE)
        # Change this to dynamically calculate the number of strips to show
        left_interp = interp_matrix(left_deconvolved, startx=0, stopx=cfg.NSAMPLES -
                                    cfg.NREMOVEDSAMPLES, starty=0, stopy=28, interpfactor=10)
        right_interp = interp_matrix(right_deconvolved, startx=0, stopx=cfg.NSAMPLES -
                                     cfg.NREMOVEDSAMPLES, starty=0, stopy=28, interpfactor=10)
        # Change this to also account for minimum distance
        leftpeaks = detect_peaks(left_interp, threshold=cfg.MINHEIGHT)
        rightpeaks = detect_peaks(right_interp, threshold=cfg.MINHEIGHT)
        # Change this to allow for min likelihood (add to config)
        pairs = match_peaks(leftpeaks, rightpeaks)
        self.deconvolved = (left_deconvolved, right_deconvolved)
        self.interped = (left_interp, right_interp)
        self.peaks = (leftpeaks, rightpeaks)
        self.pairs = pairs

    def plot_side(self, side, show=True) -> None:
        # Voltage is inverted!
        if side == "right":
            @ np.vectorize
            def get_sample(strip, sample):
                return self.stripevents[strip].rightwaveform[sample] * -1
        elif side == "left":
            @ np.vectorize
            def get_sample(strip, sample):
                return self.stripevents[strip].leftwaveform[sample] * -1
        else:
            print("Did not recognise side.")
            return
        strips = [strip for strip in self.stripevents.keys()]
        samples = [i for i in range(1014)]
        times = [i*cfg.NS_PER_SAMPLE for i in samples]
        x, y = np.meshgrid(strips, samples)
        # z = get_sample(x, y)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.view_init(elev=25.0, azim=-135)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Stripnumber")
        ax.set_zlabel("Amplitude (mV)")
        # ax.plot_surface(x, y, z, cmap='binary', cstride=8, rstride=1014)
        for strip in strips:
            ax.plot3D(times, [strip for i in range(1014)],
                      get_sample(strip, [i for i in range(1014)]), linestyle="-", marker="o", markersize=1, linewidth=1)
        plt.title(f"{side} side. Event {self.event_no}")
        if show:
            plt.show()

    def plot_both(self) -> None:
        # Voltage is inverted!
        self.plot_side("left", show=False)
        self.plot_side("right", show=False)
        plt.show()

    def plot_average(self):
        # Voltage is inverted!
        @ np.vectorize
        def get_right(strip, sample):
            return self.stripevents[strip].rightwaveform[sample] * -1

        @ np.vectorize
        def get_left(strip, sample):
            return self.stripevents[strip].leftwaveform[sample] * -1
        strips = [strip for strip in self.stripevents.keys()]
        samples = [i for i in range(1014)]
        times = [i*cfg.NS_PER_SAMPLE for i in samples]
        x, y = np.meshgrid(strips, samples)
        # z = get_sample(x, y)
        # fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.view_init(elev=25.0, azim=-135)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Stripnumber")
        ax.set_zlabel("Amplitude (mV)")
        # ax.plot_surface(x, y, z, cmap='binary', cstride=8, rstride=1014)
        for strip in strips:
            ax.plot3D(times, [strip for i in range(1014)],
                      (get_left(strip, [i for i in range(1014)]) +
                       get_right(strip, [i for i in range(1014)])) / 2.0,
                      linestyle="-", marker="o", markersize=1, linewidth=1)
        plt.title(f"LR averaged. Event {self.event_no}")
        plt.show()

    def set_rootfile(self, rootfile):
        self.rootfile = rootfile

    def write_root(self):
        """Not working"""
        breakpoint()
        if self.rootfile is None:
            print("Root file not set, returning")
            return
        if os.path.exists(self.rootfile):
            print("File exists")
            f = root.TFile(self.rootfile, "UPDATE")
        else:
            print("Creating new file")
            f = root.TFile(self.rootfile, "CREATE")
        if not f.GetListOfKeys().Contains("lappd"):
            print("Creating new tree")
            tree = root.TTree("lappd", "LAPPD Event Data")
            tree.Branch("offset", self.stripevents[13].rawleftwaveform.wave)
            tree.Fill()
        else:
            print("Tree exists")
            tree = f.Get("lappd")
            tree.SetBranchAddress(
                "offset", self.stripevents[13].rawleftwaveform.wave)
            tree.Fill()
        tree.Write()
        f.Close()


# /////////////////////////////////////////////////////////////////////////////


class LAPPDPulse():
    """A container for multiple StripPulses corresponding to an entire LAPPD pulse"""

    def __init__(self,
                 lappdevent: LAPPDEvent,
                 slice_point: int,
                 event_no: int = None,
                 lookback: float = cfg.LOOKBACK,
                 lookforward: float = cfg.LOOKFORWARD
                 ) -> None:
        self.lappdevent = lappdevent
        self.slice_point = slice_point
        self.event_no = event_no
        self.lookback = lookback
        self.lookforward = lookforward
        self.strippulses = {}
        for strip in self.lappdevent.stripevents:
            self._add_pulse(
                strip, slice_point, lookback=self.lookback, lookforward=self.lookforward)
        self.leftmatrix, self.rightmatrix = self._build_matrix()

    @ classmethod
    def search(cls, lappdevent: LAPPDEvent, strip: int = 13):
        # Build the whole LAPPDPulse, depending on the search.
        # search = stripnumber or "all" etc
        # Needs to return multiple LAPPDPulses if necessary
        searchevent: StripEvent = lappdevent.stripevents[strip]

    def _add_pulse(self,
                   strip: int,
                   slice_point: int,
                   lookback: float = cfg.LOOKBACK,
                   lookforward: float = cfg.LOOKFORWARD
                   ) -> None:
        """Adds a StripPulse centered around slice_point to the LAPPDPulse"""
        slicedleft = caenreader.sliceAroundPeak(
            self.lappdevent.stripevents[strip].rawleftwaveform, slice_point, lookback, lookforward)
        slicedright = caenreader.sliceAroundPeak(
            self.lappdevent.stripevents[strip].rawrightwaveform, slice_point, lookback, lookforward)
        strippulse = StripPulse(slicedleft, slicedright,
                                slice_point=slice_point, event_no=self.event_no)
        self.strippulses[strip] = strippulse

    def _build_matrix(self):
        lefts = [sp.left.wave for sp in self.strippulses.values()]
        rights = [sp.right.wave for sp in self.strippulses.values()]
        leftmatrix = np.vstack(lefts)
        rightmatrix = np.vstack(rights)
        return leftmatrix, rightmatrix

    def centroid_height(self, plot=False):
        leftvals = []
        rightvals = []
        strips = []
        strippos = []
        for strip in self.strippulses.keys():
            strippulse = self.strippulses[strip]
            leftvals.append(strippulse.left.height)
            rightvals.append(strippulse.right.height)
            strips.append(strip)
            strippos.append(-3.4 + (6.9*strip))
        leftcent = sum(leftvals) / len(leftvals)
        rightcent = sum(rightvals) / len(rightvals)
        gr = root.TGraphErrors(len(leftvals),
                               np.asarray(strippos, dtype="float64"),
                               np.asarray(leftvals, dtype="float64"),
                               np.asarray(
                                   [1.7 for i in range(len(strippos))], dtype="float64"),
                               np.asarray([cfg.VOLTAGEJITTER for i in range(len(strippos))], dtype="float64"))
        tf1 = root.TF1("mygaus", "-1 * gaus(0)", strippos[0], strippos[-1])
        tf1.SetParameter(0, 12)
        tf1.SetParameter(1, strippos[np.argmin(leftvals)])
        tf1.SetParameter(2, 2)
        tf1.SetParLimits(0, min(leftvals)/2.0, min(leftvals) * 3.0)
        tf1.SetParLimits(2, 0.05, 3.5)
        gr.Fit("mygaus")
        gr.SetMarkerSize(3)
        gr.SetTitle("Vertical position; Transverse position; Amplitude (mV)")
        gr.Draw("AP")
        breakpoint()
        if plot:
            plt.errorbar(strippos, leftvals, yerr=cfg.VOLTAGEJITTER, xerr=1.7,
                         marker=".", capsize=5, label="left", c="green")
            plt.plot(strippos, rightvals,
                     "x", label="right", c="purple")
            plt.xlabel = "Strip number"
            plt.ylabel = "Amplitude (mV)"
            plt.xticks(strippos)
            plt.grid(True, which="major", linestyle="--")
            plt.legend()
            plt.show()
        return leftcent, rightcent

    def centroid_integral(self, plot=False):
        leftints = []
        rightints = []
        for strip in self.strippulses.keys():
            strippulse = self.strippulses[strip]
            leftints.append(caenreader.integratedCharge(
                strippulse.left.rawpulse.wave))
            rightints.append(caenreader.integratedCharge(
                strippulse.right.rawpulse.wave))
        leftcent = sum(leftints) / len(leftints)
        rightcent = sum(rightints) / len(rightints)
        if plot:
            plt.plot(self.strippulses.keys(), leftints)
            plt.plot(self.strippulses.keys(), rightints)
            plt.show()
        return leftcent, rightcent

    def plot(self) -> None:
        @ np.vectorize
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
        # ax.plot_surface(x, y, z, cmap='binary', cstride=8, rstride=1014)
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
    cfpeaks = []
    for levent in LAPPDEvent.search_all(base_dir):
        if levent.event_no % 100 == 0:
            print(f"Analysing {levent.event_no}")
        for strip, sevent in levent.stripevents.items():
            for pulse in sevent.pulses:
                print(f"Strip {strip}: Position: {pulse.position}")
        levent.plot_both()
        breakpoint()
    for levent, lpulses in LAPPDEvent.search_strip(stripnumber, base_dir):
        # psf = lpulses[0].leftmatrix[0:3, :]
        # levent.set_rootfile("testfile.root")
        # levent.write_root()
        lpulses[0].centroid_height(plot=True)
        for npulse in lpulses[0].strippulses:
            pulse = lpulses[0].strippulses[npulse].left
            if pulse.peak_present:
                cfpeaks.append(abs(pulse.cfpeak - pulse.peaktime))
    myhist = root.TH1D("hist", "hist", 35, 0.5, 1.6)
    for value in cfpeaks:
        myhist.Fill(value)
    myhist.Draw()
    breakpoint()
