import os
import sys
from abc import ABC
from typing import Dict, Generator, Tuple

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from matplotlib import cm
from mpl_toolkits import mplot3d

from lappd.matching import (Hit, RecoHit, cfd_timestamp, detect_peaks, get_centroid, match_peaks,
                            transverse_position, x_to_t, y_to_loc)
from lappd.strip import StripEvent, StripPulse
from lappd.utils import gimmedatwave as gdw
from lappd.utils import lappdcfg as cfg
from lappd.utils.cxxbindings import caenreader
from lappd.utils.interpolation import interp_matrix
from lappd.utils.lappdcfg import allstrips
from lappd.utils.roothist import roothist
from lappd.utils.wiener import do_wiener


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

    def __init__(self,
                 stripevents: Dict[int, StripEvent],
                 event_no: int = None,
                 trigevent=None
                 ) -> None:
        self.stripevents = stripevents
        self.event_no = event_no
        self.leftmatrix, self.rightmatrix = self._build_matrix()
        self.trigevent = trigevent
        self.rootfile = None
        # Available after running reconstruct()
        self.deconvolved = None
        self.interped = None
        self.peaks = None
        self.pairs = None
        self.hits = None
        self.hiterrors = None

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
    def build_raw(cls,
                  stripfiles: Dict[int, Tuple[str, str]],
                  event_no: int,
                  trigger=None
                  ) -> "LAPPDEvent":
        """build_raw builds StripEvents without performing immediate analysis on waveforms"""
        stripevents = {}
        trigevent = None
        for strip in stripfiles:
            leftfile = stripfiles[strip][0]
            rightfile = stripfiles[strip][1]
            stripevent = StripEvent.build_raw(leftfile, rightfile, event_no)
            stripevents[strip] = stripevent
        if trigger:
            trigevent = StripEvent.build(trigger, trigger, event_no)
        return cls(stripevents, event_no=event_no, trigevent=trigevent)

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

    @classmethod
    def itr_all_raw(cls, dir: str, trigger=None):
        stripfiles = cls.get_stripfiles(dir)
        trig_file = None
        if trigger:
            trig_file = gdw.get_filename(trigger, base=dir, prefix="TR_0")
            n_entries = caenreader.getNumberEntries(trig_file)
        else:
            n_entries = caenreader.getNumberEntries(
                list(stripfiles.values())[0][0])
        print(f"Found {n_entries} events in this directory.")
        for event_no in range(n_entries):
            levent = cls.build_raw(stripfiles, event_no, trigger=trig_file)
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

    def reconstruct(self, plot=False):
        left_deconvolved = do_wiener(self.leftmatrix, cfg.TEMPLATE)
        right_deconvolved = do_wiener(self.rightmatrix, cfg.TEMPLATE)
        # Change this to dynamically calculate the number of strips to show
        left_interp = interp_matrix(left_deconvolved, startx=0, stopx=cfg.NSAMPLES -
                                    cfg.NREMOVEDSAMPLES, starty=0, stopy=28, interpfactor=cfg.INTERPFACTOR)
        right_interp = interp_matrix(right_deconvolved, startx=0, stopx=cfg.NSAMPLES -
                                     cfg.NREMOVEDSAMPLES, starty=0, stopy=28, interpfactor=cfg.INTERPFACTOR)
        # Change this to also account for minimum distance
        leftpeaks = detect_peaks(left_interp, threshold=cfg.MINHEIGHT)
        rightpeaks = detect_peaks(right_interp, threshold=cfg.MINHEIGHT)
        # Change this to allow for min likelihood (add to config)
        pairs, (left_unmatched, right_unmatched) = match_peaks(
            leftpeaks, rightpeaks, left_interp, right_interp)
        hits = []
        hiterrs = []
        for pair in pairs:
            leftcfd, leftstatus = cfd_timestamp(left_interp, pair.left)
            rightcfd, rightstatus = cfd_timestamp(right_interp, pair.right)
            if leftstatus is False or rightstatus is False:
                print("Skipped due to CFD failure")
                continue
            xpos, xposerr = transverse_position(
                x_to_t(leftcfd, interpfactor=cfg.INTERPFACTOR),
                x_to_t(rightcfd, interpfactor=cfg.INTERPFACTOR))
            ypos = y_to_loc(get_centroid(left_interp, pair.left) +
                            get_centroid(right_interp, pair.right) / 2.0, interpfactor=cfg.INTERPFACTOR)
            recohit = RecoHit(xpos, ypos, (pair.left.x + pair.right.x) /
                              2.0, (pair.left.y + pair.right.y) / 2.0, (pair.left.z + pair.right.z) / 2.0)
            hiterr = Hit(xposerr, 5)
            hits.append(recohit)
            hiterrs.append(hiterr)
        self.deconvolved = (left_deconvolved, right_deconvolved)
        self.interped = (left_interp, right_interp)
        self.peaks = (leftpeaks, rightpeaks)
        self.pairs = pairs
        self.hits = hits
        self.hiterrors = hiterrs
        if plot:
            x = np.arange(0, 1014, 1)
            y = np.arange(0, 8, 1)
            xx, yy = np.meshgrid(x, y)
            thismin = min((np.min(self.leftmatrix), np.min(self.rightmatrix)))
            thismax = max((np.max(self.leftmatrix), np.max(self.rightmatrix)))
            fig = plt.figure()
            fig.suptitle(f"Event {levent.event_no}")
            ax0 = fig.add_subplot(231, projection="3d")
            ax1 = fig.add_subplot(232)
            ax2 = fig.add_subplot(233)
            ax3 = fig.add_subplot(234, projection="3d")
            ax4 = fig.add_subplot(235)
            ax5 = fig.add_subplot(236)
            base_img = ax1.imshow(
                self.leftmatrix[0:8], aspect="auto", interpolation="none", vmin=thismin, vmax=thismax)
            ax4.imshow(self.rightmatrix[0:8], aspect="auto",
                       interpolation="none", vmin=thismin, vmax=thismax)
            ax2.imshow(left_interp[0:80],
                       aspect="auto", interpolation="none", vmin=thismin, vmax=thismax)
            ax5.imshow(right_interp[0:80],
                       aspect="auto", interpolation="none", vmin=thismin, vmax=thismax)
            surf1 = ax0.plot_surface(yy, xx, self.leftmatrix[0:8],
                                     rstride=1, cstride=1, vmin=thismin, vmax=thismax, cmap=cm.coolwarm)
            surf2 = ax3.plot_surface(yy, xx, self.rightmatrix[0:8],
                                     rstride=1, cstride=1, vmin=thismin, vmax=thismax, cmap=cm.coolwarm)
            ax0.set_title("Left, observed signal")
            ax3.set_title("Right, observed signal")
            ax1.set_title("Left, observed signal")
            ax4.set_title("Right, observed signal")
            ax2.set_title("Left, deconvolved, interpolated")
            ax5.set_title("Right, deconvolved, interpolated")
            # To set the axis ticks:
            ax1.set_xticks(np.arange(0, 1014, 100))
            ax1.set_xticklabels(np.arange(0, 210, 20))
            ax1.set_yticks(np.arange(0, 8, 1))
            ax1.set_yticklabels(np.arange(14, 6, -1))
            ax4.set_xticks(np.arange(0, 1014, 100))
            ax4.set_xticklabels(np.arange(0, 210, 20))
            ax4.set_yticks(np.arange(0, 8, 1))
            ax4.set_yticklabels(np.arange(14, 6, -1))
            ax2.set_xticks(np.arange(0, 10140, 1000))
            ax2.set_xticklabels(np.arange(0, 210, 20))
            ax2.set_yticks(np.arange(0, 80, 10))
            ax2.set_yticklabels(np.arange(14, 6, -1))
            ax5.set_xticks(np.arange(0, 10140, 1000))
            ax5.set_xticklabels(np.arange(0, 210, 20))
            ax5.set_yticks(np.arange(0, 80, 10))
            ax5.set_yticklabels(np.arange(14, 6, -1))
            ax1.set_xlabel("Time (ns)")
            ax1.set_ylabel("Stripnumber")
            ax4.set_xlabel("Time (ns)")
            ax4.set_ylabel("Stripnumber")
            ax2.set_xlabel("Time (ns)")
            ax2.set_ylabel("Stripnumber")
            ax5.set_xlabel("Time (ns)")
            ax5.set_ylabel("Stripnumber")
            fig.colorbar(base_img, ax=[ax0, ax1, ax2, ax3, ax4, ax5],
                         orientation="horizontal", fraction=0.016)
            for i, pair in enumerate(pairs):
                ax2.scatter(pair.left.x, pair.left.y,
                            c="pink", s=25, marker="o")
                ax5.scatter(pair.right.x, pair.right.y,
                            c="pink", s=25, marker="o")
                # , c="green", s=4, marker="x")
                ax2.annotate(str(i), xy=(pair.left.x, pair.left.y), c="white")
                # , c="green", s=4, marker="x")
                ax5.annotate(str(i), xy=(
                    pair.right.x, pair.right.y), c="white")
            for lhit in left_unmatched:
                ax2.scatter(lhit.x, lhit.y, c="red", s=25, marker="x")
            for rhit in right_unmatched:
                ax5.scatter(rhit.x, rhit.y, c="red", s=25, marker="x")
            plt.show()

    def plot_hits(self):
        if not self.hits:
            print("Need to run reconstruct() first")
            return
        strip_positions = y_to_loc(np.arange(0, 28, 1), interpfactor=1)
        plt.errorbar([hit.recox for hit in self.hits], [hit.recoy for hit in self.hits], xerr=[hiterr.x for hiterr in self.hiterrors],
                     yerr=[hiterr.y for hiterr in self.hiterrors], marker="o", markersize=2.5, linestyle="None", capsize=2.5)
        for i in strip_positions:
            plt.axhline(i, c="purple", alpha=0.2)
            plt.ylim(0, 200)
            plt.xlim(-150, 150)
            plt.xlabel("Horizontal position (mm)")
            plt.ylabel("Vertical position (mm)")
        for i, (xpos, ypos, _, _, _) in enumerate(self.hits):
            plt.annotate(str(i), xy=(xpos+2.5, ypos+2.5), c="black")
        plt.show()

    def calculate_trigger_times(self):
        times = []
        if not self.hits:
            # print("No reconstructed hits")
            return times
        if not self.trigevent:
            print("No trigger event")
            return times
        if not self.trigevent.pulses:
            print("Did not find any pulses in this trigger event")
            return times
        try:
            trigger_time = self.trigevent.pulses[0].left.cfpeak
        except IndexError:
            breakpoint()
        for hit in self.hits:
            event_time = x_to_t(hit.x) - trigger_time
            times.append(event_time)
        return times

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
    times = []
    for levent in LAPPDEvent.itr_all_raw(base_dir, trigger="_0"):
        if np.max(levent.leftmatrix) > 70 or np.max(levent.rightmatrix) > 70 or np.max(levent.leftmatrix) < 4 or np.max(levent.rightmatrix) < 4:
            continue
        # if levent.event_no > 10000:
        #     break
        # for strip, sevent in levent.stripevents.items():
        #     for pulse in sevent.pulses:
        #         print(f"Strip {strip}: Position: {pulse.position}")
        levent.reconstruct()
        times += levent.calculate_trigger_times()
    # for levent, lpulses in LAPPDEvent.search_strip(stripnumber, base_dir):
    #     # psf = lpulses[0].leftmatrix[0:3, :]
    #     # levent.set_rootfile("testfile.root")
    #     # levent.write_root()
    #     lpulses[0].centroid_height(plot=True)
    #     for npulse in lpulses[0].strippulses:
    #         pulse = lpulses[0].strippulses[npulse].left
    #         if pulse.peak_present:
    #             cfpeaks.append(abs(pulse.cfpeak - pulse.peaktime))
    # myhist = root.TH1D("hist", "hist", 35, 0.5, 1.6)
    # for value in cfpeaks:
    #     myhist.Fill(value)
    # myhist.Draw()
    roothist(times, 1)
    breakpoint()
