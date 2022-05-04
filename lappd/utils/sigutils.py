import operator
import os
import pdb
from array import array
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from scipy import interpolate, signal
from scipy.optimize import curve_fit
from scipy.stats import moyal

from . import gimmedatwave as gdw


def mv_to_adc(millivolts: float, bit_r: int = 12) -> float:
    resolution = (2**bit_r) - 1
    return resolution * millivolts/1000.0


def adc_to_mv(adc: float, bit_r: int = 12) -> float:
    resolution = (2**bit_r) - 1
    return adc / resolution * 1000


def get_truth(inp, relate, cut):
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '==': operator.eq}
    return ops[relate](inp, cut)


def ns_to_sample(ns, ns_per_sample):
    return int(ns/ns_per_sample)


def reduced_mean(wave: np.ndarray, plot: bool = False) -> float:
    dev = np.std(wave)
    avg = np.mean(wave)
    reduced = wave[np.abs(wave) < avg + dev]
    value = np.mean(reduced)
    if plot:
        plt.plot(range(len(wave)), wave)
        plt.axhline(value + dev, c='orange', alpha=0.4)
        plt.axhline(value - dev, c='orange', alpha=0.4)
        plt.axhline(value, c='r', alpha=1)
        plt.axhline(avg, c='purple', alpha=1)
        plt.show()
    return value


def slice_around_peak(wave, lookback, lookforward, ns_per_sample=0.2):
    peak = np.argmin(wave)
    lookback_samples = peak-int(lookback / ns_per_sample)
    lookforward_samples = peak+int(lookforward / ns_per_sample)
    if lookback_samples < 0:
        lookback_samples = 0
    if lookforward_samples > len(wave):
        lookforward_samples = -1
    return wave[lookback_samples:lookforward_samples]


def slice_around_sample(wave, peak_sample, lookback, lookforward, ns_per_sample=0.2):
    lookback_samples = peak_sample-int(lookback / ns_per_sample)
    lookforward_samples = peak_sample+int(lookforward / ns_per_sample)
    if lookback_samples < 0:
        lookback_samples = 0
    if lookforward_samples > len(wave):
        lookforward_samples = -1
    return wave[lookback_samples:lookforward_samples]


def threshold(wave: np.ndarray, thresh: float, polarity: str = "negative") -> bool:
    """Any sample in wave must pass below threshold to return True
       Use combination of threshold and !threshold for min/max voltage cuts"""
    if polarity == "positive":
        wave = -wave
        thresh = -thresh
    return True if np.any(wave < thresh) else False


def pulse_width_cf(wave: np.ndarray, width: float, ns_per_sample: float = 0.2,
                   fraction: float = 0.2, condition: str = "less", plot: bool = False
                   ) -> bool:
    """Constant fraction pulse width cut

    Args:
        wave (np.ndarray): Waveform or reduced waveform
        width (float): Pulse width threshold value in ns
        ns_per_sample (float, optional): Digitiser sampling rate 
                                         (200ps for V1742). Defaults to 0.2.
        fraction (float, optional): Fraction of total pulse height at which to 
                                    calculate pulse width. Defaults to 0.2.
        condition (str, optional): "greater" or "less" than width.
                                   Defaults to "less".
        plot (bool, optional): Plot waveform with pulse width and CF indicated.
                               Defaults to False.

    Returns:
        bool: True if passed cut, False otherwise
    """
    # pass in a cut wave and pulse width threshold
    peak_sample = np.argmin(wave)
    thresh = wave[peak_sample] * fraction
    # samples either side of peak
    width_samples = int(round(width / ns_per_sample))
    if (peak_sample + width_samples > len(wave)) or (peak_sample - width_samples < 0):
        return False
    back_samples = wave[peak_sample-width_samples: peak_sample]
    forward_samples = wave[peak_sample: peak_sample+width_samples]
    back_crossing = np.argmax(back_samples[::-1] > thresh)
    forward_crossing = np.argmax(forward_samples > thresh)
    breakpoint()
    crossing_width = (back_crossing + forward_crossing) * ns_per_sample
    if plot:
        x = range(len(wave))
        plt.plot(x, wave)
        plt.axvline(peak_sample - back_crossing, c='purple', alpha=0.2)
        plt.axvline(peak_sample + forward_crossing, c='purple', alpha=0.2)
        plt.axhline(thresh, c='r', alpha=0.2)
        plt.show()
    if condition.lower() == "greater":
        crossing_width *= -1
        width *= -1
    print("measured width:", crossing_width)
    print("Comparing to: ", width)
    if crossing_width < width:
        return True
    else:
        return False


def pulse_width_new(wave: np.ndarray, width: float, ns_per_sample: float = 0.2,
                    fraction: float = 0.2, condition: str = "less", plot: bool = False
                    ) -> bool:
    """Constant fraction pulse width cut

    Args:
        wave (np.ndarray): Waveform or reduced waveform
        width (float): Pulse width threshold value in ns
        ns_per_sample (float, optional): Digitiser sampling rate 
                                         (200ps for V1742). Defaults to 0.2.
        fraction (float, optional): Fraction of total pulse height at which to 
                                    calculate pulse width. Defaults to 0.2.
        condition (str, optional): "greater" or "less" than width.
                                   Defaults to "less".
        plot (bool, optional): Plot waveform with pulse width and CF indicated.
                               Defaults to False.

    Returns:
        bool: True if passed cut, False otherwise
    """
    # pass in a cut wave and pulse width threshold
    if condition not in ["greater", "less"]:
        raise ValueError
    peak_sample = np.argmin(wave)
    thresh = wave[peak_sample] * fraction
    # samples either side of peak
    width_samples = int(round(width / ns_per_sample))
    if (peak_sample + width_samples > len(wave)) or (peak_sample - width_samples < 0):
        return False
    back_samples = wave[peak_sample-width_samples: peak_sample]
    forward_samples = wave[peak_sample: peak_sample+width_samples]
    back_crossing = np.argmax(back_samples[::-1] > thresh)
    forward_crossing = np.argmax(forward_samples > thresh)
    crossing_width = (back_crossing + forward_crossing) * ns_per_sample
    print("measured width:", crossing_width)
    print("Comparing to: ", width)
    if plot:
        x = range(len(wave))
        plt.plot(x, wave)
        plt.axvline(peak_sample - back_crossing, c='purple', alpha=0.2)
        plt.axvline(peak_sample + forward_crossing, c='orange', alpha=0.2)
        plt.axhline(thresh, c='r', alpha=0.2)
        plt.show()
    # if wave has not reached threshold within pulse width, value will be zero
    if (back_crossing == 0) or (forward_crossing == 0):
        if condition == "greater":
            return True
        else:
            return False
    else:
        if condition == "less":
            return True
        else:
            return False


def rise_time(wave: np.ndarray, time: float, condition: str = "less",
              ns_per_sample: float = 0.2, fraction1: float = 0.2,
              fraction2: float = 0.8, plot: bool = False) -> bool:
    peak_sample = np.argmin(wave)
    peak = wave[peak_sample]
    rise_samples = wave[:peak_sample]
    if len(rise_samples) == 0:
        return False
    frac1_sample = peak_sample - \
        np.argmax(rise_samples[::-1] > peak * fraction1)
    frac2_sample = peak_sample - \
        np.argmax(rise_samples[::-1] < peak * fraction2)
    value = (frac2_sample - frac1_sample) * ns_per_sample
    if plot:
        x = range(len(wave))
        plt.plot(x, wave)
        plt.axvline(frac1_sample, c='purple', alpha=0.2)
        plt.axvline(frac2_sample, c='purple', alpha=0.2)
        plt.axhline(peak*fraction1, c='r', alpha=0.2)
        plt.axhline(peak*fraction2, c='r', alpha=0.2)
        plt.show()
    if condition.lower() == "greater":
        value *= -1
        time *= -1
    if value < time:
        return True
    else:
        return False


def cfd(wave, fraction, times=None, userpeak=None, plot=False, samplesabovethresh=1):
    peak_sample = np.argmin(wave) if not userpeak else userpeak
    # print(peak_sample)
    peak_value = wave[peak_sample]
    threshold = peak_value * fraction
    beforepeak = wave[:peak_sample]
    thresh_sample = None
    for i, sample in enumerate(beforepeak[::-1]):
        if sample > threshold:
            if np.all(beforepeak[i:i+samplesabovethresh] > threshold):
                thresh_sample = peak_sample - i
                # thresh_sample = ((peak_sample - i) +
                #                  (peak_sample - (i+1))) / 2.0
                break
    if thresh_sample is None:
        return None
    if plot and thresh_sample:
        plt.plot(wave, "o", linestyle="-", markersize=1)
        plt.axvline(thresh_sample, c="g")
        plt.axvline(peak_sample, c="purple")
        plt.axhline(threshold, c="r")
        plt.show()
    if times is None:
        return thresh_sample
    else:
        return (times[thresh_sample] + times[thresh_sample + 1]) / 2.0


def minmax(wave: np.ndarray, max_minus_min: float, condition: str = "greater") -> bool:
    """Calculate and cut on max and min of baseline subtracted wave"""
    wave_max = np.max(wave)
    wave_min = np.min(wave)
    difference = np.abs(wave_max - wave_min)
    if condition.lower() == "less":
        difference *= -1
        max_minus_min *= -1
    if (difference) > max_minus_min:
        return True
    else:
        return False


def draw_fft(wave: np.ndarray, ns_per_sample: float = 0.2) -> None:
    x = np.fft.fftfreq(len(wave), ns_per_sample*1e-9)
    y = np.fft.fft(wave)
    y = np.fft.fftshift(y / len(y))
    plt.plot(x, np.abs(y))
    plt.show()


def draw_rfft(wave: np.ndarray, ns_per_sample: float = 0.2) -> None:
    x = np.fft.rfftfreq(len(wave), ns_per_sample*1e-9)
    y = np.fft.rfft(wave)
    #y = np.fft.fftshift(y / len(y))
    plt.plot(x, np.abs(y))
    plt.show()


def butterworth_lowpass(wave: np.ndarray, freq: float, ns_per_sample: float = 0.2) -> np.ndarray:
    nyquist = 0.5 * 1 / (ns_per_sample * 1e-9)
    b, a = signal.butter(2, freq/nyquist, btype='low', analog=False)
    filt_wave = signal.filtfilt(b, a, wave)
    return filt_wave


def fit_landau(wave):
    # for some reason, only python Array objects can be cast to std::vectors
    # with root python bindings?
    #wavevec = root.std.vector("float")()
    #xvec = root.std.vector("float")()
    #wavevec = root.TVector("float")
    #xvec = root.TVector("float")
    root.gROOT.SetMacroPath(os.getcwd() + "/lappd/utils/")
    root.gROOT.LoadMacro("langau.C")
    wave = wave * -1
    wavevec = array('d')
    xvec = array('d')
    for i, value in enumerate(wave):
        # wavevec.push_back(value)
        # xvec.push_back(i)
        wavevec.append(value)
        xvec.append(i)
    canvas = root.TCanvas("c1", "A Simple Graph Example", 200, 10, 500, 300)
    # breakpoint()
    graph = root.TGraph(len(wave), xvec, wavevec)
    graph.Fit("landau", "B")
    graph.Draw("AC*")
    #plt.plot(range(len(wave)), wave)
    # plt.show()
    # root.langau.langaus(wave)
    breakpoint()


def fmoyal(x, a, b):
    return -moyal.pdf(x, loc=a, scale=b)


def scan_other(otherwave: np.ndarray,
               peakx: int,
               peaky: float,
               scanrange: int = 10,
               ns_per_sample: float = 0.2,
               peak_height_req: float = 0.8,
               plot: bool = False
               ) -> bool:
    """Scan another wave for a corresponding peak that reaches a fraction of
    at least peak_height_req of height of original peak. scanrange is in ns"""
    scansamples = int(scanrange / ns_per_sample)
    leftscan = peakx - scansamples
    rightscan = peakx + scansamples + 1
    if leftscan < 0:
        leftscan = 0
    if rightscan > len(otherwave):
        rightscan = len(otherwave)
    scanarea = otherwave[leftscan:rightscan]
    if np.min(scanarea) < peak_height_req * peaky:
        return True
    else:
        return False


def peak_find(wave, peakparams, ns_per_sample=0.2):
    # Call scipy.find_peaks
    peaks, peakinfo = signal.find_peaks(
        wave*-1.0,
        height=((peakparams.getfloat('minheight')),
                peakparams.getfloat('maxheight')),
        distance=peakparams.getfloat('mindistance') / ns_per_sample,
        width=(peakparams.getfloat('minwidth') / ns_per_sample,
               peakparams.getfloat('maxwidth') / ns_per_sample),
        rel_height=peakparams.getfloat('relheight'))
    return peaks, peakinfo


def coarse_correlate(leftwave,
                     leftpeaks,
                     rightwave,
                     rightpeaks,
                     max_offset,
                     minrelpulseheight=0.5,
                     maxrelpulseheight=1.5,
                     ns_per_sample=0.2):
    # An attempt at correlating peaks on two waveform produced from either side of a stripline
    # Looks for pulses within max_offset (ns) and chooses the one with the largest pulse height
    # Probably some edge cases that have been missed, but scipy.find_peaks ensures that there
    # should not be multiple identified peaks within a given distance anyway
    peak_pairs = []
    offsets = []
    leftpeakheights = leftwave[leftpeaks]
    rightpeakheights = rightwave[rightpeaks]
    # Loop through the shorter of the two lists
    first_itr = leftpeaks if len(leftpeaks) <= len(rightpeaks) else rightpeaks
    first_info = leftpeakheights if first_itr is leftpeaks else rightpeakheights
    sec_itr = leftpeaks if first_itr is not leftpeaks else rightpeaks
    sec_info = leftpeakheights if first_info is not leftpeaks else rightpeakheights
    for i, ipeak in enumerate(first_itr):
        potential_partners = []
        partner_offsets = []
        partner_heights = []
        for j, jpeak in enumerate(sec_itr):
            offset = (ipeak - jpeak) * ns_per_sample
            # if first_itr is rightpeaks then ipeak is right end of strip
            # we want to do left - right
            offset = offset * -1 if first_itr is rightpeaks else offset
            # check that pulses are within max_offset of each other in time
            if abs(offset) < max_offset:
                # check that the pulse height is within rel_threshold
                relative_height = first_info[i] / sec_info[j]
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


def correlate_peaks(leftpeaks,
                    leftpeakinfo,
                    rightpeaks,
                    rightpeakinfo,
                    max_offset,
                    minrelpulseheight=0.5,
                    maxrelpulseheight=1.5,
                    ns_per_sample=0.2):
    # An attempt at correlating peaks on two waveform produced from either side of a stripline
    # Looks for pulses within max_offset (ns) and chooses the one with the largest pulse height
    # Probably some edge cases that have been missed, but scipy.find_peaks ensures that there
    # should not be multiple identified peaks within a given distance anyway
    peak_pairs = []
    offsets = []
    # Loop through the shorter of the two lists
    first_itr = leftpeaks if len(leftpeaks) <= len(rightpeaks) else rightpeaks
    first_info = leftpeakinfo if first_itr is leftpeaks else rightpeakinfo
    sec_itr = leftpeaks if first_itr is not leftpeaks else rightpeaks
    sec_info = leftpeakinfo if first_info is not leftpeakinfo else rightpeakinfo
    for i, ipeak in enumerate(first_itr):
        potential_partners = []
        partner_offsets = []
        partner_heights = []
        for j, jpeak in enumerate(sec_itr):
            offset = (ipeak - jpeak) * ns_per_sample
            # if first_itr is rightpeaks then ipeak is right end of strip
            # we want to do left - right
            offset = offset * -1 if first_itr is rightpeaks else offset
            # check that pulses are within max_offset of each other in time
            if abs(offset) < max_offset:
                # check that the pulse height is within rel_threshold
                relative_height = first_info['peak_heights'][i] / \
                    sec_info['peak_heights'][j]
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


def correlate_peaks_new(leftwave, rightwave,
                        leftpeaks,
                        leftpeakinfo,
                        rightpeaks,
                        rightpeakinfo,
                        max_offset,
                        minrelpulseheight=0.5,
                        maxrelpulseheight=1.5,
                        ns_per_sample=0.2):
    peak_pairs = []
    offsets = []
    potential_pairs = []
    potential_offsets = []
    potential_heights = []
    # Loop through left peaks first
    for i, lpeak in enumerate(leftpeaks):
        for j, rpeak in enumerate(rightpeaks):
            offset = (lpeak - rpeak) * ns_per_sample
            if abs(offset) < max_offset:
                # check that the pulse height is within rel_threshold
                relative_height = leftpeakinfo['peak_heights'][i] / \
                    rightpeakinfo['peak_heights'][j]
                if minrelpulseheight < relative_height < maxrelpulseheight:
                    potential_pairs.append((lpeak, rpeak))
                    potential_offsets.append(offset)
                    potential_heights.append(relative_height)
    conflicts = []
    # if len(potential_pairs) == 1:
    #    return peak_pairs, offsets
    for pairpair in combinations(potential_pairs, 2):
        print(pairpair)
        if abs(pairpair[0][0] - pairpair[1][0]) < max_offset or abs(pairpair[0][1] - pairpair[1][1]) < max_offset:
            conflicts.append(pairpair)
    if len(conflicts) > 1:
        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, sharey=True)
        fig.set_size_inches(14.5, 6.5)
        ax1.plot(gdw.ns(n_samples=len(leftwave)), leftwave)
        ax1.plot(leftpeaks*0.2, leftwave[leftpeaks], "x")
        ax2.plot(gdw.ns(n_samples=len(rightwave)), rightwave)
        ax2.plot(rightpeaks*0.2, rightwave[rightpeaks], "x")
        plt.show()
        pdb.set_trace()
    # for i, ipair in enumerate(potential_pairs):
    #     for j, jpair in enumerate(potential_pairs):
    #         if i == j:
    #             continue
    #         else:
    #             if abs(ipair[0] - jpair[0]) < max_offset or abs(ipair[1] - jpair[1]) < max_offset:
    #                 pdb.set_trace()
    #     print("cool")
    #     pdb.set_trace()
    # return
    return peak_pairs, offsets


def cubic_spline(x, y, ratio=3, plot=False):
    newx = np.linspace(x[0], x[-1], ratio * len(x))
    interpfunc = interpolate.interp1d(x, y, kind="cubic")
    if plot:
        plt.plot(x, y, 'o', newx, interpfunc(newx), "-")
        plt.show()
    return newx, interpfunc(newx)
