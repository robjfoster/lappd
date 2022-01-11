import operator
import os
import pdb
from array import array

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from scipy import signal
from scipy.optimize import curve_fit
from scipy.stats import moyal


def mv_to_adc(millivolts: float, bit_r: int = 12) -> float:
    resolution = (2**bit_r) - 1
    return resolution * millivolts/1000.0


def get_truth(inp, relate, cut):
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '==': operator.eq}
    return ops[relate](inp, cut)


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
    crossing_width = (back_crossing + forward_crossing) * ns_per_sample
    if plot:
        x = range(len(wave))
        plt.plot(x, wave)
        plt.axvline(peak_sample - back_crossing, c='purple')
        plt.axvline(peak_sample + forward_crossing, c='purple')
        plt.axhline(thresh, c='r')
        plt.show()
    if condition.lower() == "greater":
        crossing_width *= -1
        width *= -1
    if crossing_width < width:
        return True
    else:
        return False


def rise_time(wave: np.ndarray, time: float, condition: str = "less",
              ns_per_sample: float = 0.2, fraction1: float = 0.2,
              fraction2: float = 0.8, plot: bool = False) -> bool:
    peak_sample = np.argmin(wave)
    peak = wave[peak_sample]
    rise_samples = wave[:peak_sample]
    frac1_sample = peak_sample - \
        np.argmax(rise_samples[::-1] > peak * fraction1)
    frac2_sample = peak_sample - \
        np.argmax(rise_samples[::-1] < peak * fraction2)
    value = (frac2_sample - frac1_sample) * ns_per_sample
    if plot:
        x = range(len(wave))
        plt.plot(x, wave)
        plt.axvline(frac1_sample, c='purple')
        plt.axvline(frac2_sample, c='purple')
        plt.axhline(peak*fraction1, c='r')
        plt.axhline(peak*fraction2, c='r')
        plt.show()
    if condition.lower() == "greater":
        value *= -1
        time *= -1
    if value < time:
        return True
    else:
        return False


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
    #canvas = root.TCanvas("c1", "A Simple Graph Example", 200, 10, 500, 300)
    # breakpoint()
    #graph = root.TGraph(len(wave), xvec, wavevec)
    #graph.Fit("landau", "B")
    # graph.Draw("AC*")
    plt.plot(range(len(wave)), wave)
    plt.show()
    root.langau.langaus(wave)
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
