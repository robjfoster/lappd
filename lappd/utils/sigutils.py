import operator
import pdb

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


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


def minmax(wave: np.ndarray, max_minus_min: float, condition="greater") -> bool:
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
