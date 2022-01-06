import pdb

import matplotlib.pyplot as plt
import numpy as np


def mv_to_adc(millivolts: float, bit_r: int = 12) -> float:
    resolution = (2**bit_r) - 1
    return resolution * millivolts/1000.0


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
    if polarity == "positive":
        wave = -wave
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
        condition (str, optional): Greater than or less than width.
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
    breakpoint()
    if condition.lower() == "greater":
        crossing_width *= -1
        width *= -1
    if crossing_width < width:
        return True
    else:
        return False
