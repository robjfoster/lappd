import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def lognormal(t, B, Q, t0, m, sigma):
    return B + (Q / ((t - t0) * np.sqrt(2 * np.pi) * sigma)) * np.exp((-0.5 * (np.log((t - t0) / m) / sigma)**2))


def lognormal2(t, tau, sigma, scale):
    return (scale / (t * sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (np.log(t / tau) / sigma)**2)


def lognormal3(t, Q, t0, tau, sigma):
    return - (Q / ((t + t0) * np.sqrt(2 * np.pi * sigma**2))) * np.exp(-(np.log((t + t0) / tau))**2 / (2 * sigma**2))


def lnfit(lpulse, debug=False):
    times = lpulse.times
    data = lpulse.wave
    bounds = ((0, min(times), 0, 0), (200, max(times), 200, 3))
    p0 = (10, times[int(len(times) / 2)], 10, 0.8)
    popt, pcov = curve_fit(lognormal3, times, data, p0=p0, bounds=bounds)
    plt.plot(times, data, "x")
    plt.plot(times, lognormal3(times, *popt))
    if debug:
        print(
            f"Q: {popt[0]}, t0: {popt[1]}, tau: {popt[2]}, sigma: {popt[3]}, ")
        print(f"Integrated: {np.trapz(data, x=times)}")
        print(f"Centre: {popt[2] - popt[1]}")
    plt.show()
    return


if __name__ == "__main__":
    tvals = np.arange(0.2, 30, 0.2)
    # breakpoint()
    #yvals = lognormal(tvals, 0, 0.4, -2, 8.5, 69)
    #yvals = lognormal2(tvals, 5, 0.1, 20)
    yvals = lognormal3(tvals, 15, 1, 25, 0.05)
    print(f"STD: {np.log(np.std(yvals))}")
    print(f"Mean: {np.mean(yvals)}")
    print(f"Median: {np.median(yvals)}")
    print(f"Integrated: {np.trapz(yvals, x=tvals)}")
    plt.plot(tvals, yvals)
    plt.show()
