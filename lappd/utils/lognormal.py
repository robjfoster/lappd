import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import ROOT as root


def lognormal(t, B, Q, t0, m, sigma):
    return B + (Q / ((t - t0) * np.sqrt(2 * np.pi) * sigma)) * np.exp((-0.5 * (np.log((t - t0) / m) / sigma)**2))


def lognormal2(t, tau, sigma, scale):
    return (scale / (t * sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (np.log(t / tau) / sigma)**2)


def lognormal3(t, Q, t0, tau, sigma):
    return - (Q / ((t + t0) * np.sqrt(2 * np.pi * sigma**2))) * np.exp(-(np.log((t + t0) / tau))**2 / (2 * sigma**2))


def lognormal4(t, Q, tau, sigma):
    return - (Q / ((t) * np.sqrt(2 * np.pi * sigma**2))) * np.exp(-(np.log((t) / tau))**2 / (2 * sigma**2))


class LogNormal:

    def __call__(self, arr, par):
        return par[3] - (par[0] / ((arr[0]) * np.sqrt(2 * np.pi * par[2]**2))) * np.exp(-(np.log((arr[0]) / par[1]))**2 / (2 * par[2]**2))

    def npcall(self, arr, par):
        return par[3] - (par[0] / ((arr) * np.sqrt(2 * np.pi * par[2]**2))) * np.exp(-(np.log((arr) / par[1]))**2 / (2 * par[2]**2))


class ScintPDF:

    def __call__(self, arr, par):
        # arr[0] = t
        # par[0] = A1
        # par[1] = tau1
        # par[2] = tauR
        # par[3] = tau2
        return (par[0] * ((np.exp(-arr[0]/par[1]) - np.exp(-arr[0]/par[2])) / (par[1] - par[2])) +
                (1 - par[0]) * ((np.exp(-arr[0]/par[3]) - np.exp(-arr[0]/par[2])) / (par[3] - par[2])))


def root_ln(x, y):
    maxtime = np.max(x)
    mintime = np.min(x)
    # c1 = root.TCanvas("c1", "Simple graph example", 200, 10, 700, 500)
    gr = root.TGraph(len(x), x.astype(
        "float64"), y.astype("float64"))
    lognormal = LogNormal()
    tf1 = root.TF1("pyln", lognormal, mintime, maxtime, 4)
    tf1.SetParLimits(0, 0, 200)
    tf1.SetParLimits(1, mintime, maxtime)
    tf1.SetParLimits(2, 0.0, 3)
    tf1.SetParLimits(3, -5, 5)
    tf1.SetParNames("Q", "Centre", "Sigma", "Baseline")
    tf1.SetParameters(25, (maxtime + mintime) / 2, 0.05, 0)
    gr.Fit("pyln", "Q")
    # gr.Draw("AP")
    # gr.SetMarkerSize(0.5)
    # gr.SetMarkerStyle(20)
    # c1.Update()
    return tf1


def fitscint(hist):
    scintpdf = ScintPDF()
    tf1gaus = root.TF1("gauss", "gaus(0)", -5, 5)
    tf1scint = root.TF1("scintpdf", scintpdf, -10, 30, 4)
    scintxgausconv = root.TF1Convolution("scintpdf", "gauss", -10, 30, True)
    scintxgausconv.SetNofPointsFFT(10000)
    tf1scintxgaus = root.TF1(
        "scintxgaus", scintxgausconv, -10, 30, scintxgausconv.GetNpar())
    tf1scintxgaus.SetParNames("A1", "tau1", "tauR",
                              "tau2", "gausscale", "mean", "sigma")
    tf1scintxgaus.SetParameters(0.68, 3.0, 1.7, 11, 1, 0, 0.25)
    tf1scintxgaus.SetParLimits(0, 0.6, 0.7)
    tf1scintxgaus.SetParLimits(1, 2.0, 4.0)
    tf1scintxgaus.SetParLimits(2, 1.0, 2.0)
    tf1scintxgaus.SetParLimits(3, 10, 12)
    tf1scintxgaus.SetParLimits(4, 1, 1)
    tf1scintxgaus.SetParLimits(5, -1, 2)
    tf1scintxgaus.SetParLimits(6, 0.1, 0.4)
    # hist.Fit("gaus")
    hist.Fit("scintxgaus")
    hist.Draw()
    breakpoint()


def lnfit(lpulse, upsample_rate=10, debug=True):
    times = lpulse.times
    new_times = np.arange(np.min(times), np.max(times),
                          (times[1] - times[0]) / 10)
    data = lpulse.wave
    bounds = ((0, np.min(times), 0, 0), (200, np.max(times), 200, 3))
    p0 = (10, times[int(len(times) / 2)], 10, 0.8)
    popt, pcov = curve_fit(lognormal3, times, data, p0=p0, bounds=bounds)
    plt.plot(times, data, "x")
    plt.plot(new_times, lognormal3(new_times, *popt))
    if debug:
        print(
            f"Q: {popt[0]}, t0: {popt[1]}, tau: {popt[2]}, sigma: {popt[3]}, ")
        print(f"Integrated: {np.trapz(data, x=times)}")
        print(f"Centre: {popt[2] - popt[1]}")
    plt.show()
    breakpoint()
    return


def lnfit4(lpulse, upsample_rate=10, debug=True):
    times = lpulse.times
    new_times = np.arange(np.min(times), np.max(times),
                          (times[1] - times[0]) / 10)
    data = lpulse.wave
    bounds = ((0, np.min(times), 0), (200, np.max(times), 3))
    p0 = (10, times[int(len(times) / 2)], 0.8)
    popt, pcov = curve_fit(lognormal4, times, data, p0=p0, bounds=bounds)
    # plt.plot(times, data, "x")
    # plt.plot(new_times, lognormal4(new_times, *popt))
    # if debug:
    #     print(
    #         f"Q: {popt[0]}, tau: {popt[1]}, sigma: {popt[2]}, ")
    #     print(f"Integrated: {np.trapz(data, x=times)}")
    #     print(f"Centre: {popt[2] - popt[1]}")
    # plt.show()
    root_ln(times, data)
    return


if __name__ == "__main__":
    tvals = np.arange(0.2, 30, 0.2)
    # breakpoint()
    # yvals = lognormal(tvals, 0, 0.4, -2, 8.5, 69)
    # yvals = lognormal2(tvals, 5, 0.1, 20)
    yvals = lognormal3(tvals, 15, 1, 25, 0.05)
    print(f"STD: {np.log(np.std(yvals))}")
    print(f"Mean: {np.mean(yvals)}")
    print(f"Median: {np.median(yvals)}")
    print(f"Integrated: {np.trapz(yvals, x=tvals)}")
    plt.plot(tvals, yvals)
    plt.show()
