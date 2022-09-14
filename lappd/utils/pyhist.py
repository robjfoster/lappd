import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from scipy.optimize import curve_fit


def rootgaus(x, constant, mean, sigma):
    return constant * np.exp(-0.5 * ((x - mean) / sigma)**2)


def pyhist(data, fit=False, binwidth=1, xlims=None):
    mindata = min(data)
    maxdata = max(data)
    bins = np.arange(mindata, maxdata + binwidth, binwidth)
    plt.hist(data, bins=bins, histtype="step")
    if fit:
        # Fit with root, draw with python
        tf1 = root.TF1("gaussian", "gaus", mindata, maxdata)
        hist = root.TH1D("pyhist", "pyhist", len(bins), mindata, maxdata)
        for value in data:
            hist.Fill(value)
        hist.Draw()
        hist.Fit(tf1)
        fitconstant = tf1.GetParameter(0)
        fitmean = tf1.GetParameter(1)
        fitsigma = tf1.GetParameter(2)
        fit = rootgaus(bins, fitconstant, fitmean, fitsigma)
        plt.plot(bins, fit)
        del hist
    if xlims:
        plt.xlim(*xlims)
    plt.show()
