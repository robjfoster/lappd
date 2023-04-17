import sys

import ROOT as root

filepath = sys.argv[1]
f = root.TFile(filepath, "READ")
c1 = f.Get("c1")
hist = c1.GetPrimitive("hist")
g1 = root.TF1("pedestal", "gaus", -2e6, 2e6)
g1.SetParameter(1, 0)
g1.SetLineColor(2)  # red
g1.SetLineStyle(2)
g2 = root.TF1("singlePE", "gaus", -2e6, 15e6)
g2.SetParameter(1, 10e6)
g2.SetParLimits(1, 2e6, 12e6)
g2.SetLineColor(3)  # green
g2.SetLineStyle(2)
hist.Fit("pedestal", "R")
hist.Fit("singlePE", "R+")
g3 = root.TF1("doublePE", "gaus", 4e6, 30e6)
g3.SetParameter(1, g2.GetParameter(1)*2)
g3.SetParameter(0, g2.GetParameter(0)/10)
g3.SetParLimits(1, g2.GetParameter(1)*2-1e6, g2.GetParameter(1)*2+1e6)
g3.SetLineColor(4)  # blue
g3.SetLineStyle(2)
hist.Fit("doublePE", "R+")
total = root.TF1("total", "gaus(0)+gaus(3)+gaus(6)", -1.5e6, 25e6)
pars = []
for i in range(3):
    if i == 1:
        pars.append(0)
    else:
        pars.append(g1.GetParameter(i))
for i in range(3):
    if i == 1:
        pars.append(7e6)
    else:
        pars.append(g2.GetParameter(i))
for i in range(3):
    if i == 1:
        pars.append(14e6)
    else:
        pars.append(g3.GetParameter(i))
total.SetParameters(*pars)
total.SetLineColor(6)
total.SetLineStyle(1)
hist.Fit("total", "R+")
hist.SetFillColorAlpha(9, 0.35)
c1.Draw()
breakpoint()
