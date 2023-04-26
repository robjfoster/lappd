import sys

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root

from lappd.strip import StripEvent
from lappd.utils import lappdcfg as cfg
from lappd.utils import sigutils as su
from lappd.utils.roothist import roothist

base_dir = sys.argv[1]
strips = [int(i) for i in sys.argv[2].split(",")]

print(strips)

charges = tuple([] for i in strips)
heights = tuple([] for i in strips)
total_pulses = [0 for i in strips]
print(charges)

sample_start = int(80 / cfg.NS_PER_SAMPLE)
sample_end = int(100 / cfg.NS_PER_SAMPLE)

for sevents in zip(*[StripEvent.itr_file_raw(i, base_dir) for i in strips]):
    for i, sevent in enumerate(sevents):
        lwave = sevent.leftwaveform
        rwave = sevent.rightwaveform
        lwindow = lwave[sample_start:sample_end]
        rwindow = rwave[sample_start:sample_end]
        lpeak_time = (np.argmax(lwave < -8))
        rpeak_time = np.argmax(rwave < -8)
        peak_time = (lpeak_time+rpeak_time) / 2 * 0.2
        if peak_time > 1:
            if lpeak_time < 1:
                peak_time = rpeak_time
            if rpeak_time < 1:
                peak_time = lpeak_time
        inwindowpeak = int((np.argmin(lwindow) + np.argmin(rwindow))/2)
        start_int = inwindowpeak - 10 if inwindowpeak > 10 else 0
        end_int = inwindowpeak + \
            20 if (inwindowpeak + 20) < (sample_end -
                                         sample_start) else (sample_end-sample_start) if start_int != 0 else 30
        if end_int == (sample_end-sample_start):
            start_int = end_int - 30
        totalcharge = 0
        leftcharge = np.trapz(
            lwindow[start_int:end_int], dx=cfg.NS_PER_SAMPLE) * -1
        rightcharge = np.trapz(
            rwindow[start_int: end_int], dx=cfg.NS_PER_SAMPLE) * -1
        totalcharge = leftcharge + rightcharge
        lheight = -1*np.min(lwindow)
        rheight = -1*np.min(rwindow)
        height = (lheight+rheight) / 2
        heights[i].append(height)
        picocoulombs = totalcharge/1e3*1e-9/50/1.6e-19
        charges[i].append(picocoulombs)
        if np.min(lwindow) < -8 or np.min(rwindow) < -8:
            total_pulses[i] += 1

balances = []
for i, strip in enumerate(heights):
    if i == 0 or i == len(heights)-1:
        continue
    thisbalance = []
    for j, event in enumerate(strip):
        try:
            balance = (event / heights[i-1][j]) - (event / heights[i+1][j])
            thisbalance.append(balance)
        except IndexError:
            breakpoint()
    balances.append(thisbalance)

pos = base_dir.split("/")[-2]

plt.bar(strips, total_pulses)
# This can be used for adding a point on the plot denoting laser position
# otherax = plt.twiny()
# otherax.set_xticks(range(0, 41))
# otherax.set_xlim(0, 40)
# otherax.scatter(int(pos[1:]), max(total_pulses), c='red')
plt.xlabel("Strip")
plt.ylabel("Total pulses")
plt.savefig(f"data/vertical_scan/total_pulses_{pos}.png")
plt.show()
breakpoint()
plt.clf()
plt.bar(strips, [sum(i) for i in charges])
plt.xlabel("Strip")
plt.ylabel("Total charge")
plt.savefig(f"data/vertical_scan/total_charge_{pos}.png")
plt.clf()
plt.bar(strips, [np.mean(i) for i in charges])
plt.xlabel("Strip")
plt.ylabel("Mean charge")
plt.savefig(f"data/vertical_scan/mean_charge_{pos}.png")
plt.clf()
plt.bar(strips, [np.mean(i) for i in heights])
plt.xlabel("Strip")
plt.ylabel("Mean height")
plt.savefig(f"data/vertical_scan/mean_height_{pos}.png")
plt.clf()
plt.bar(strips[1:-1], [np.mean(i) for i in balances])
plt.xlabel("Strip")
plt.ylabel("'Balance'")
plt.savefig(f"data/vertical_scan/balance_{pos}.png")
plt.clf()
