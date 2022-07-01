import sys

import matplotlib.pyplot as plt
import numpy as np

from lappd.lappdevent import LAPPDEvent

base_dir = sys.argv[1]
stripnumber = int(sys.argv[2])

avg_wave = np.zeros(30)
num_waves = 0
peak_value = -20

for levent in LAPPDEvent.search_all(base_dir):
    stripevent = levent.stripevents[stripnumber]
    for strippulse in stripevent.pulses:
        leftpulse = strippulse.left
        rightpulse = strippulse.right
        leftpeakval = leftpulse.wave[leftpulse.rawpeak]
        rightpeakval = rightpulse.wave[rightpulse.rawpeak]
        leftscaledwave = leftpulse.wave * (peak_value / leftpeakval)
        rightscaledwave = rightpulse.wave * (peak_value / rightpeakval)
        avg_wave += leftscaledwave
        avg_wave += rightscaledwave
        num_waves += 1
avg_wave /= (num_waves * 2)
breakpoint()
plt.plot(avg_wave)
plt.title(f"Average of {num_waves} waveforms")
plt.xlabel("Time (ns)")
plt.ylabel("Amplitude (mV)")
plt.show()
