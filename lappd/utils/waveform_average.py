import sys

import matplotlib.pyplot as plt
import numpy as np

from ..lappdevent import LAPPDEvent

base_dir = sys.argv[1]
stripnumber = int(sys.argv[2])

avg_wave = np.zeros(150)
num_waves = 0

for levent in LAPPDEvent.search_all(base_dir):
    stripevent = levent.stripevents[stripnumber]
    for strippulse in stripevent.pulses:
        leftpulse = strippulse.left
        rightpulse = strippulse.right
        avg_wave += leftpulse.smoothedwave
        avg_wave += rightpulse.smoothedwave
        num_waves += 1
avg_wave /= (num_waves * 2)
plt.plot(avg_wave)
plt.show()
