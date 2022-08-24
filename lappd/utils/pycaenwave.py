# struct CAENHeader
# {
#     unsigned int eventSize;  // Should be 4120 for x1742 models (6*4byte header + 1024*4byte record)
#     unsigned int boardID;    // Not used
#     unsigned int pattern;    // Not used
#     unsigned int channel;    // This is the channel within the group i.e. ch8 is 0 but in group2
#     unsigned int eventCount; // Not exact event count as it rolls over whenever the limit is reached
#     unsigned int eventTime;  // Time when event is created in digitiser memory, not trigger time
# };

# struct CAENWave
# {
#     CAENHeader header;
#     std::vector<float> wave;
#     std::vector<float> times;
#     int eventNo;
# };
from dataclasses import dataclass

import numpy as np


@dataclass
class pyCAENHeader:
    eventSize: int
    boardID: int
    pattern: int
    channel: int
    eventCount: int
    eventTime: int

# This is just a class to emulate the CAENWave class from the cxx code to make
# simulated waveforms work with the analysis. Eventually the analysis should be
# refactored to make this redundant


@dataclass
class pyCAENWave:
    header: pyCAENHeader
    wave: np.ndarray
    times: np.ndarray
    eventNo: int
