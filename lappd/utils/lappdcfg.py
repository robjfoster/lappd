import configparser
import re
import sys

import numpy as np

filepath = "config/L104setup.cfg"
config = configparser.ConfigParser()
config.read(filepath)
if not config.sections():
    sys.exit("Config file is empty (probably wrong file path)")

peakparams = config['PEAKPARAMS']
daqconfig = config['DAQCONFIG']
analysis = config['ANALYSIS']
constants = config['CONSTANTS']

NS_PER_SAMPLE = daqconfig.getfloat("nspersample")
NSAMPLES = daqconfig.getint("nsamples")
TIMEJITTER = daqconfig.getfloat("timejitter")
VOLTAGEJITTER = daqconfig.getfloat("voltagejitter")
NREMOVEDSAMPLES = daqconfig.getint("nremovedsamples")

MINHEIGHT = peakparams.getfloat("minheight")  # mV
MINDISTANCE = peakparams.getfloat("mindistance")
RELHEIGHT = peakparams.getfloat("relheight")
MINRELPULSEHEIGHT = peakparams.getfloat("minrelpulseheight")
MAXRELPULSEHEIGHT = peakparams.getfloat("maxrelpulseheight")
MAXOFFSET = peakparams.getfloat("maxoffset")
INTERPFACTOR = peakparams.getint("interpfactor")
LOOKBACK = peakparams.getfloat("slicelookback")
LOOKFORWARD = peakparams.getfloat("slicelookforward")

TEMPLATE = np.load(analysis.get("template"))
TIMEPDFSIGMA = analysis.getfloat("timepdfsigma")
TIMEPDFSTART = analysis.getfloat("timepdfstart")
TIMEPDFEND = analysis.getfloat("timepdfend")
LOCPDFSIGMA = analysis.getfloat("locpdfsigma")
LOCPDFSTART = analysis.getfloat("locpdfstart")
LOCPDFEND = analysis.getfloat("locpdfend")
AMPLPDFSIGMA = analysis.getfloat("amplpdfsigma")
AMPLPFDFSTART = analysis.getfloat("amplpdfstart")
AMPLPDFEND = analysis.getfloat("amplpdfend")
PROPERROR = analysis.getfloat("properror")
MINLIKE = analysis.getfloat("minlike")

SOL = constants.getfloat("sol")
STRIPLENGTH = constants.getfloat("striplength")
STRIPVELOCITY = constants.getfloat("stripvelocity")
STRIPWIDTH = constants.getfloat("stripwidth")
STRIPSPACING = constants.getfloat("stripspacing")

MAXDELTA = 200 / (STRIPVELOCITY * SOL)  # Max LR offset in ns
MAXDELTAERROR = (PROPERROR / STRIPVELOCITY) * MAXDELTA
TIMES = np.arange(NS_PER_SAMPLE, NS_PER_SAMPLE *
                  (NSAMPLES - NREMOVEDSAMPLES), NS_PER_SAMPLE)
# The average radius of electron charge cloud exiting the bottom MCP (mm)
CLOUDRADIUS = 4

allstrips = [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,
             -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14]

STRIPMAP = {
    14: 0,
    13: 1,
    12: 2,
    11: 3,
    10: 4,
    9: 5,
    8: 6,
    7: 7,
    6: 8,
    5: 9,
    4: 10,
    3: 11,
    2: 12,
    1: 13,
    -1: 14,
    -2: 15,
    -3: 16,
    -4: 17,
    -5: 18,
    -6: 19,
    -7: 20,
    -8: 21,
    -9: 22,
    -10: 23,
    -11: 24,
    -12: 25,
    -13: 26,
    -14: 27,
}


def get_channel_config(configfile):
    cfg_channels = [config['DAQTOSTRIP'][i] for i in config['DAQTOSTRIP']]
    channels = [int(s[0] + s[1]) if len(s) >
                2 else None for s in [re.split('(\d+)', s) for s in cfg_channels]]
    return channels


def readout_status(channelconfig, strip):
    count = channelconfig.count(strip)
    if count == 2:
        return "DS"
    elif count == 1:
        return "SS"
    else:
        return None


print(f"Loaded config file: {filepath}")
