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

SOL = constants.getfloat("sol")
STRIPLENGTH = constants.getfloat("striplength")
STRIPVELOCITY = constants.getfloat("stripvelocity")
STRIPWIDTH = constants.getfloat("stripwidth")
STRIPSPACING = constants.getfloat("stripspacing")

MAXDELTA = 200 / (STRIPVELOCITY * SOL)  # Max LR offset in ns
MAXDELTAERROR = (PROPERROR / STRIPVELOCITY) * MAXDELTA

allstrips = [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, -
             1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14]


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
