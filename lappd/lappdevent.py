import os
import pdb
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from scipy import signal

from utils.lappdcfg import config as lcfg
import utils.gimmedatwave as gdw
import utils.sigutils as su

root.gSystem.Load(os.path.dirname(os.path.realpath(__file__)) 
                  + "/utils/caenreader_cpp.so")
caenreader = root.CAENReader

base_dir = "/run/media/watchman/4d90ac40-7247-42ac-9d73-a76f7f121c5d/LAPPD/data/ledtest/125mm-L104-825VMCP-50VPC/"

# Attempt at event structure
# Each LAPPDEvent contains 28 StripEvents
# StripEvent contains the raw waveforms for each side of that strip
# Also contains Pulses, sliced waveforms around each identified peak, these
# pulses *should* be correlated

class LAPPDEvent():
    
    def __init__(self) -> None:
        pass

class StripEvent():
    
    def __init__(self) -> None:
        pass

class Pulse():
    
    def __init__(self, waveform, ns_per_sample=0.2) -> None:
        pass

# "An event" = two pulses with calculated offset and amplitude

def build_strip_event(stripnumber, event_no, dir):
    # Take stripnumber, find associated files from cfg
    # Load event from each channel
    # Preprocess each wave, storing raw wave?
    # Perform "event building" - find offset, centroid, amplitude
    # Multiple events may be created from a single file - create subevents?
    leftchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"L"]
    leftfile = gdw.get_filename(leftchannel, base=dir)
    leftwave = caenreader.readCAENWave(leftfile, event_no)
    rightchannel = lcfg['STRIPTODAQ'][str(stripnumber)+"R"]
    rightfile = gdw.get_filename(rightchannel, base=dir)
    rightwave = caenreader.readCAENWave(rightfile, event_no)
    caenreader.preprocessWave(leftwave.wave)
    caenreader.preprocessWave(rightwave.wave)
    leftpeaks = caenreader.findPeaks(leftwave.wave, -4, 15)
    rightpeaks = caenreader.findPeaks(rightwave.wave, -4, 15)
    peaks, offsets = su.coarse_correlate(leftwave, leftpeaks, rightwave, rightpeaks, 3)
    for peak in peaks:
        slicedleft = caenreader.sliceAroundPeak(leftwave, peak[0], 3, 3)
        pdb.set_trace()
    pdb.set_trace()

if __name__ == "__main__":
    build_strip_event(9, 280 , base_dir)