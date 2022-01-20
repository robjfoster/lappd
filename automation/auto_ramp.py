import os
import sys
import subprocess

from auto_utils import hvparse, monitor, cmd_sethv, check_hv, ramp_and_check


init_hvinfo = monitor()
channels = [1,2,3,4,5]

detector = sys.argv[1]
target_voltage = sys.argv[2]

if detector == 104 or detector == 96:
    print(f"L{detector} selected.")
    if target_voltage >2450:
        sys.exit("Voltage too high for chosen LAPPD.")
else:
    sys.exit("The selected detector doesn't have ramp settings yet.")