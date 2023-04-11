'''
File to automatically ramp up the 5 LAPPD voltages automatically through the CAEN V6533N module utilising the sethv command wrappers.
Useage: "python3 auto_ramp.py [detector] [target_mcp_voltage]"
e.g. "python3 auto_ramp.py 104 825" for ramping L104 up to 825 V across each MCP.

!!!!To kill all channels immediately, give -666 as the [target_mcp_voltage]!!!!
Usage: "python3 auto_ramp.py 104 -666"
'''

import sys
import subprocess

from auto_utils import monitor, cmd_sethv, ramp_and_check

init_hvinfo = monitor()


def init_hv_channels():
    '''Defines HV channel settings to match the CAEN crate.'''
    ch_XoX = 0
    ch_EoX = 1
    ch_XoE = 2
    ch_EoE = 3
    ch_PC = 4

    channels = [ch_XoX, ch_EoX, ch_XoE, ch_EoE, ch_PC]
    return channels


def hv_power(onoff):
    for i, chan in enumerate(hv_channels):
        if onoff:
            if init_hvinfo[chan]['power'] == 'ON':
                print('Channel ' + chan + ' power already on.')
            else:
                reset_cmd = f'/home/watchman/Documents/LAPPDCrate/sethv/sethv --power {chan} On'
                subprocess.call(reset_cmd, shell=True)
        else:
            if init_hvinfo[chan]['power'] == 'OFF':
                print('Channel ' + chan + ' power already off.')
            else:
                reset_cmd = f'/home/watchman/Documents/LAPPDCrate/sethv/sethv --power {chan} Off'
                subprocess.call(reset_cmd, shell=True)


hv_channels = init_hv_channels()
detector = int(sys.argv[1])
target_mcp_voltage = int(sys.argv[2])
mcp_to_anode_voltage = 200
mcp_to_mcp_voltage = 200


############## KILL SWITCH ###############
if target_mcp_voltage == -666:
    hv_power(0)
    sys.exit("All channels turned off")
##########################################


# Check voltage won't be too high for the Photocathode for the selected LAPPD.
v_max = mcp_to_anode_voltage + mcp_to_mcp_voltage + (target_mcp_voltage*2)

if detector in (104, 96):
    print(f"L{detector} selected.")
    if v_max > 2450:
        sys.exit("Voltage too high for chosen LAPPD.")
else:
    sys.exit("The selected detector doesn't have ramp settings yet.")


current_mcp_voltage = int(init_hvinfo[1]['vmon'] - init_hvinfo[0]['vmon'])
print("current mcp volts = "+str(current_mcp_voltage))


if abs(target_mcp_voltage - current_mcp_voltage) < 200:
    v_steps = [target_mcp_voltage]
elif target_mcp_voltage < current_mcp_voltage:
    v_steps = list(range(current_mcp_voltage, int(target_mcp_voltage), -200))
else:
    v_steps = list(range(current_mcp_voltage, target_mcp_voltage, 200))

if v_steps[0] == current_mcp_voltage:
    v_steps.pop(0)

if v_steps[-1] != target_mcp_voltage:
    v_steps.append(target_mcp_voltage)
    print("v_steps: ", v_steps)


# Set up the 5 voltages for the LAPPD. If over 200 V then ramp in stages.

hv_power(1)

for v_mcp in v_steps:
    if v_mcp in (0, 200):
        v_XoX = v_EoX = v_XoE = v_EoE = v_PC = v_mcp
        print(f"Setting all voltages to {v_mcp}.")
    else:
        v_XoX = mcp_to_anode_voltage
        v_EoX = mcp_to_anode_voltage + v_mcp
        v_XoE = mcp_to_anode_voltage + v_mcp + mcp_to_mcp_voltage
        v_EoE = mcp_to_anode_voltage + mcp_to_mcp_voltage + (2 * v_mcp)
        v_PC = mcp_to_anode_voltage + mcp_to_mcp_voltage + (2 * v_mcp)
    chan_volts = [v_XoX, v_EoX, v_XoE, v_EoE, v_PC]

    print(f"Setting voltages for {v_mcp} across the MCPs.")

    for i, chan in enumerate(hv_channels):
        cmd_sethv(chan, chan_volts[i])

    for i, chan in enumerate(hv_channels[::-1]):
        ramp_and_check(chan)

# Check voltages after ramping is complete. If voltages are all 0, turn off power.
final_hvinfo = monitor()

if final_hvinfo[4]['vmon'] < 2:
    print("Voltage is at zero on all channels. Turning power off.")
    hv_power(0)

else:
    print("Voltage has been reached on all channels. The detector is ready")
