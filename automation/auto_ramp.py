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

# Check voltage won't be too high for the Photocathode for the selected LAPPD.
v_max = 400+(target_mcp_voltage*2)

if detector in (104, 96):
    print(f"L{detector} selected.")
    if v_max > 2450:
        sys.exit("Voltage too high for chosen LAPPD.")
else:
    sys.exit("The selected detector doesn't have ramp settings yet.")


v_steps = list(range(200, target_mcp_voltage, 200))
print("v_steps:", v_steps)
if v_steps[-1] < target_mcp_voltage:
    v_steps.append(target_mcp_voltage)

# Set up the 5 voltages for the LAPPD. If over 200 V then ramp in stages.
# TODO: Set it to ramp back down in stages too.
# Currently can just set all to 0 V and ramp from there.

hv_power(1)

for v_mcp in v_steps:
    if v_mcp in (0, 200):
        v_XoX = v_EoX = v_XoE = v_EoE = v_PC = v_mcp
        print(f"Setting all voltages to {v_mcp}.")
    else:
        v_XoX = 200
        v_EoX = 200 + v_mcp
        v_XoE = 400 + v_mcp
        v_EoE = 400 + (2 * v_mcp)
        v_PC = 400 + (2 * v_mcp)
    chan_volts = [v_XoX, v_EoX, v_XoE, v_EoE, v_PC]

    print(f"Setting voltages for {v_mcp} across the MCPs.")

    for i, chan in enumerate(hv_channels):
        cmd_sethv(chan, chan_volts[i])

    for i, chan in enumerate(hv_channels[::-1]):
        ramp_and_check(chan)

    print("Voltage has been reached on all channels. The detector is ready")
