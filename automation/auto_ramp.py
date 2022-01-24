import os
import sys
import subprocess

from auto_utils import hvparse, monitor, cmd_sethv, check_hv, ramp_and_check

init_hvinfo = monitor()

def init_hv_channels():
    ch_XoX = 0
    ch_EoX = 1
    ch_XoE = 2
    ch_EoE = 3
    ch_PC = 4

    channels = [ch_XoX, ch_EoX, ch_XoE, ch_EoE, ch_PC]
    return channels
 
    
hv_channels = init_hv_channels()
detector = sys.argv[1]
target_mcp_voltage = sys.argv[2]

v_max = 400*(target_mcp_voltage*2)

if detector == 104 or detector == 96:
    print(f"L{detector} selected.")
    if v_max >2450:
        sys.exit("Voltage too high for chosen LAPPD.")
else:
    sys.exit("The selected detector doesn't have ramp settings yet.")


v_steps = list(range(200,target_mcp_voltage, 200)
if v_steps[-1] < target_mcp_voltage: 
    v_steps.append(target_mcp_voltage)

for v_mcp in v_steps:
    if v_mcp == 0 or v_mcp == 200:
        v_XoX = v_EoX = v_XoE = v_EoE = v_PC = v_mcp
    else:
        v_XoX = 200
        v_EoX = 200 + v_mcp
        v_XoE = 400 + v_mcp
        v_EoE = 400 + (2 * v_mcp)
        v_PC = 400 + (2 * v_mcp)
    chan_volts = [v_XoX, v_EoX, v_XoE, v_EoE, v_PC]
    
    for i,chan in enumerate(hv_channels):
        cmd_sethv(hv_channels[i], chan_volts[i])
    
    for i,chan in enumerate(hv_channels[::-1]):
        ramp_and_check[chan]          
             
    


    
    
    
    
