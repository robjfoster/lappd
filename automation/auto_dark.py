import pdb
import os
import sys
import subprocess

from auto_utils import hvparse, monitor, cmd_wavedump, cmd_sethv, check_hv, ramp_and_check


time = sys.argv[1]
try:
    config_path = sys.argv[2]
except IndexError:
    config_path = "/run/media/watchman/4d90ac40-7247-42ac-9d73-a76f7f121c5d/LAPPD/LAPPDConfig_X742.txt"

PC_CHANNEL = 4
#MONITOR_CMD = '/home/watchman/Documents/LAPPD/CAENCrate/sethv/sethv --monitor'
startdir = os.getcwd()

class cd:
    # context manager for cd command
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


init_hvinfo = monitor()
ORIGINAL_PC_VOLTAGE = init_hvinfo[PC_CHANNEL]['vset']
if ORIGINAL_PC_VOLTAGE > 2055: sys.exit("Has the voltage been reset?")

voltages = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250]

# Main command loop
for setv in voltages:
    dirname = f"pc_{setv}"
    dirpath = os.path.join(startdir, dirname)
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
    with cd(dirpath):
#    subprocess.call(f"mkdir -p pc_{setv}", shell=True)
#    subprocess.call(f"cd pc_{setv}", shell=True)
        print("\n\n\n\n\n\n\n\n\n\n")
        pc_voltage = ORIGINAL_PC_VOLTAGE + setv
        print(f"###\n###\n###\n Setting pc_voltage to {pc_voltage}")
        cmd_sethv(PC_CHANNEL, pc_voltage)
        print("Voltage set")
        ramp_and_check(PC_CHANNEL)
        print("Waiting for voltage to stabilise...")
        subprocess.call("sleep 5s", shell=True)
        cmd_wavedump(time, config_path)
        subprocess.call("cd ..", shell=True)
# Now reset voltages back to beginning
cmd_sethv(PC_CHANNEL, ORIGINAL_PC_VOLTAGE)
ramp_and_check(PC_CHANNEL)
print("PC voltage reset")
print("Done.")

# wavedump_cmd = ""
# wavedump_cmd += '(sleep 3s && echo "s" && sleep 3s && echo "W" && '
# wavedump_cmd += f'sleep {time}s && echo "W" && sleep 3s && echo "s" && sleep 3s && echo "q")'
# wavedump_cmd += f' | wavedump {config_path}'

# for v in range(0, 101, 10):
#     # maybe change this
#     print("\n\n\n\n\n\n\n\n\n\n")
#     pc_voltage = ORIGINAL_PC_VOLTAGE + v
#     print(f"###\n###\n###\n Setting pc_voltage to {pc_voltage}")
#     vset_cmd = f'/home/watchman/Documents/LAPPD/CAENCrate/sethv/sethv --VSet {PC_CHANNEL} {pc_voltage}'
#     subprocess.call(vset_cmd, shell=True)   
#     print("Voltage set")
#     # monitor stuff
#     ramp_and_check()
#     # loop_count = 0
#     # while True:
#     #     if check_hv():
#     #         print("Voltage reached")
#     #         break
#     #     else:
#     #         print("Ramping voltage...")
#     #         loop_count += 1
#     #         if loop_count > 15: sys.exit("Monitoring took a long time, something's wrong")
#     #         subprocess.call("sleep 10s", shell=True)
#     print("Waiting for voltage to stabilise...")
#     subprocess.call("sleep 5s", shell=True)
#     subprocess.call(wavedump_cmd, shell=True)
# # Power down now
# powerdown_cmd = f'/home/watchman/Documents/LAPPD/CAENCrate/sethv/sethv --VSet {PC_CHANNEL} {ORIGINAL_PC_VOLTAGE}'
# subprocess.call(powerdown_cmd, shell=True)
# ramp_and_check()
# # loop_count = 0
# # while True:
# #         if check_hv():
# #             print("Voltage reached")
# #             break
# #         else:
# #             print("Ramping voltage...")
# #             loop_count += 1
# #             if loop_count > 15: sys.exit("Monitoring took a long time, something's wrong")
# #             subprocess.call("sleep 10s", shell=True)
# print("PC voltage reset")
# print("Done.")
