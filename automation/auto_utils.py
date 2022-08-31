import sys
import subprocess


def hvparse(output: str):
    '''Parsing information from CAEN HV module.'''
    output = output.decode("utf-8").split("\n")
    error_msg = output[0]
    header = output[1].split(",")
    hvinfo = {}
    for row in output[2:-1]:
        row = row.split(",")
        channel = int(row[1])
        values = {}
        values["model"] = row[0]
        values["iset"] = float(row[2].split("uA")[0])
        values['vset'] = float(row[3].split("V")[0])
        values['vmon'] = float(row[4].split("V")[0])
        values['vmax'] = float(row[5].split("V")[0])
        values["imon"] = float(row[6].split("uA")[0])
        values["rup"] = float(row[7].split("Vps")[0])
        values["rdown"] = float(row[8].split("Vps")[0])
        values["power"] = row[9]
        values["err"] = row[10]
        hvinfo[channel] = values
    return hvinfo

def monitor(cmd='/home/watchman/Documents/LAPPD/CAENCrate/sethv/sethv --monitor'):
    '''Getting monitoring infomation from HV.'''
    # TODO: Rampup/Rampdown is considered an error message
    try:
        output = subprocess.check_output(cmd, shell=True)
        hvinfo = hvparse(output)
    except subprocess.CalledProcessError as e:
        print(f"Error message: {e.output}")
        sys.exit("Monitoring error")
    return hvinfo

def cmd_wavedump(time, config_path):
    '''Running wavedump.'''
    wavedump_cmd = ""
    wavedump_cmd += '(sleep 3s && echo "s" && sleep 3s && echo "W" && '
    wavedump_cmd += f'sleep {time}s && echo "W" && sleep 3s && echo "s" && sleep 3s && echo "q")'
    wavedump_cmd += f' | wavedump {config_path}'
    subprocess.call(wavedump_cmd, shell=True)

def cmd_sethv(channel, voltage):
    '''Setting HV to desired voltage.'''
    reset_cmd = f'/home/watchman/Documents/LAPPD/CAENCrate/sethv/sethv --VSet {channel} {voltage}'
    subprocess.call(reset_cmd, shell=True)

def check_hv(channel: int, tolerance:float = 1) -> bool:
    '''Check vset and vmon are with tolerance'''
    # try:
    #     monitoroutput = subprocess.check_output(MONITOR_CMD, shell=True).decode("utf-8")
    # except subprocess.CalledProcessError as e:
    #     print(e.output)
    #     sys.exit("dead")
    hvinfo = monitor()
    vset = hvinfo[channel]['vset']
    vmon = hvinfo[channel]['vmon']
    # parsed_output = monitoroutput.split("\n")
    # header = parsed_output[:2]
    # hv_data = parsed_output[2:]
    # hv_channel = hv_data[PC_CHANNEL]
    # hv_channel_data = hv_channel.split(",")
    # if int(hv_channel_data[1]) != PC_CHANNEL: sys.exit("WRONG CHANNEL")
    # vmon = float(hv_channel_data[4].split("V")[0])
    return(bool(abs(vmon-vset) <= tolerance))


def ramp_time(channel:int, tolerance:float = 1) -> float:
    '''Estimate time needed to ramp to target voltage'''
    hvinfo = monitor()
    vmon = hvinfo[channel]['vmon']
    vset = hvinfo[channel]['vset']

    if vset < vmon-tolerance:
        rampupdown = 'rup'
    elif vset > vmon+tolerance:
        rampupdown = 'rdown'
    elif abs(vset-vmon) <= tolerance:
        ramptime = 0
        return ramptime

    ramprate = hvinfo[channel][rampupdown]
    vdiff = abs(vset-vmon)

    ramptime = vdiff/ramprate
    return ramptime


def ramp_and_check(channel:int, tolerance:float = 1, sleeptime=None, loopcount=15):
    '''Ramp voltages and monitor along the way.'''
    loop_count = 0
    ramptime = ramp_time(channel, tolerance) + 1
    stime = ramptime if not sleeptime else sleeptime

    while True:
        if check_hv(channel):
            print("Voltage reached")
            break
        print("Ramping voltage...")
        loop_count += 1
        if loop_count > loopcount:
            sys.exit("Monitoring took a long time, something's wrong")
        subprocess.call(f"sleep {stime}s", shell=True)
