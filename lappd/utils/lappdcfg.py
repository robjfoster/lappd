import configparser
import re
import sys

config = configparser.ConfigParser()
config.read("config/L104setup.cfg")
if not config.sections():
    sys.exit("Config file is empty (probably wrong file path)")

def get_channel_config(configfile):
    cfg_channels = [config['DAQTOSTRIP'][i] for i in config['DAQTOSTRIP']]
    channels =  [int(s[0]+ s[1]) if len(s)>2 else None for s in [re.split('(\d+)',s) for s in cfg_channels]]
    return channels

def readout_status(channelconfig, strip):
    count = channelconfig.count(strip)
    if count == 2:
        return "DS"
    elif count == 1:
        return "SS"
    else:
        return None