import configparser
import sys

config = configparser.ConfigParser()
config.read("config/L104setup.cfg")
if not config.sections():
    sys.exit("Config file is empty (probably wrong file path)")
