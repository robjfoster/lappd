import os
import sys

time = sys.argv[1]
config_loc = sys.argv[2]

shellstring = ""
shellstring += '(sleep 3s && echo "s" && sleep 3s && echo "W" && '
shellstring += f'sleep {time}s && echo "W" && sleep 3s && echo "s" && sleep 3s && echo "q")'
shellstring += f' | wavedump {config_loc}'

with open(f"rundark_{time}.sh", "w") as f:
    f.write(shellstring)
