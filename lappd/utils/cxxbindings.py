import os

import ROOT as root

root.gSystem.Load(os.path.dirname(os.path.realpath(__file__))
                  + "/../cxx/caenreader_cxx.so")
caenreader = root.CAENReader

# root.gSystem.Load(os.path.dirname(os.path.realpath(__file__))
#                  + "/../cxx/lognormal_cxx.so")
#lognormal = root.lognormal
