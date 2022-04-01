import os

import ROOT as root

root.gSystem.Load(os.path.dirname(os.path.realpath(__file__))
                  + "/caenreader_cpp.so")
caenreader = root.CAENReader
