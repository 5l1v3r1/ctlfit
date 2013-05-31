#!/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/bin/python
'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  just a simple wrapper for fitting from simulated data
'''
import os
import sys

JOBDIR="/ebio/ag-neher/share/users/rneher/CTL_Fitting/model_fit/"
command = "/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/bin/python "+ JOBDIR+"simple_fit.py "
for arg in sys.argv[1:]:
    command=command + " " + arg

print command
os.system(command)
