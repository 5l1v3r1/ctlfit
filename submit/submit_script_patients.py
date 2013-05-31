#!/usr/bin/python                                                                                                          
'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  wrapper for performing posterior estimating on the cluster
'''

import os
import sys
JOBDIR="/ebio/ag-neher/share/users/rneher/CTL_Fitting/patient_fit/"
command = "/ebio/ag-neher/share/programs/EPD/bin/python "+ JOBDIR+"fit_patients_posterior_only.py "
for arg in sys.argv[1:]:
    command+=" "+arg

print command
os.system(command)
