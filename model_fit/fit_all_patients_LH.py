'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  wrapper script for fitting all patients locally, with a variety of parameter estimates.
'''

import pylab as plt
import numpy as np
import sys
import subprocess as sp
import fit_patients as p

p = reload(p)

Nlist = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]
taulist=[0, 20]
patients = ['CH40', 'CH58', 'CH77']
Flist = [2,5,10]

for N in Nlist:
    for F in Flist:
        for pat in patients:
            for tau in taulist:
                p.fit_patient_LH(N, 1e-5, tau, pat, F)
                plt.close('all')
