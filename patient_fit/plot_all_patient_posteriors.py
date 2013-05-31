'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  wrapper script to plot patient posteriors for range of F, tau, and N values.
'''

import pylab as plt
import numpy as np
import sys
import subprocess as sp
import plot_patient_posterior as p

p = reload(p)

Nlist = [1e5,1e6,1e7,1e8]
taulist=[0, 20]
patients = ['CH40', 'CH58', 'CH77']
Flist = [1, 5, 10]

for F in Flist:
    for pat in patients:
        for tau in taulist:
            p.plot_patient_posterior(Nlist, 1e-5, tau, pat, F)
            plt.close('all')

