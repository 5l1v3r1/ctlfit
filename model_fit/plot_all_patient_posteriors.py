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
taulist=[10,  20]
patients = ['CH40', 'CH58', 'CH77']
Flist = [2,5,10]
mu = 1e-5
for F in Flist:
    for pat in patients:
        for tau in taulist:
            p.plot_patient_posterior(Nlist, mu, tau, pat, F)
            plt.close('all')

