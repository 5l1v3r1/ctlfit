'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  batch submit script to compute posterior distributions for patients.
'''

import pylab as plt
import sys
import os


Nlist = [1e5, 1e6,1e7,1e8]
taulist=[0,20]
patients = ['CH40', 'CH58', 'CH77']
Flist = [1,5,10]

jobno = 0

for F in Flist:
    for pat in patients:
        for N in Nlist:
            for tau in taulist:
                command='qsub -l h_rt=8:55:0 -cwd submit/submit_script_patients.py --pop '+str(N)+' --tau '+str(tau)+' --pat '+pat+' --F '+str(F) +' --runno ' + str(jobno)
                print command
                os.system(command)
