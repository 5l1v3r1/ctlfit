'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  simulate and fit data based on varying estimated mu and N and simulated r
'''

#!/usr/bin/python                                                                                                          
import os
import numpy as np
import sys
import subprocess as sp

maxjobs = 200
jobno = 0
#varied N
Nlist = [1e5, 1e6, 1e7, 1e8]
taulist = [0]
for N in Nlist:
    for tau in taulist:
        for jobno in range(maxjobs):
            call = 'qsub -p 50 -cwd submit/simple_wrap.py --N ' + str(N) + ' --tau ' + str(tau) + ' --runno ' + str(jobno) + ' --runname varN'
            call = call.split(' ')
            sp.call(call)

#varied mu
mulist = [1e-4, 1e-5, 1e-6]
taulist = [0]
Nlist = [1e7]
for mu in mulist:
    for tau in taulist:
        for jobno in range(maxjobs):
            call = 'qsub -p 50 -cwd submit/simple_wrap.py --mu ' + str(mu) + ' --tau ' + str(tau) + ' --runno ' + str(jobno) + ' --runname varmu_logN' + str(np.log10(N))
            call = call.split(' ')
            sp.call(call)


#varied r
rlist = [0, 1e-3, 1e-2, 1e-1]
taulist = [0]
Nlist = [1e7]
for r in rlist:
    for tau in taulist:
        for jobno in range(maxjobs):
            call = 'qsub -p 50 -cwd submit/simple_wrap.py --r ' + str(r) + ' --tau ' + str(tau) + ' --runno ' + str(jobno) + ' --runname varr_logN' + str(np.log10(N))
            call = call.split(' ')
            sp.call(call)



