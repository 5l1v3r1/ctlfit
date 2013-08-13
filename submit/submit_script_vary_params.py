#!/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/bin/python
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

maxjobs = 100
jobno = 0

simMu = 1e-5

#varied N
Nlist = [1e4, 1e5, 1e6, 1e7, 1e8]
taulist = [0]
mu=simMu
for N in Nlist:
    for tau in taulist:
        for jobno in range(maxjobs):
            call = 'qsub  -l h_rt=00:59:59  -cwd submit/simple_wrap.py --N ' + str(N) +' --simMu '+str(simMu)+' --mu '+str(mu) + ' --tau ' + str(tau) + ' --runno ' + str(jobno) + ' --runname varN'
            os.system(call)
            #call = call.split(' ')
            #sp.call(call)

#varied mu
mulist = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
taulist = [0]
Nlist = [1e7]
for N in Nlist:
    for mu in mulist:
        for tau in taulist:
            for jobno in range(maxjobs):
                call = 'qsub -l h_rt=00:59:59  -cwd submit/simple_wrap.py --simN ' + str(N) +' --simMu '+str(simMu)+ ' --N ' + str(N) + ' --mu ' + str(mu) + ' --tau ' + str(tau) + ' --runno ' + str(jobno) + ' --runname varmu_logN' + str(np.log10(N))
                os.system(call)
#                call = call.split(' ')
#                sp.call(call)


#varied r
rlist = [0, 1e-3, 1e-2, 1e-1]
taulist = [0]
Nlist = [1e7]
mu=simMu
for N in Nlist:
    for r in rlist:
        for tau in taulist:
            for jobno in range(maxjobs):
                call = 'qsub -l h_rt=00:59:59 -cwd submit/simple_wrap.py --simN ' + str(N) +' --simMu '+str(simMu)+' --mu '+str(mu)+ ' --N ' + str(N) + ' --r ' + str(r) + ' --tau ' + str(tau) + ' --runno ' + str(jobno) + ' --runname varr_logN' + str(np.log10(N))
                os.system(call)
#                call = call.split(' ')
#               sp.call(call)

