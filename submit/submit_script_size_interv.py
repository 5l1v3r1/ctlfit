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

#varied size
sizelist = [5, 10,20,50,100]
taulist = [0]
Nlist = [1e7]
simMu=1e-5
mu=simMu
F=10.0
T=''
for N in Nlist:
    for size in sizelist:
        for tau in taulist:
            for jobno in range(maxjobs):
                call = 'qsub  -l h_rt=00:10:59 -cwd submit/simple_wrap.py --size ' + str(size) + ' --simN ' + str(N) + ' --N ' + str(N) +' --simMu '+str(simMu)+' --mu '+str(mu)+\
                    ' --tau ' + str(tau) +' --F ' + str(F) + ' --runno ' + str(jobno) + ' --runname varsize_'+T+'logN' + str(np.log10(N)) + '_size_' + str(size)
                print call
                os.system(call)
#                call = call.split(' ')
#                sp.call(call)

#varied frequency
freqlist = [5, 10, 20, 40, 70,100]
taulist = [0]
Nlist = [1e7]
ssize = 20
for N in Nlist:
    for freq in freqlist:
        for tau in taulist:
            for jobno in range(maxjobs):
                call = 'qsub -l h_rt=00:10:59 -cwd submit/simple_wrap.py --interv ' + str(freq) + ' --simN ' + str(N) + ' --N ' + str(N) +' --simMu '+str(simMu)+' --mu '+str(mu)+\
                    ' --tau ' + str(tau)+' --F ' + str(F) + ' --runno ' + str(jobno) + ' --runname varfreq_'+T+'logN' + str(np.log10(N)) + '_freq_' + str(freq) + '_size_'+str(ssize)
                os.system(call)
#                call = call.split(' ')
#                sp.call(call)
