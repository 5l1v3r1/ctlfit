'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  illustrates the concept of "dominant genotypes" as well as the lesser genotypes.
          Log and linear scale plots are produced.
'''

import numpy as np
import pylab as plt
import sys
sys.path.append('src/')
import ctl_fit
from ctlutils import test_data
from ctlutils import col
from ctlutils import params as default_params
import os
import FFPopSim as h
params = {'backend': 'ps',  
          'axes.labelsize': 36, 
          'axes.titlesize': 36,
          'text.fontsize': 36,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 25,
'xtick.labelsize': 30,
'ytick.labelsize': 30,
'lines.linewidth': 3,
'text.usetex': True}
plt.rcParams.update(params)

plt.ion()


#set up parameters
S=1
F=10

N = default_params['N']                        # Population size
mu = default_params['mu']                        # Mutation rate
r = default_params['r'] 
f =default_params['f']        # Fitness (main/additive effects)
sample_size=default_params['sample_size'] 
tp = default_params['tp'] 

dynrange=20
L=len(default_params['f'])

#generate test data
gt, gt_traj, pop,all_gt_traj = test_data(N,mu,r,f,tp,sample_size, True)

fig = plt.figure(figsize = (8,8))
margins = [0.18, 0.14, 0.77, 0.8]
axes = plt.Axes(fig, margins) 
fig.add_axes(axes)
#plot test data
considered_loci = [2**x-1 for x in range(L+1)]
for locus in xrange(2**L):
    if locus in considered_loci:
        if locus==0:
            plt.plot(all_gt_traj[:,locus], c=col[int(np.log2(locus+1))], label = 'Founder')
        else:
            plt.plot(all_gt_traj[:,locus], c=col[int(np.log2(locus+1))], label = r'$'+''.join(['+']*int(np.log2(locus+1)))+''.join(['-']*(L-int(np.log2(locus+1))))+'$')
        plt.plot(default_params['tp'], 1.0*gt[:,int(np.log2(locus+1))]/sample_size, ls='None',marker='o', c=col[int(np.log2(locus+1))], markersize=20)
    else:
        plt.plot(all_gt_traj[:,locus], ls='-', c='y')    

#adjust plot
plt.xlabel(r'time ($\textrm{days}$)')
plt.ylabel('genotype frequencies')
#plt.title(r'$N=10^{'+str(int(np.log10(N)))+'}$, $r='+str(r)+r'$, $\mu=10^{'+str(int(np.log10(mu)))+'}$')
fig.text(0.03,0.93,'A', fontsize=48)
axes.set_xlim([0, tp[-1]*1.05])
axes.set_ylim([0, 1.05])
#plt.legend(loc=4)
plt.show()
plt.savefig('figures/dominant_genotypes_demonstration.pdf')
plt.savefig('figures/dominant_genotypes_demonstration.svg')
axes.set_yscale('log')
axes.set_ylim([1.0/N*.5, 2])
axes.legend().set_visible(False)
del fig.texts[-1]
fig.text(0.03,0.93,'B', fontsize=48)
plt.yticks([1e-6, 1e-4, 1e-2, 1])
plt.savefig('figures/dominant_genotypes_demonstration_log.pdf')
plt.savefig('figures/dominant_genotypes_demonstration_log.svg')
plt.show()
