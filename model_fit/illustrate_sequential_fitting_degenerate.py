'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  illustrates the effect of fitting loci sequentially;
          produces plots and likelihood surfaces to this effect.
          
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
sys.path.append('/Users/richard/Projects/FFPopSim/FFPopSim/pkg/python/')
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
ms=20
plt.rcParams.update(params)

def get_yticks(vals):
    n=len(vals)
    if max(vals)-min(vals) < 0.15:
        dyticks = round(max(vals)/5, 2)
    elif max(vals)-min(vals) < .2:
        dyticks = 0.05
    elif max(vals)-min(vals) < .4:
        dyticks = .1
    elif max(vals)-min(vals) < 1.0:
        dyticks = .2
    else:
        dyticks = .5
        
    scale_factor = 1.0*(n-1)/(max(vals)-min(vals))
    dytickvals = np.arange(0,max(vals),dyticks)
    dyticklocs = (dytickvals-min(vals))*scale_factor
    return dytickvals,dyticklocs

def get_xticks(vals):
    n=len(vals)
    if max(vals)-min(vals) < 15:
        dxticks = round((max(vals)-min(vals))/5)
    elif max(vals)-min(vals) < 40:
        dxticks = 10
    elif max(vals)-min(vals) < 70:
        dxticks = 20
    elif max(vals)-min(vals) < 100:
        dxticks = 25
    else:
        dxticks = 50
   
    scale_factor = 1.0*(n-1)/(max(vals)-min(vals))
    dxtickvals = np.arange(min(vals)-min(vals)%dxticks,max(vals),dxticks, dtype='int')
    dxticklocs = (dxtickvals-min(vals))*scale_factor
    
    if dxticklocs[0]<0:
        dxticklocs=dxticklocs[1:]
        dxtickvals=dxtickvals[1:]
    return dxtickvals, dxticklocs



plt.ion()

#set up parameters
N = default_params['N']                        # Population size
mu = default_params['mu']                        # Mutation rate
r = default_params['r'] 
f =default_params['f']        # Fitness (main/additive effects)
sample_size=default_params['sample_size'] 
tp = [0,15,80,150,250]
S=1
F=10
dynrange=20
L=len(f)

#generate test data
gt, gt_traj,pop= test_data(N,mu,r,f,tp,sample_size)
#initialize ctl_fit instance and perform initial guesses
ctl = ctl_fit.ctl_fit()
ctl.setup(N,mu,gt,tp,sample_size*np.ones_like(gt),F,S)
ctl.initial_guesses_one_parameter_ST()
uptolocus=L


##############################################################################
# PRODUCE FIGURE WITH LH SURFACES AFTER INTIAL FITTING
##############################################################################

fig=plt.figure(figsize = (8,6))
fig.subplots_adjust(left=0.12, bottom=0.14, top=0.9, right=0.9, hspace = 0.3)
ax=fig.add_subplot('111')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

#generate likelihood surface for initial guesses.
#note that the surfaces are based on a function (ML) other than the one used to generate the guesses.
#therefore, there is no reason to expect the guesses to be any good.
for locus in xrange(uptolocus):
    ds=np.round(ctl.initial_fitness[locus]/20,3)
    dtau=2
    frange = np.arange(ctl.fitness[locus]-20*ds, ctl.fitness[locus]+20*ds, ds)
    strange = np.arange(max(0,ctl.seedtimes[locus]-15*dtau), ctl.seedtimes[locus]+15*dtau, dtau)
    LH,minLH,s,t = ctl.LH_slice(locus, frange, strange, locus+1, False)
    if (locus==0):
        ctl.fitness[locus]=frange[np.argmin(LH[:,0])]
        ctl.seedtimes[locus]=0

    if locus>0:
        fig.add_subplot(2,2,locus)
        findex = np.where(frange==ctl.fitness[locus])
        sindex = np.where(strange==ctl.seedtimes[locus])
        #LH[findex, sindex] = minLH + 0.75*dynrange
        plt.title("Epitope "+str(locus+1)) 
        plt.imshow(np.ma.filled(LH,minLH+dynrange), interpolation='nearest',aspect='auto', vmin=minLH, vmax=minLH+dynrange)

        dytickvals, dyticklocs = get_yticks(frange)
        dxtickvals, dxticklocs = get_xticks(strange)
        
        plt.yticks(dyticklocs, map(str, dytickvals))
        plt.xticks(dxticklocs, map(str, dxtickvals))

ax.set_ylabel(r'escape rate ($\textrm{day}^{-1}$)')
ax.set_xlabel(r'seed times ($\textrm{day}$)')
plt.savefig('figures/sequential_degen_LH_initial.pdf')

#perform one round of multilocus fitting and save the resultant trajectory.
ctl.multi_locus_fit()
ctl.likelihood_c(uptolocus, 1, 1)


##############################################################################
# PRODUCE FIGURE WITH TRAJECTORIES AFTER MULTILCOUS FITTING INCLUDING LEGEND WITH DATA
##############################################################################
fig = plt.figure(figsize = (8,8))
margins = [0.18, 0.14, 0.77, 0.8]
axes = plt.Axes(fig, margins) 
fig.add_axes(axes)
plt.title(r'$N=10^{'+str(int(np.log10(N)))+'}$, $r='+str(r)+r'$, $\mu=10^{'+str(int(np.log10(mu)))+'}$')
model_traj=np.loadtxt('src/temp_test/traj.dat')
LH=np.loadtxt('src/temp_test/data.dat')
for locus in xrange(uptolocus+1):
    plt.plot(gt_traj[:,locus], ls='--', c=col[locus])
    plt.plot(tp, 1.0*gt[:,locus]/sample_size, ls='None',marker='o', c=col[locus], markersize=ms)    
    #if locus>0:  plt.plot(model_traj[:,locus], c=col[locus], label = r'$\epsilon_'+str(locus)+'='+str(round(ctl.fitness[locus-1],2))+r',\;\tau_'+str(locus)+'='+str(int(ctl.seedtimes[locus-1]))+'$')
    
    if locus==0:
        plt.plot(model_traj[:,locus], c=col[int(np.log2(locus+1))], label = 'Founder')
    #if locus>0:  plt.plot(model_traj[:,locus], c=col[locus], label = r'$'+''.join(['+']*locus)+''.join(['-']*(L-locus))+'$' )
    else: plt.plot(model_traj[:,locus], c=col[locus], label = r'$'+''.join(['+']*locus)+''.join(['-']*(L-locus))+'$' )
    #else:  plt.plot(model_traj[:,locus], c=col[locus])

ax=plt.gca()
fig.text(0.03,0.93,'A', fontsize=44)
plt.ylim([1.0/N, 1])
plt.subplots_adjust(bottom=0.15)
plt.xlabel(r'time ($\textrm{days}$)')
plt.ylabel('genotype frequencies')
plt.legend(loc=4)
plt.savefig('figures/sequential_traj.pdf')
ax.set_yscale('log')
plt.ylim([1.0/N, 2])
plt.yticks([1e-8, 1e-6, 1e-4, 1e-2, 1])
plt.show()
plt.savefig('figures/sequential_degen_traj_log.pdf')


#perform additional rounds of multilocus fitting
ctl.multi_locus_fit(n_iter=5)
plt.figure(figsize=(8,6))
plt.title(r'$N=10^{'+str(int(np.log10(N)))+'}$, $r='+str(r)+r'$, $\mu=10^{'+str(int(np.log10(mu)))+'}$')
model_traj=np.loadtxt('src/temp_test/traj.dat')
LH=np.loadtxt('src/temp_test/data.dat')
for locus in xrange(uptolocus+1):
    plt.plot(gt_traj[:,locus], ls='--', c=col[locus])
    plt.plot(tp, 1.0*gt[:,locus]/sample_size, ls='None',marker='o', c=col[locus], markersize=ms)    
    if locus>0:  plt.plot(model_traj[:,locus], c=col[locus], label = r'$\epsilon_'+str(locus)+'='+str(round(ctl.fitness[locus-1],2))+r',\;\tau_'+str(locus)+'='+str(int(ctl.seedtimes[locus-1]))+'$')
    else:  plt.plot(model_traj[:,locus], c=col[locus])

ax=plt.gca()
plt.ylim([1.0/N, 1])
plt.xlabel(r'time ($\textrm{days}$)')
plt.ylabel('genotype frequencies')
plt.subplots_adjust(bottom=0.15)
plt.show()
plt.legend(loc=4)
plt.savefig('figures/sequential_traj_refined.pdf')
plt.ylim([1.0/N, 2])
ax.set_yscale('log')
plt.savefig('figures/sequential_degen_traj_refined_log.pdf')

##############################################################################
# PRODUCE FIGURE WITH LH SURFACES AFTER MULTILCOUS FITTING
##############################################################################


fig=plt.figure(figsize = (8,6))
fig.subplots_adjust(left=0.12, bottom=0.14, top=0.9, right=0.9, hspace = 0.3)
ax=fig.add_subplot('111')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
for locus in xrange(uptolocus):
    ds=np.round(ctl.fitness[locus]/20,3)
    dtau=2
    
    frange = np.arange(ctl.fitness[locus]-20*ds, ctl.fitness[locus]+20*ds, ds)
    strange = np.arange(max(0,ctl.seedtimes[locus]-15*dtau), ctl.seedtimes[locus]+15*dtau, dtau)
    LH,minLH,s,t = ctl.LH_slice(locus, frange, strange, locus+1, False)

    if locus>0:
        fig.add_subplot(2,2,locus)
        findex = np.where(frange==ctl.fitness[locus])
        sindex = np.where(strange==ctl.seedtimes[locus])
        #LH[findex, sindex] = minLH + 0.75*dynrange
        plt.title("Epitope "+str(locus+1)) 
        plt.imshow(np.ma.filled(LH,minLH+dynrange), interpolation='nearest',aspect='auto', vmin=minLH, vmax=minLH+dynrange)
        dytickvals, dyticklocs = get_yticks(frange)
        dxtickvals, dxticklocs = get_xticks(strange)
        
        plt.yticks(dyticklocs, map(str, dytickvals))
        plt.xticks(dxticklocs, map(str, dxtickvals))

ax.set_ylabel(r'escape rate ($\textrm{day}^{-1}$)')
ax.set_xlabel(r'seed times ($\textrm{day}$)')

plt.savefig('figures/sequential_degen_LH_final_refined.pdf')

##############################################################################
# PRODUCE FIGURE WITH TRAJECTORIES AFTER MULTILCOUS FITTING
##############################################################################
plt.figure(figsize = (8,6))

LH=np.loadtxt('src/temp_test/data.dat')
model_traj=np.loadtxt('src/temp_test/traj.dat')
for locus in xrange(uptolocus+1):
    plt.plot(gt_traj[:,locus], ls='--', c=col[locus])
    plt.plot(tp, 1.0*gt[:,locus]/sample_size, ls='None',marker='o', c=col[locus], markersize=15)    
    if locus>0:  plt.plot(model_traj[:,locus], c=col[locus], label = r'$'+''.join(['+']*locus)+''.join(['-']*(L-locus))+'$' )
    else:  plt.plot(model_traj[:,locus], c=col[locus], label = 'Founder')
ax=plt.gca()
plt.xlabel(r'time ($\textrm{days}$)')
plt.ylabel('genotype frequencies')
plt.subplots_adjust(bottom=0.15)
plt.show()
plt.legend(loc=4)
ctl.ctl_clean()
plt.savefig('figures/sequential_degen_traj_illustration.pdf')
