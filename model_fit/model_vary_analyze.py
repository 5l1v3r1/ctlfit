'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Produces plots based on varying N, mu, and r.
          Requires that submits/submit_script_vary_params.py has already been run.
          
'''

import numpy as np
import matplotlib.pyplot as plt
from ctlutils import params as default_params
from ctlutils import col
import sys
import copy
import glob
import cPickle as pickle
from time import strftime

params = {'backend': 'ps',  
          'axes.labelsize': 24, 
          'axes.titlesize': 24,
          'text.fontsize': 24,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 20,
'xtick.labelsize': 20,
'ytick.labelsize': 20,
'lines.linewidth': 2,
'text.usetex': True}
plt.rcParams.update(params)

plt.ion()

#initialize global parameters
tau = 0
F = 10.0
S = 1
today = '130418'

f=default_params['f']

#varied mu; initialize parameters
runname = 'varmu'
mulist = [1e-4, 1e-5, 1e-6]
N=default_params['N']
simN=default_params['N']
r=default_params['r']
numruns = 100

#LOAD DATA
fitnesses_mu = np.zeros([len(mulist), numruns, len(f)])
dirname = 'figures/' + today + '_varmu_logN'+ str(np.log10(N))  + '_tau_' + str(tau) + '_prior_F_'+str(F)+'_S_'+str(S)
for i in range(len(mulist)):
    for k in range(numruns):
        fitnesses_mu[i,k,:] = pickle.load(open(dirname + '/logN_' + str(int(np.log10(N))) + '_logmu_' + str(int(np.log10(mulist[i]))) + '_r_' + str(r) + '_runno_' + str(k) +'.pickle', 'r'))[0]


runname = 'varN'
mu=default_params['mu']
r=default_params['r']
Nlist = [1e5, 1e6, 1e7,1e8]
#LOAD DATA
fitnesses_N = np.zeros([len(Nlist), numruns, 5])
dirname = 'figures/' + today + '_varN_tau_' + str(tau) + '_prior_F_'+str(F)+'_S_'+str(S)
for j in range(len(Nlist)):
    for k in range(numruns):
        fitnesses_N[j,k,:] = pickle.load(open(dirname + '/logN_' + str(int(np.log10(Nlist[j]))) + '_logmu_' + str(int(np.log10(mu))) + '_r_' + str(r) + '_runno_' + str(k) +'.pickle', 'r'))[0]


#LOAD DATA
runname = 'varr'
mu=default_params['mu']
N=default_params['N']
rlist = [0.0, 0.001, 0.01, 0.1]

fitnesses_r = np.zeros([len(rlist), numruns, 5])
dirname = 'figures/' + today + '_varr_logN'+ str(np.log10(N))+'_tau_' + str(tau) + '_prior_F_'+str(F)+'_S_'+str(S)
for j in range(len(rlist)):
    for k in range(numruns):
        fitnesses_r[j,k,:] = pickle.load(open(dirname + '/logN_' + str(int(np.log10(N))) + '_logmu_' + str(int(np.log10(mu))) + '_r_' + str(rlist[j]) + '_runno_' + str(k) +'.pickle', 'r'))[0]


###################################################################
#MAKE COMBINED PLOT PN VARIATON WITH WRONG ASSUMPTIONS ABOUT THE MODEL PARAMETER
###################################################################

fig=plt.figure(figsize=(16,6))
fig.subplots_adjust(left=0.12, bottom=0.15, top=0.95, right=0.95, wspace = 0.02)

#PLOT MEANS AND STANDARD DEVIATIONS_ MUTATION RATE

plt.subplot(1,3,2)
for l in range(len(f)):
    x = mulist
    plt.errorbar(np.log10(np.array(x)), np.average(fitnesses_mu[:,:,l], axis=1)/f[l], 
                 np.std(fitnesses_mu[:,:,l], axis=1)/f[l], 
                 label='$\epsilon_' + str(l+1) + ' = ' + str(f[l]) + '$', elinewidth=2, capsize=7, mew=2, c=col[l+1])

ax=plt.gca()
ax.set_xlim([np.log10(min(x))-0.5, np.log10(max(x))+0.5])
ax.set_ylim([0.0,2])
ax.set_xlabel('model mutation rate $(\mathrm{days}^{-1})$')
ax.set_yticks([])
ticks = ['$10^{' + str(int(np.log10(z))) + '}$' for z in mulist]
ax.set_xticklabels(ticks)
xloc = maxloc = plt.MultipleLocator(1)
ax.xaxis.set_major_locator(xloc)
plt.axhline(y=1, color='k', ls='--')

#PLOT MEANS AND STANDARD DEVIATIONS_ POPULATION SIZE
plt.subplot(1,3,1)
for l in range(len(f)):
    x = Nlist
    plt.errorbar(np.log10(np.array(x)), np.average(fitnesses_N[:,:,l], axis=1)/f[l], 
                 np.std(fitnesses_N[:,:,l], axis=1)/f[l], 
                 label='$\epsilon_' + str(l+1) + ' = ' + str(f[l]) + '$', elinewidth=2, capsize=7, mew=2, c=col[l+1])

ax=plt.gca()
ax.set_xlim([np.log10(min(x))-0.5, np.log10(max(x))+0.5])
ax.set_ylim([0.0,2])
ax.set_xlabel('model population size')
ax.set_ylabel('normalized escape rate')
ticks = ['$10^{' + str(int(np.log10(z))) + '}$' for z in Nlist]
ax.set_xticklabels(ticks)
xloc = maxloc = plt.MultipleLocator(1)
ax.xaxis.set_major_locator(xloc)
plt.legend(loc=3,ncol=2, columnspacing=0.2, handletextpad=0.2)
plt.axhline(y=1, color='k', ls='--')

#PLOT MEANS AND STANDARD DEVIATIONS_ RECOMBINATION RATE
plt.subplot(1,3,3)
for l in range(len(f)):
    x = copy.deepcopy(rlist)
    x[0]+=1e-4
    plt.errorbar(np.log10(np.array(x)), np.average(fitnesses_r[:,:,l], axis=1)/f[l], 
                 np.std(fitnesses_r[:,:,l], axis=1)/f[l], 
                 label='$\epsilon_' + str(l+1) + ' = ' + str(f[l]) + '$', elinewidth=2, capsize=7, mew=2, c=col[l+1])

ax=plt.gca()
ax.set_xlim([np.log10(min(x))-0.5, np.log10(max(x))+0.5])
ax.set_ylim([0.0,2])
ax.set_xlabel('recombination rate $(\mathrm{days}^{-1})$')
ticks = ['$0$']
ticks += ['$10^{' + str(int(np.log10(z))) + '}$' for z in rlist[1:]]
ax.set_xticklabels(ticks)
ax.set_yticks([])
xloc = maxloc = plt.MultipleLocator(1)
ax.xaxis.set_major_locator(xloc)
plt.axhline(y=1, color='k', ls='--')

fig.text(0.13,0.86,'A', fontsize=32)
fig.text(0.42,0.86,'B', fontsize=32)
fig.text(0.7,0.86,'C', fontsize=32)

plt.show()
plt.savefig('figures/model_variation_logN_' + str(int(np.log10(N))) + '_tau_'+str(tau)+'_F_'+str(F)+'_S_'+str(S)+'.pdf')
plt.savefig('figures/model_variation_logN_' + str(int(np.log10(N))) + '_tau_'+str(tau)+'_F_'+str(F)+'_S_'+str(S)+'.svg')


#PLOT HISTOGRAMS
bins = np.linspace(0,.8,81)
for i,mu in enumerate(mulist):
    plt.figure()
    plt.title('mutation rate $10^{' + str(int(np.log10(mu))) + '}$, $N = 10^' + str(int(np.log10(N))) + '$, $S = ' + str(S) + '$, $F = ' + str(F) + '$')
    for l in xrange(len(f)):
        y,x=np.histogram(fitnesses_mu[i,:,l], bins, normed='True')
        xbins = 0.5*(x[1:]+x[:-1])
        plt.plot(xbins,y, ls='-', c=col[l+1])
        plt.plot([f[l],f[l]], [0,15], lw=2, c=col[l+1])
    ax=plt.gca()
    ax.set_xlim([0,.8])
    plt.xlabel('escape rate')
    plt.ylabel('frequency')
    plt.savefig('figures/model_variation_mu_logN_' + str(int(np.log10(N))) + '_logmu_'+str(int(np.log10(mu)))+'_histogram.pdf') 
    plt.savefig('figures/model_variation_mu_logN_' + str(int(np.log10(N))) + '_logmu_'+str(int(np.log10(mu)))+'_histogram.svg')
    plt.close()


###################################################################
#PLOT VARIATON WITH WRONG ASSUMPTIONS ABOUT THE POPULATION SIZE
###################################################################

bins = np.linspace(0,.6,51)
for i,N in enumerate(Nlist):
    plt.figure()
    plt.title('est. pop. size $10^{' + str(int(np.log10(N))) + '}$, $N = 10^{' + str(int(np.log10(simN))) + '}$, $S = ' + str(S) + '$, $F = ' + str(F) + '$')
    for l in xrange(len(f)):
        y,x=np.histogram(fitnesses_N[i,:,l], bins, normed='True')
        xbins = 0.5*(x[1:]+x[:-1])
        plt.plot(xbins,y, ls='-', c=col[l+1])
        plt.plot([f[l],f[l]], [0,25], lw=2, c=col[l+1])
    ax=plt.gca()
    ax.set_xlim([0,.6])
    plt.xlabel('escape rate')
    plt.ylabel('frequency')
    plt.savefig('figures/model_variation_N_logN_6_modelN_'+str(int(np.log10(N)))+'_histogram.pdf')



###################################################################
#PLOT VARIATON WITH WRONG ASSUMPTIONS ABOUT THE POPULATION SIZE
###################################################################

bins = np.linspace(0,.6,51)
for i,r in enumerate(rlist):
    plt.figure()
    plt.title('recombination rate $' + str(rlist[i]) + '$, $N = 10^' + str(int(np.log10(N))) + '$, $S = ' + str(S) + '$, $F = ' + str(F) + '$')
    for l in xrange(len(f)):
        y,x=np.histogram(fitnesses_r[i,:,l], bins, normed='True')
        xbins = 0.5*(x[1:]+x[:-1])
        plt.plot(xbins,y, ls='-', c=col[l+1])
        plt.plot([f[l],f[l]], [0,25], lw=2, c=col[l+1])
    ax=plt.gca()
    ax.set_xlim([0,.6])
    plt.xlabel('escape rate')
    plt.ylabel('frequency')
    plt.savefig('figures/model_variation_r_logN_' + str(int(np.log10(N))) + '_r_'+str(r)+'_histogram.pdf')
    plt.close()
