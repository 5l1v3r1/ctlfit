'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Produces plots based on varying N, mu, and r.
          Requires that submits/submit_script_vary_params.py has already been run.
          
'''

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import sys
import os
import copy
import glob
from ctlutils import col
from ctlutils import params as default_params
import cPickle as pickle
from time import strftime


params = {'backend': 'ps',  
          'axes.labelsize': 24, 
          'axes.titlesize': 24,
          'text.fontsize': 24,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 17,
'xtick.labelsize': 20,
'ytick.labelsize': 20,
'lines.linewidth': 2,
'text.usetex': True}
plt.rcParams.update(params)



plt.ion()
#initialize global parameters
S=1
F=10.0

N=default_params['N']
#mu=default_params['mu']
simMu=1e-4
mu=simMu
r=default_params['r']
f=default_params['f']
sample_size = default_params['sample_size']
tau=0
#today = '130418'
today = '130727'


#varied sample size; initialize parameters
sizelist = [5, 10,20,50,100]
freqlist = [5, 10, 20, 40, 70, 100]
numruns = 100
#####################################################################
# LOAD DATA
#####################################################################
T=''

fitnesses_size = np.zeros([len(sizelist), numruns, len(f)])
fitnesses_freq = np.zeros([len(freqlist), numruns, len(f)])
fitnesses_size_2p = np.zeros([len(sizelist), numruns, len(f)])
fitnesses_freq_2p = np.zeros([len(freqlist), numruns, len(f)])

for i in range(len(sizelist)):
    dirname = today + '_varsize_'+T+'logN' + str(np.log10(N)) + '_size_' + str(sizelist[i]) +'_simMu_'+str(simMu)+ '_tau_' + str(tau) + '_prior_F_'+str(F)+'_S_'+str(S)
    for k in range(numruns):
        fitnesses_size[i,k,:],fitnesses_size_2p[i,k,:] = pickle.load(open('figures/' + dirname + '/logN_' + str(int(np.log10(N))) + '_logmu_' + str(int(np.log10(mu))) + '_r_'+str(r)+'_runno_' + str(k) +'.pickle', 'r'))[0::2]

for i in range(len(freqlist)):
    dirname = today + '_varfreq_'+T+'logN' + str(np.log10(N)) + '_freq_' + str(freqlist[i]) + '_size_'+str(sample_size) +'_simMu_'+str(simMu)+'_tau_' + str(tau) + '_prior_F_'+str(F)+'_S_'+str(S)
    for k in range(numruns):
        fitnesses_freq[i,k,:],fitnesses_freq_2p[i,k,:] = pickle.load(open('figures/' + dirname + '/logN_' + str(int(np.log10(N))) + '_logmu_' + str(int(np.log10(mu))) + '_r_'+str(0.0)+'_runno_' + str(k) +'.pickle', 'r'))[0::2]

#####################################################################
# PLOTS SHOWING MEDIAN AND STANDARD DEVIATIONS FOR DIFFERENT SAMPLE SIZES AND SAMPLING FREQUENCIES
#####################################################################
lowerp = 10
upperp = 90
for M in ['multi', '2p']:
    fig = plt.figure(figsize=(8,5))
    fig.subplots_adjust(left=0.12, bottom=0.15, top=0.95, right=0.95, wspace = 0.02)
    if M=='multi': 
        fitsize  = fitnesses_size
        fitfreq  = fitnesses_freq
    elif M=='2p':
        fitsize  = fitnesses_size_2p
        fitfreq  = fitnesses_freq_2p
    
    #plt.suptitle('$N = 10^' + str(int(np.log10(N))) + '$, $S = ' + str(S) + '$, $F = ' + str(F) + '$', fontsize=24)
    plt.subplot(1,2,1)
    for l in range(len(f)):
        med = np.mean(fitsize[:,:,l], axis=1)/f[l]
        yerr = [np.std(fitsize[:,:,l],axis=1)/f[l], np.std(fitsize[:,:,l],axis=1)/f[l]]
        plt.errorbar(sizelist, med, yerr = yerr, 
                     label='$\epsilon_' + str(l+1) + ' = ' + str(f[l]) + '$', elinewidth=2, capsize=7, mew=2, c=col[l+1])
    
    ax=plt.gca()
    ax.set_xlim([min(sizelist)-10, max(sizelist)+10])
    ax.set_ylim([0.0,2])
    ax.set_xlabel('sample size')
    ax.set_ylabel('normalized escape rate')
    ax.text(5,1.8,'A', fontsize=32)
    plt.axhline(y=1, color='k', ls='--')
    
    plt.subplot(1,2,2)
    for l in range(len(f)):
        med = np.mean(fitfreq[:,:,l], axis=1)/f[l]
        yerr = [np.std(fitfreq[:,:,l],axis=1)/f[l], np.std(fitfreq[:,:,l],axis=1)/f[l]]
        plt.errorbar(freqlist, med, yerr =yerr, 
                     label='$\epsilon_' + str(l+1) + ' = ' + str(f[l]) + '$', elinewidth=2, capsize=7, mew=2, c=col[l+1])

    ax=plt.gca()
    ax.set_xlim([min(freqlist)-10, max(freqlist)+10])
    ax.set_ylim([0.0,2])
    ax.set_xlabel(r'sampling interval ($\textrm{days}$)')
    ax.set_yticks([])
    ax.text(5,1.8,'B', fontsize=32)
    plt.axhline(y=1, color='k', ls='--')
    
    plt.subplot(1,2,1)
    plt.legend(loc=1,ncol=2, columnspacing=0.2, handletextpad=0.2)
    
    plt.show()
    if M=='multi': 
        plt.savefig('figures/model_sampling_variation_'+T+'ML_logN_' + str(int(np.log10(N))) +'_simMu_'+str(simMu)+ '_F_'+str(int(F))+'_S_'+str(S)+'.pdf')
        plt.savefig('figures/model_sampling_variation_'+T+'ML_logN_' + str(int(np.log10(N))) +'_simMu_'+str(simMu)+ '_F_'+str(int(F))+'_S_'+str(S)+'.svg')
    elif M=='2p':
        plt.savefig('figures/model_sampling_variation_'+T+'2p_logN_' + str(int(np.log10(N))) +'_simMu_'+str(simMu)+ '_F_'+str(int(F))+'.pdf')




#####################################################################
# HISTOGRAMS
#####################################################################


if False:
    bins = np.linspace(0,.6,51)
    for M in ['multi', '2p']:
        fig = plt.figure(figsize=(8,5))
        fig.subplots_adjust(left=0.12, bottom=0.15, top=0.95, right=0.95, wspace = 0.02)
        if M=='multi': 
            fitsize  = fitnesses_size
            fitfreq  = fitnesses_freq
        elif M=='2p':
            fitsize  = fitnesses_size_2p
            fitfreq  = fitnesses_freq_2p

        for si,size in enumerate(sizelist):
            plt.figure()
            plt.title('sample size $' + str(size) + '$, $N = 10^' + str(int(np.log10(N))) + '$, $S = ' + str(S) + '$, $F = ' + str(F) + '$')
            for l in xrange(len(f)):
                y,x=np.histogram(fitsize[si,:,l], bins, normed='True')
                xbins = 0.5*(x[1:]+x[:-1])
                plt.plot(xbins,y, ls='-', c=col[l+1])
                plt.plot([f[l],f[l]], [0,25], lw=2, c=col[l+1])
            ax=plt.gca()
            ax.set_xlim([0,.6])
            plt.xlabel('escape rate')
            plt.ylabel('frequency')
            if M=='multi': 
                plt.savefig('figures/model_variation_size_ML_logN_' + str(int(np.log10(N))) +'_simMu_'+str(simMu)+ '_size_'+str(size)+'_histogram.pdf')
                plt.savefig('figures/model_variation_size_ML_logN_' + str(int(np.log10(N))) +'_simMu_'+str(simMu)+ '_size_'+str(size)+'_histogram.svg')
            elif M=='2p':
                plt.savefig('figures/model_variation_size_2p_logN_' + str(int(np.log10(N))) +'_simMu_'+str(simMu)+ '_size_'+str(size)+'_histogram.pdf')

            plt.close()


        #varied sampling interval  
        size = 20
        bins = np.linspace(0,.6,51)
        for si,sfreq in enumerate(freqlist):
            plt.figure()
            plt.title('sampling interval $' + str(sfreq) + '$, $N = 10^' + str(int(np.log10(N))) + '$, $S = ' + str(S) + '$, $F = ' + str(F) + '$')
            for l in xrange(len(f)):
                y,x=np.histogram(fitfreq[si,:,l], bins, normed='True')
                xbins = 0.5*(x[1:]+x[:-1])
                plt.plot(xbins,y, ls='-', c=col[l+1])
                plt.plot([f[l],f[l]], [0,25], lw=2, c=col[l+1])
            ax=plt.gca()
            ax.set_xlim([0,.6])
            plt.xlabel('escape rate')
            plt.ylabel('frequency')
            if M=='multi': 
                plt.savefig('figures/model_variation_freq_ML_logN_' + str(int(np.log10(N))) + '_freq_' + str(sfreq) +'_simMu_'+str(simMu)+ '_size_'+str(size)+'_histogram.pdf')
                plt.savefig('figures/model_variation_freq_ML_logN_' + str(int(np.log10(N))) + '_freq_' + str(sfreq) +'_simMu_'+str(simMu)+ '_size_'+str(size)+'_histogram.svg')
            elif M=='2p':
                plt.savefig('figures/model_variation_freq_2p_logN_' + str(int(np.log10(N))) + '_freq_' + str(sfreq) +'_simMu_'+str(simMu)+ '_size_'+str(size)+'_histogram.pdf')
                plt.close()
