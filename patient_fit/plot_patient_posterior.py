'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  plots patient posteriors for range of F, tau, and N values.
'''

import numpy as np
import pylab as plt
import glob
import pickle
import argparse
params = {'backend': 'ps',  
          'axes.labelsize': 24, 
          'axes.titlesize': 24,
          'text.fontsize': 24,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 20,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'lines.linewidth': 2,
'text.usetex': True}
plt.rcParams.update(params)

def plot_patient_posterior(Nlist, mu=1e-5, tau=0, patient='CH40', F=10):

    plt.ion()
    S = 1

    col=['k', 'r', 'b', 'g', 'c', 'm','y']
    
    #genotype names are input manually
    if patient == 'CH40':
        gts = ['nef185', 'gag113, gag389, vpr74', 'vif161', 'env145']
        vitaly = {'nef185':[0.22,0.112,0.536], 'gag113': [0.17, 0.011, 0.468],
                  'gag389':[0.17, 0.012, 0.443], 'vpr74':[0.16 ,0.14,0.516], 
                  'vif161':[0.37 ,0.115,0.379], 'env145':[0,0,0]}
    elif patient == 'CH58':
        gts = ['env581', 'env830', 'nef105', 'gag236']
        vitaly = {'env581':[0.1,0.005,0.808], 'env830':[0.12 ,0.074,0.424], 'nef105':[0.07,0.006,0.132], 'gag236': [0.08, 0.014, 0.15]}
    elif patient == 'CH77':
        gts = ['tat55', 'env350', 'nef17', 'nef73']
        vitaly = {'tat55':[0.42, 0.367, 1.177], 'env350':[0.36 ,0.244,0.743], 'nef17':[0.3 ,0.139,0.458], 'nef73':[0.29,0.236,0.834]}
    outdir = 'gt_data/' + patient
    
    #load posterior data
    pickle_list = glob.glob(outdir+'/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_*_posterior.pickle')
    fitness_posterior = {}
    seedtime_posterior = {}
    
    #get initial tau, as well as N and mu, from the posterior data
    for filename in pickle_list:
        pfile = open(filename,'r')
        A,B,C = pickle.load(pfile)
        pfile.close()
        N,tau1,mu1 = C
        fitness_posterior[N]=A
        seedtime_posterior[N]=B
    
    
    fig=plt.figure(figsize = (12,9))
    for N in Nlist:
        L=fitness_posterior[N].shape[1]
        print "number of posterior distributions: ", L
        if patient=='CH40':
            L=L-2
        for locus in xrange(L):
            #take pains to avoid plotting the same thing three times, since CH40 has degenerate loci
            plt.subplot(2,L/2,locus+1)
            plt.title("Epitope "+str(locus+1)+': '+gts[locus])
            if patient == 'CH40' and locus == 1:
                plt.title("Epitopes 2-4: " + gts[locus])
            elif patient == 'CH40' and locus > 1:
                plt.title("Epitope " + str(locus+3)+': '+ gts[locus])
                locus=locus+2
            D=fitness_posterior[N][:,locus]
            #choose sensible bin sizes
            y,x=np.histogram(D, bins=np.arange(0,min(max(D), 1),0.04*max(D)), normed='True')
            if locus==len(gts)-3 or (patient == 'CH40' and locus == len(gts)-1):
                plt.plot(0.5*(x[1:]+x[:-1]), y, lw=2,label=r'$N=10^{'+str(int(np.log10(N)))+'}$')
            else:
                plt.plot(0.5*(x[1:]+x[:-1]), y, lw=2)
            if (locus%(L/2) == 0):
                plt.ylabel('escape rate posterior')
            if locus+1>L/2:
                plt.xlabel(r'escape rate ($\textrm{day}^{-1}$)')
            ax=plt.gca()
            ax.set_xlim([0,1])

    for locus in xrange(L):
        epi=gts[locus].split(',')[0]
        ax=plt.subplot(2,L/2,locus+1)
        plt.plot([vitaly[epi][0],vitaly[epi][0]], [0,3], c='k', lw=3)
        #plt.plot( [vitaly[epi][1],vitaly[epi][2]], [0.2, 0.2], c='k', lw=3)
        
        ax.annotate('', ( vitaly[epi][1], 1),
                    (vitaly[epi][2], 1),
                #xycoords="figure fraction", textcoords="figure fraction",
                #ha="right", va="center",
                arrowprops=dict(arrowstyle='|-|',
                                fc="k", ec="k",lw=3
                                ),
                    )

    plt.subplot(2,L/2,2)
    plt.legend()
    #fig.suptitle('patient '+patient, fontsize = 28)
    fig.text(.02, .96, 'patient '+patient, fontsize = 36)
    plt.show()
    plt.savefig('figures/patient_posterior/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_posterior.pdf')
    plt.savefig('figures/patient_posterior/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_posterior.svg')

#if this is the main script, run the above function
if __name__=="__main__":
    plot_patient_posterior([1e6])

