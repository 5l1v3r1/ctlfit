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
import sys
sys.path.append('./src/')
import ctl_fit
params = {'backend': 'ps',  
          'axes.labelsize': 30, 
          'axes.titlesize': 24,
          'text.fontsize': 24,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 24,
'xtick.labelsize': 24,
'ytick.labelsize': 24,
'lines.linewidth': 2,
'text.usetex': True}
plt.rcParams.update(params)

# list of color codes used for plotting
col=['r', 'b', 'g', 'c','k', 'm','y']

def plot_patient_posterior(Nlist, mu=1e-5, tau=0, patient='CH40', F=10):
    plt.ion()
    S = 1

    ### genotype names are input manually
    ### resu;ts from Vitaly's paper are put in manually to compare to our posteriors
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

    ### load posterior data from pickle files produced by patient_fit/fit_patients_posterior_only.py
    pickle_list = glob.glob(outdir+'/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_mu_'+str(mu)+'_logN_*_posterior.pickle')
        
    # dictionarys that hold the MCMC trajectory for different population sizes
    fitness_posterior = {}
    seedtime_posterior = {}
    LH={}
    #get initial tau, as well as N and mu, from the posterior data
    for filename in pickle_list:
        pfile = open(filename,'r')
        A,B,C,D = pickle.load(pfile)
        pfile.close()
        N,tau1,mu1 = D
        LH[N]=C
        fitness_posterior[N]=A
        seedtime_posterior[N]=B

    #### plot the posterior distribution of escape rates.
    fig=plt.figure(figsize = (12,9))
    L=fitness_posterior[Nlist[0]].shape[1]
    if patient=='CH40':
        L=L-2
    print "number of posterior distributions: ", L
    for locus in xrange(L):
        #take pains to avoid plotting the same thing three times, since CH40 has degenerate loci
        plt.subplot(2,L/2,locus+1)
        plt.title("Epitope "+str(locus+1)+': '+gts[locus])
        if patient == 'CH40' and locus == 1:
            plt.title("Epitopes 2-4: " + gts[locus])
        elif patient == 'CH40' and locus > 1:
            plt.title("Epitope " + str(locus+3)+': '+ gts[locus])
            locus=locus+2
        for ni,N in enumerate(Nlist):
            fitPostLocus=fitness_posterior[N][:,locus]
            #choose sensible bin sizes
            try:
                y,x=np.histogram(fitPostLocus, bins=np.arange(0,min(max(fitPostLocus), 1),0.04*max(fitPostLocus)), normed='True')
                if locus==len(gts)-1 or (patient == 'CH40' and locus == len(gts)+1):
                    plt.plot(0.5*(x[1:]+x[:-1]), y, lw=2,label=r'$N=10^{'+str(int(np.log10(N)))+'}$', c=col[ni])
                else:
                    plt.plot(0.5*(x[1:]+x[:-1]), y, lw=2, c=col[ni])
            except:
                print patient, N, mu,F,tau,L,locus
        ### add axis labels
        if (locus%(L/2) == 0):
            plt.ylabel('escape rate posterior')
        if locus+1>L/2:
            plt.xlabel(r'escape rate ($\textrm{day}^{-1}$)')
        ax=plt.gca()
        ax.set_xlim([0,1])

        ### add vitaly's estimates to the figure
    for locus in xrange(L):
        print gts, locus
        epi=gts[locus].split(',')[0]
        ax=plt.subplot(2,L/2,locus+1)
        plt.plot([vitaly[epi][0],vitaly[epi][0]], [0,3], c='k', lw=3)
        #plt.plot( [vitaly[epi][1],vitaly[epi][2]], [0.2, 0.2], c='k', lw=3)
        
        ax.annotate('', ( vitaly[epi][1], 1),
                    (vitaly[epi][2], 1),
                #xycoords="figure fraction", textcoords="figure fraction",
                #ha="right", va="center",
                arrowprops=dict(arrowstyle='|-|', fc="k", ec="k",lw=3))

    ####### add the legend to one particular subplot of patient CH40
    if patient=='CH40':
        plt.subplot(2,L/2,L)
        plt.legend()
    #fig.suptitle('patient '+patient, fontsize = 28)
    fig.text(.02, .96, 'patient '+patient, fontsize = 36)
    plt.show()
    plt.savefig('figures/patient_posterior/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_mu_'+str(mu)+'_tau_'+str(tau)+'_posterior.pdf')
    #plt.savefig('figures/patient_posterior/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_mu_'+str(mu)+'_tau_'+str(tau)+'_posterior.svg')

    fig = plt.figure(figsize = (12,9))
    for N in Nlist:
        y,x=np.histogram(LH[N], normed='True')
        plt.plot(0.5*(x[1:]+x[:-1]), y, lw=2, label=r'$N=10^{'+str(int(np.log10(N)))+'}$')
    plt.legend()
    plt.xlabel('Likelihood')
    plt.savefig('figures/patient_posterior/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_mu_'+str(mu)+'_tau_'+str(tau)+'_LH_dis.pdf')
    plt.close()

    ####### plot the most likely trajectory
    seq_data_file = 'gt_data/'+patient+'_gts.txt'
    dummy_founder_count = 10 #makes no difference
    for N in Nlist:
        # read the genotype data
        with open(seq_data_file, 'r') as gt_file:
            gts = gt_file.readline().strip().split()[1:]
        A = np.loadtxt(seq_data_file, skiprows=1)
        gt = np.zeros((A.shape[0]+1, A.shape[1]-1))
        tp = np.zeros((A.shape[0]+1))
        gt[1:,:] = A[:,1:]            #pull out the genotype counts 
        tp[1:] = A[:,0]+tau           #pull out the time points and correct for time of infection
        gt[0,0] = dummy_founder_count #set the number of founder sequences initially, makes no difference

        sample_sizes = dummy_founder_count*np.ones_like(gt)
        sample_sizes[1:,:] = np.loadtxt('gt_data/'+patient+'_sample_sizes.txt', skiprows=1)[:,1:]
        L=gt.shape[1]-1
        
        #### set up a ctl fitting instance to integrate the trajectories
        ctl = ctl_fit.ctl_fit()
        ctl.setup(N,mu1,gt,tp,sample_sizes,F,S)

        #### determine the most likely parameters and copy them to ctl fit
        most_likely = np.argmin(LH[N])
        print "most likely parameters at ", most_likely, 'LH:',LH[N][most_likely], fitness_posterior[N][most_likely,:], seedtime_posterior[N][most_likely,:]
        ctl.seedtimes = seedtime_posterior[N][most_likely,:]
        ctl.fitness = fitness_posterior[N][most_likely,:]
        ctl.likelihood_c(ctl.L, 1, 1) #generate trajectory data

        #### plot the trajectories
        fig = plt.figure(figsize = (12,9))
        model_traj=np.loadtxt('src/temp_test/traj.dat')
        traj_LH=np.loadtxt('src/temp_test/data.dat')
        modelT= np.arange(model_traj.shape[0])-tau
        for locus in xrange(L+1):
            plt.plot(tp-tau, 1.0*gt[:,locus]/sample_sizes[:,locus], ls='None',marker='o', markersize=15, c=col[locus])
            if locus>0:  plt.plot(modelT, model_traj[:,locus], c=col[locus], label = gts[locus]+r', $\epsilon='+str(round(ctl.fitness[locus-1],2))+r',\;\tau='+str(int(ctl.seedtimes[locus-1]-tau))+'$')
            else:  plt.plot(modelT, model_traj[:,locus], c=col[locus], label = 'founder')

        ax=plt.gca()
        ax.set_yscale('log')
        plt.ylim([1.0/N, 2])
        plt.xlabel(r'time ($\textrm{days}$)')
        plt.ylabel('genotype frequencies')
        fig.text(.02, .96, 'patient '+patient, fontsize = 36)
        plt.subplots_adjust(bottom=0.15)
        plt.legend(loc=4)
        plt.show()
        plt.savefig('figures/patient_traj/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_mu_'+str(mu)+'_logN_'+str(int(np.log10(N)))+'_most_likely_traj.pdf')
        plt.close()


#if this is the main script, run the above function
if __name__=="__main__":
    plot_patient_posterior([1e6])

