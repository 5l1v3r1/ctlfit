'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Performs MCMC sampling to determine the posterior distribution for the escape rates
          (selection coefficients). Outputs the resultant posteriors as pickle files.
'''

#!/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/bin/python
import numpy as np
import pylab as plt
import sys
sys.path.append('src/')
import pickle
import ctl_fit
import os
import argparse

#parse the command line arguments
parser = argparse.ArgumentParser(description="Fit sequential sweeps to genotype data")
parser.add_argument('--pop', default=1e8, type=float, help='Population size (N)')
parser.add_argument('--mu', default=1e-5,type=float, help='Mutation rate u')
parser.add_argument('--tau', default=0,type=int, help='time offset to correct for deviation from first time point and start of CTL selection, positive tau, earlier CTLS')
parser.add_argument('--pat',  default='CH40', help='Patient to be fitted')
parser.add_argument('--F',  default=5,type = int, help='prior against high escape rates')
parser.add_argument('--runno', default=0, type=int, help='number of run (needed to avoid conflicts)')
params=parser.parse_args()

# Globals
N = params.pop                        # Population size
mu = params.mu                       # Mutation rate
patient=params.pat
tau=params.tau
r = 0.00
S=1
F=params.F
runno = params.runno
dynrange=20
dummy_founder_count = 10             # number of founder sequences initially, makes no difference   

#location of the patient data
outdir = 'gt_data/'+patient
seq_data_file = 'gt_data/'+patient+'_gts.txt'
try:
    os.mkdir(outdir)
except:
    print "cannot create directory",outdir


try:
    with open(seq_data_file, 'r') as gt_file:
        gts = gt_file.readline().strip().split()[1:]
except IOError, e:
    print e
    gts=[]
    print "File ", seq_data_file, "not found"
except:
    gts=[]
    print "Error opening ", seq_data_file

if len(gts)>0:
    #load data and pull out genotype counts    
    A = np.loadtxt(seq_data_file, skiprows=1)
    gt = np.zeros((A.shape[0]+1, A.shape[1]-1))
    tp = np.zeros((A.shape[0]+1))
    gt[1:,:] = A[:,1:]            #pull out the genotype counts 
    tp[1:] = A[:,0]+tau           #pull out the time points and correct for time of infection
    gt[0,0] = dummy_founder_count #set the number of founder sequences initially, makes no difference

    sample_sizes = dummy_founder_count*np.ones_like(gt)
    sample_sizes[1:,:] = np.loadtxt('gt_data/'+patient+'_sample_sizes.txt', skiprows=1)[:,1:]
    L=gt.shape[1]-1

    #initialize ctl_fit
    ctl = ctl_fit.ctl_fit()
    ctl.setup(N,mu,gt,tp,sample_sizes,F,S, runno)
    ctl.initial_guesses_one_parameter_ST()
    uptolocus=L

    #perform fitting
    ctl.multi_locus_fit(ctl.L, 5)
    ctl.likelihood_c(ctl.L, 1)

    #perform MCMC sampling; 1000000 iterations, 1000 samples
    ctl.MCMC(1,1000000, 1000)
    Fsample = np.array(ctl.MCMC_fitness)
    STsample = np.array(ctl.MCMC_seedtimes)

    #dump the posterior data
    with open(outdir+'/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_'+str(int(np.log10(N)))+'_posterior.pickle', 'w') as pfile:
        pickle.dump((Fsample, STsample, (N,tau, mu)), pfile)

    ctl.ctl_clean()
