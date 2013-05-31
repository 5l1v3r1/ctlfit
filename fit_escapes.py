'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Performs MCMC sampling to determine the posterior distribution for the escape rates
          (selection coefficients). Outputs the resultant posteriors as pickle files.
'''

#!/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/bin/python
import numpy as np
import sys
sys.path.append('src/')
import pickle
import ctl_fit
import os
import argparse

#parse the command line arguments
parser = argparse.ArgumentParser(description="Fit sequential escapes to genotype data")
parser.add_argument('--pop', default=1e7, type=float, help='Population size (N)')
parser.add_argument('--mu', default=1e-5,type=float, help='Mutation rate u')
parser.add_argument('--tau', default=0,type=int, help='time offset to correct for deviation from first time point and start of CTL selection, positive tau, earlier CTLS')
parser.add_argument('--input',  default='gt_data/CH40_gts.txt', help='file with table of variant counts')
parser.add_argument('--samplesizes',  default='', help='file with the number of sequences per time point')
parser.add_argument('--output',  default='fit_escape_output', help='directory into which the output is written')
parser.add_argument('--Phi',  default=10,type = int, help='prior against high escape rates')
parser.add_argument('--nsamples', default=1000, type=int, help='number of samples taken from the MC chain')
parser.add_argument('--step', default=1000, type=int, help='number of moves between consecutive samples')
params=parser.parse_args()

# Globals
N = params.pop                        # Population size
mu = params.mu                       # Mutation rate
seq_data_file=params.input
patient_id = seq_data_file.split('/')[-1].split('.')[0]
tau=params.tau
r = 0.00
S=1
Phi=params.Phi
dummy_founder_count = 10             # number of founder sequences initially, makes no difference

#location of the patient data
outdir = params.output
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

    #if no sample sizes are provided, sum the genotype counts, other wise, read from file
    if (params.samplesizes==''):
        sample_sizes = np.repeat([np.sum(gt, axis=1)], gt.shape[1], axis=0).T
    else:
        sample_sizes = dummy_founder_count*np.ones_like(gt)
        sample_sizes[1:,:] = np.loadtxt(params.samplesizes, skiprows=1)[:,1:]

    L=gt.shape[1]-1

    #initialize ctl_fit
    ctl = ctl_fit.ctl_fit()
    ctl.setup(N,mu,gt,tp,sample_sizes,Phi,S)
    ctl.initial_guesses_one_parameter_ST()
    uptolocus=L

    #perform fitting
    ctl.multi_locus_fit(ctl.L, 5)
    ctl.likelihood_c(ctl.L, 1)

    with open(outdir+'/'+patient_id+'_Phi_'+str(Phi)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_'+str(int(np.log10(N)))+'_ML_escape_rates.dat', 'w') as f:
        f.write('\t'.join(['epitope','escape rate','seed time'])+'\n')
        for ii,epi in enumerate(gts[1:]):
            print epi,"escape rate:", ctl.fitness[ii], "seed time", ctl.seedtimes[ii]
            f.write("\t".join(map(str,[epi, ctl.fitness[ii],  ctl.seedtimes[ii]]))+'\n')

    #perform MCMC sampling
    ctl.MCMC(1,params.nsamples*params.step, params.step)
    Fsample = np.array(ctl.MCMC_fitness)
    STsample = np.array(ctl.MCMC_seedtimes)

    #write the posterior data to file
    header_str = "\t".join([escape+' rate' for escape in gts]+[escape+' seed time' for escape in gts])
    try:
        np.savetxt(
            outdir+'/'+patient_id+'_Phi_'+str(Phi)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_'+str(int(np.log10(N)))+'_posterior.dat',
            np.concatenate((Fsample, STsample), axis=1), header=header_str)
    except:
        np.savetxt(
            outdir+'/'+patient_id+'_Phi_'+str(Phi)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_'+str(int(np.log10(N)))+'_posterior.dat',
            np.concatenate((Fsample, STsample), axis=1))

    ctl.ctl_clean()
