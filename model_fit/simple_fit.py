'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Generates test data and performs a single fit.
          Used for cluster applications.
          
'''

import numpy as np
import sys
sys.path.append('src/')
import cPickle as pickle
import ctl_fit
import os
import errno
from ctlutils import test_data
from ctlutils import params as default_params
from time import strftime
import FFPopSim as h

import argparse
parser = argparse.ArgumentParser(description="Fit sequential sweeps to genotype data")
parser.add_argument('--simN', default=default_params['N'], type=float, help='simulated population size (N)')
parser.add_argument('--N', default=default_params['N'], type=float, help='Population size (N)')
parser.add_argument('--mu', default=default_params['mu'],type=float, help='Mutation rate u')
parser.add_argument('--simMu', default=default_params['mu'],type=float, help='simulated mutation rate u')
parser.add_argument('--size', default = default_params['sample_size'], type=float, help = 'sample size')
parser.add_argument('--interv', default = 0, type=float, help = 'recombination rate r')
parser.add_argument('--r', default = default_params['r'], type=float, help = 'recombination rate r')
parser.add_argument('--F', default = default_params['F'], type=float, help = 'Fitness prior')

parser.add_argument('--tau', default=0,type=int, help='time offset to correct for deviation from first time point and start of CTL selection; positive tau, earlier CTLs')
parser.add_argument('--runno', default=0, type=int, help='number of run (needed to avoid conflicts)')
parser.add_argument('--runname', default = 'test', type = str)
params=parser.parse_args()

simN = params.simN
N = params.N
r = params.r
mu = params.mu
simMu = params.simMu
tau = params.tau
runno = params.runno
runname = params.runname
sample_size = params.size
sampling_interval = params.interv

F = params.F
S = 1
f = default_params['f']       # Fitness (main/additive effects)
if sampling_interval == 0:
    tp = default_params['tp'] 
    tp = np.array(tp)
else:
    tp = np.arange(0,401,sampling_interval)
tp[1:]+=tau
tp = list(tp)

today = strftime('%y%m%d')
dirname = 'figures/' + today + '_' + runname + '_simMu_'+str(simMu)+'_tau_' + str(tau) + '_prior_F_'+str(F)+'_S_'+str(S)
try:
    #os.stat(dirname)
    os.makedirs(dirname)
except OSError, e:
    if e.errno != errno.EEXIST:
        raise

gt, gt_traj, pop = test_data(simN,simMu,r,f,tp,sample_size)
ctl = ctl_fit.ctl_fit()
ctl.setup(N,mu,gt,tp,sample_size*np.ones_like(gt),F,S, runname + '_logN_' + str(int(np.log10(N))) + '_logmu_' + str(int(np.log10(mu))) + '_r_' + str(r) + '_runno_' + str(runno))

ctl.initial_guesses_two_parameter()
fit_2p = np.array(ctl.initial_fitness)
st_2p = np.array(ctl.initial_seedtimes)

ctl.initial_guesses_one_parameter_ST()
ctl.multi_locus_fit(n_iter=5)


pickle_file = open(dirname+'/logN_' + str(int(np.log10(N))) + '_logmu_' + str(int(np.log10(mu))) + '_r_' + str(r) + '_runno_' + str(runno) +'.pickle', 'w')
pickle.dump((ctl.fitness,ctl.seedtimes, fit_2p, st_2p), pickle_file)
pickle_file.close()

ctl.ctl_clean()
