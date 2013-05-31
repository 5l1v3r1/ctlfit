'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Estimates parameters based on varied sample size and prior (F, S) values.
          A prerequisite for make_samplesize_figures.py.
          
'''
import numpy as np
import sys
sys.path.append('src/')
import cPickle as pickle
import ctl_fit
import ctlutils
from ctlutils import params
import os
import errno
from time import strftime
col=ctlutils.col

nrun=10
ssList = [5, 10, 20,50,100]
today = strftime('%y%m%d')

dirname = 'figures/' + today + '_sample_size_logN_'+ str(int(np.log10(params['N']))) + '_logmu_' + str(np.log10(params['mu'])) + '_r_' + str(params['r']) + '_coeffs_'
dirname+="_".join(map(str,params['f']))

try:
    #os.stat(dirname)
    os.makedirs(dirname)
except OSError, e:
    if e.errno != errno.EEXIST:
        raise

command = 'cp *.py '+dirname
os.system(command)

S = float(sys.argv[1])
F = float(sys.argv[2])

estimates_ML = {}
estimates_1p = {}
estimates_2p = {}
for sample_size in ssList:
    estimates_1p[sample_size] = []
    estimates_2p[sample_size] = []
    estimates_ML[sample_size] = []
    for ri in xrange(nrun):
        gt, gt_traj, pop = ctlutils.test_data(params['N'],params['mu'],params['r'],params['f'],params['tp'],params['sample_size'])
        ctl = ctl_fit.ctl_fit()
        ctl.setup(params['N'],params['mu'],gt,params['tp'],params['sample_size']*np.ones_like(gt),F,S, str(np.random.randint(1000000)))
        ctl.initial_guesses_one_parameter_ST()
        estimates_1p[sample_size].append(np.array(ctl.initial_fitness))
        ctl.initial_guesses_two_parameter()
        estimates_2p[sample_size].append(np.array(ctl.initial_fitness))
        ctl.multi_locus_fit(n_iter=5)
        estimates_ML[sample_size].append(np.array(ctl.fitness))
        ctl.ctl_clean()

ML_mean = np.zeros((len(ssList), len(f)))
ML_iqd = np.zeros((len(ssList), len(f)))
p1_mean = np.zeros((len(ssList), len(f)))
p1_iqd = np.zeros((len(ssList), len(f)))
p2_mean = np.zeros((len(ssList), len(f)))
p2_iqd = np.zeros((len(ssList), len(f)))
for si,sample_size in enumerate(ssList):
    estimates_1p[sample_size] = np.array(estimates_1p[sample_size])
    estimates_2p[sample_size] = np.array(estimates_2p[sample_size])
    p1_mean[si,:]  = np.median(estimates_1p[sample_size], axis=0)
    p1_iqd[si,:] = np.percentile(estimates_1p[sample_size], 75, axis=0)-np.percentile(estimates_1p[sample_size], 25, axis=0)
    p2_mean[si,:]  = np.median(estimates_2p[sample_size], axis=0)
    p2_iqd[si,:] = np.percentile(estimates_2p[sample_size], 75, axis=0)-np.percentile(estimates_2p[sample_size], 25, axis=0)

    estimates_ML[sample_size] = np.array(estimates_ML[sample_size])
    ML_mean[si,:]  = np.median(estimates_ML[sample_size], axis=0)
    ML_iqd[si,:] = np.percentile(estimates_ML[sample_size], 75, axis=0)-np.percentile(estimates_ML[sample_size], 25, axis=0)
    

pickle_file = open(dirname+'/sample_size_prior_F_'+str(F)+'_S_'+str(S)+'_means.pickle', 'w')
pickle.dump((ML_mean, ML_iqd, p1_mean, p1_iqd, p2_mean, p2_iqd), pickle_file)
pickle_file.close()
pickle_file = open(dirname+'/sample_size_prior_F_'+str(F)+'_S_'+str(S)+'_full.pickle', 'w')
pickle.dump((estimates_ML, estimates_1p, estimates_2p), pickle_file)
pickle_file.close()
