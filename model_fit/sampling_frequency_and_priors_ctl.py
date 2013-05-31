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
sfreqlist = [5, 10, 20,50,100]
today = strftime('%y%m%d')

dirname = 'figures/' + today + '_sampling_frequency_logN_'+ str(int(np.log10(params['N']))) + '_logmu_' + str(np.log10(params['mu'])) + '_r_' + str(params['r']) + '_coeffs_'

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
for sfreq in sfreqlist:
    estimates_1p[sfreq] = []
    estimates_2p[sfreq] = []
    estimates_ML[sfreq] = []
    tp = np.arange(0,401,sfreq)
    for ri in xrange(nrun):
        gt, gt_traj, pop = ctlutils.test_data(params['N'],params['mu'],params['r'],params['f'],params['tp'],params['sample_size'])
        ctl = ctl_fit.ctl_fit()
        ctl.setup(params['N'],params['mu'],gt,params['tp'],params['sample_size']*np.ones_like(gt),F,S, str(np.random.randint(1000000)))
        ctl.initial_guesses_one_parameter_ST()
        estimates_1p[sfreq].append(np.array(ctl.initial_fitness))
        ctl.initial_guesses_two_parameter()
        estimates_2p[sfreq].append(np.array(ctl.initial_fitness))
        ctl.multi_locus_fit(n_iter=5)
        estimates_ML[sfreq].append(np.array(ctl.fitness))
        ctl.ctl_clean()

ML_mean = np.zeros( (len(sfreqlist), len(f)))
ML_iqd = np.zeros((len(sfreqlist), len(f)))
p1_mean = np.zeros((len(sfreqlist), len(f)))
p1_iqd = np.zeros((len(sfreqlist),len(f)))
p2_mean = np.zeros((len(sfreqlist), len(f)))
p2_iqd = np.zeros((len(sfreqlist), len(f)))
for si,sfreq in enumerate(sfreqlist):
    estimates_1p[sfreq] = np.array(estimates_1p[sfreq])
    estimates_2p[sfreq] = np.array(estimates_2p[sfreq])
    p1_mean[si,:]  = np.median(estimates_1p[sfreq], axis=0)
    p1_iqd[si,:] = np.percentile(estimates_1p[sfreq], 75, axis=0)-np.percentile(estimates_1p[sfreq], 25, axis=0)
    p2_mean[si,:]  = np.median(estimates_2p[sfreq], axis=0)
    p2_iqd[si,:] = np.percentile(estimates_2p[sfreq], 75, axis=0)-np.percentile(estimates_2p[sfreq], 25, axis=0)

    estimates_ML[sfreq] = np.array(estimates_ML[sfreq])
    ML_mean[si,:]  = np.median(estimates_ML[sfreq], axis=0)
    ML_iqd[si,:] = np.percentile(estimates_ML[sfreq], 75, axis=0)-np.percentile(estimates_ML[sfreq], 25, axis=0)

pickle_file = open(dirname+'/sampling_interval_prior_F_'+str(F)+'_S_'+str(S)+'_means.pickle', 'w')
pickle.dump((ML_mean, ML_iqd, p1_mean, p1_iqd, p2_mean, p2_iqd,f), pickle_file)
pickle_file.close()
pickle_file = open(dirname+'/sampling_interval_prior_F_'+str(F)+'_S_'+str(S)+'_full.pickle', 'w')
pickle.dump((estimates_ML, estimates_1p, estimates_2p,f), pickle_file)
pickle_file.close()
