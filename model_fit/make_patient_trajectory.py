'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  produces the trajectories for a given set of escape rates
'''

#!/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/bin/python
import numpy as np
import matplotlib
matplotlib.use('PDF')
import pylab as plt
import sys
sys.path.append('src/')
import cPickle as pickle
import ctl_fit
import os

params = {'backend': 'ps',  
          'axes.labelsize': 20, 
          'axes.titlesize': 20,
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 16,
'xtick.labelsize': 18,
'ytick.labelsize': 18,
'lines.linewidth': 2,
'text.usetex': True}
plt.rcParams.update(params)

def make_patient_trajectories(N=1e6, mu=1e-5, tau=20, patient='CH40', F=10, escape_rates=[], seed_times=[]):
  
    #globals
    S=1
    dynrange=20
    dummy_founder_count = 10             # number of founder sequences initially, makes no difference
    
    plt.ion()
    
    col=['k', 'r', 'b', 'g', 'c', 'm','y']
    #location of the genotype data
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
        ctl.setup(N,mu,gt,tp,sample_sizes,F,S)
        ctl.initial_guesses_one_parameter_ST()
        uptolocus=L

        ctl.multi_locus_fit(ctl.L, 5)
        ctl.likelihood_c(ctl.L, 1, 1) #generate trajectory data
        #plot the trajectories
        plt.figure()
        #plt.title(r'$N=10^{'+str(int(np.log10(N)))+r'}$,  $\mu=10^{'+str(int(np.log10(mu)))+r'}$, $\tau='+str(int(tau))+'$')
        model_traj=np.loadtxt('src/temp_test/traj.dat')
        LH=np.loadtxt('src/temp_test/data.dat')
        modelT= np.arange(model_traj.shape[0])-tau
        for locus in xrange(uptolocus+1):
            plt.plot(tp-tau, 1.0*gt[:,locus]/sample_sizes[:,locus], ls='None',marker='o', markersize=15, c=col[locus])
            if locus>0:  plt.plot(modelT, model_traj[:,locus], c=col[locus], label = gts[locus]+r', $\epsilon='+str(round(ctl.fitness[locus-1],2))+r',\;\tau='+str(int(ctl.seedtimes[locus-1]-tau))+'$')
            else:  plt.plot(modelT, model_traj[:,locus], c=col[locus], label = 'founder')

        ax=plt.gca()
        ax.set_yscale('log')
        plt.ylim([1.0/N, 2])
        plt.xlabel(r'time ($\textrm{days}$)')
        plt.ylabel('genotype frequencies')
        plt.subplots_adjust(bottom=0.15)
        plt.legend(loc=4)
        plt.show()
        plt.savefig('figures/patient_traj/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_'+str(int(np.log10(N)))+'_final_traj.pdf')
