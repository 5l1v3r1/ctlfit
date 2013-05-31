'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Fits a single patient and produces likelihood plots based on one set of parameters.
          from fit_patients import fit_patient_LH is needed.
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
import argparse

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

def fit_patient_LH(N=1e6, mu=1e-5, tau=20, patient='CH40', F=10):
    
    
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

        #do the actual fitting    
        ctl.multi_locus_fit(ctl.L, 5)
        ctl.likelihood_c(ctl.L, 1, 1) #generate trajectory data
        #plot the trajectories
        plt.figure()
        plt.title(r'$N=10^{'+str(int(np.log10(N)))+r'}$,  $\mu=10^{'+str(int(np.log10(mu)))+r'}$, $\tau='+str(int(tau))+'$')
        model_traj=np.loadtxt('src/temp_test/traj.dat')
        LH=np.loadtxt('src/temp_test/data.dat')
        modelT= np.arange(model_traj.shape[0])-tau
        for locus in xrange(uptolocus+1):
            plt.plot(tp-tau, 1.0*gt[:,locus]/sample_sizes[:,locus], ls='None',marker='o', c=col[locus])
            if locus>0:  plt.plot(modelT, model_traj[:,locus], c=col[locus], label = gts[locus]+r', $\epsilon='+str(round(ctl.fitness[locus-1],2))+r',\;\tau='+str(int(ctl.seedtimes[locus-1]-tau))+'$')
            else:  plt.plot(modelT, model_traj[:,locus], c=col[locus], label = gts[locus])

        ax=plt.gca()
        ax.set_yscale('log')
        plt.ylim([1.0/N, 2])
        plt.xlabel(r'time ($\textrm{days}$)')
        plt.ylabel('genotype frequencies')
        plt.subplots_adjust(bottom=0.15)
        plt.legend(loc=4)
        plt.show()
        plt.savefig('figures/patient_traj/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_'+str(int(np.log10(N)))+'_final_traj.pdf')

        #make LH surfaces for final fit
        plt.figure(figsize = (12,9))
        for locus in xrange(uptolocus):
            ds=ctl.fitness[locus]/20
            dtau=2

            frange = np.linspace(ctl.fitness[locus]-20*ds, ctl.fitness[locus]+20*ds, 41)
            strange = np.linspace(max(0,ctl.seedtimes[locus]-15*dtau), ctl.seedtimes[locus]+15*dtau, 1+ctl.seedtimes[locus]+15*dtau - max(0,ctl.seedtimes[locus]-15*dtau))
            LH,minLH,s,t = ctl.LH_slice(locus, frange, strange, locus+1, False)

            if locus>-1:
                plt.subplot(2,L/2,locus+1)

                findex = np.where(frange==ctl.fitness[locus])
                sindex = np.where(strange==ctl.seedtimes[locus])
                LH[findex, sindex] = minLH + 0.75*dynrange

                plt.title("Epitope "+str(locus+1)+': '+gts[locus+1]) 
                plt.imshow(np.ma.filled(LH,minLH+dynrange), interpolation='nearest',aspect='auto', vmin=minLH, vmax=minLH+dynrange)
                if max(frange)-min(frange) < 0.30:
                    dyticks = round(max(frange)/5, 2)
                elif max(frange)-min(frange) < .4:
                    dyticks = 0.05
                elif max(frange)-min(frange) < .7:
                    dyticks = .1
                elif max(frange)-min(frange) < 1.3:
                    dyticks = .2
                else:
                    dyticks = .3

                ii = (len(frange)-1)*dyticks/(max(frange)-min(frange))
                dyticklocs = np.arange(0,len(frange),ii)

                if min(frange) != 0:
                    dyticklocs += (dyticks-min(frange))%(dyticks)*1.0
                if max(dyticklocs) >= len(frange):
                    dyticklocs = dyticklocs[:-1]

                dytickvals = np.round(np.array(dyticklocs/(len(frange)-1)*(max(frange)-min(frange))+min(frange)), 2)
                if dytickvals[0]<=0:
                    dytickvals[0] = 0.0
                plt.yticks(dyticklocs, map(str, dytickvals))
                print dyticks, dyticklocs, dytickvals

                if max(strange)-min(strange) < 25:
                    dxticks = round((max(strange)-min(strange))/5)
                elif max(strange)-min(strange) < 40:
                    dxticks = 5
                elif max(strange)-min(strange) < 70:
                    dxticks = 10
                elif max(strange)-min(strange) < 100:
                    dxticks = 20
                else:
                    dxticks = 30

                ii = (len(strange)-1)*dxticks/(max(strange)-min(strange))
                dxticklocs = np.arange(0,len(strange),ii)

                if min(strange) > 0:
                    dxticklocs += (dxticks-min(strange))%(dxticks)*1.0
                if max(dxticklocs) >= len(strange):
                    dxticklocs = dxticklocs[:-1]
                dxtickvals = np.around(np.array(1.0*dxticklocs/(len(strange)-1)*(max(strange)-min(strange))+min(strange))).astype(int)
                print dxtickvals, np.array(1.0*dxticklocs/(len(strange)-1)*(max(strange)-min(strange))+min(strange))
                plt.xticks(dxticklocs, map(str, dxtickvals))
                if ((locus)%(L/2) == 0):
                    plt.ylabel(r'escape rate ($\textrm{day}^{-1}$)')
                if locus>=L/2:
                    plt.xlabel(r'seed time ($\textrm{days}$)')

        plt.savefig('figures/patient_LH/'+patient+'_F_'+str(F)+'_S_'+str(S)+'_tau_'+str(tau)+'_logN_'+str(int(np.log10(N)))+'_final_LH.pdf')
        ctl.ctl_clean()
