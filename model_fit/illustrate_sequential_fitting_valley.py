
import numpy as np
import pylab as plt
import sys
sys.path.append('src/')
import ctl_fit
from ctlutils import test_data
from ctlutils import col
from ctlutils import params as default_params
import os
sys.path.append('/Users/richard/Projects/FFPopSim/FFPopSim/pkg/python/')
import FFPopSim as h
params = {'backend': 'ps',  
          'axes.labelsize': 36, 
          'axes.titlesize': 36,
          'text.fontsize': 36,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 25,
'xtick.labelsize': 30,
'ytick.labelsize': 30,
'lines.linewidth': 3,
'text.usetex': True}
ms=20
plt.rcParams.update(params)

#set up parameters
N = 1e9                      # Population size
mu = default_params['mu']                        # Mutation rate
r = default_params['r'] 
sample_size=default_params['sample_size'] 
tp = [0,20,40,80,150,250]
S=1
F=10
dynrange=20

#########################################
# SET UP FITNESS EFFECTS. f_1>f_2 strong interaction between 3 and 4
#########################################
f = np.array([0.15, 0.05, 0.25,  0.2,  1.0])
#f = 0.5*np.array([0.15, 0.1, 0.15,  0.08,  0.25])
gt = [0b100, 0b1000, 0b1, 0b10, 0b1100]
L=4

# PRODUCE DATA IN EPISTATIC MODE (FITNESS IS PROVIDED AS A LIST OF GENERAL COEFFICIENTS, NOT ONLY ADDITIVE)
gt, gt_traj,pop, agt= test_data(N,mu,r,(gt,f,L),tp,sample_size, all_gts=True, epistatic=True)
#gt, gt_traj,pop, agt= test_data(N,mu,r,f,tp,sample_size, all_gts=True)

# PULL OUT GENOTYPES THAT ARE ABOVE FREQUENCY 0.1 AT LEAST ONCE
gt_freq = [agt[:,x] for x in  range(2**L)]
gt =[]
gts = []
for xi, x in enumerate(gt_freq):
    #print xi, np.max(x)
    if (np.max(x)>0.1):
        gt.append(x)
        gts.append(xi)
gt = np.asarray(gt)

# PRODUCE SAMPLES
gt_samples =np.zeros(( len(tp),len(gt)))
for ti,t in enumerate(tp):
    remaining_samples = sample_size
    remaining_freq = 1.0
    for locus in xrange(len(gt)):
        if gt[locus,t]/remaining_freq<=1:
            gt_samples[ti,locus] = np.random.binomial(remaining_samples, gt[locus,t]/remaining_freq)
        else:
            print genotype_frequencies[t,2**locus-1],remaining_freq
            gt_samples[ti,locus]= remaining_samples
        remaining_samples-=gt_samples[ti,locus]
        remaining_freq-=gt[locus,t]
        if remaining_samples<=0: break


#PLOT SAMPLES AND ACTUAL TRAJECTORIES
fig = plt.figure(figsize = (8,8))
margins = [0.18, 0.14, 0.77, 0.8]
count=0
for traj,gi in zip(gt,gts):
    label_str = L*['-']
    for bi,b in enumerate(str(bin(gi))[2:]):
        if b=='1': label_str[bi]='+'
    plt.plot(tp, gt_samples[:,count]/sample_size, c=col[count], ls='none', marker='o', markersize=ms)
    if gi>0: plt.plot(traj,label=r"$"+"".join(label_str)+"$", c=col[count], lw=3)
    else: plt.plot(traj,label=r"Founder", c=col[count], lw=3)
    count+=1

plt.legend()
ax=plt.gca()
plt.savefig('figures/sequential_valley_traj.pdf')
#ax.set_yscale('log')


######################################
#define mutation rates
#####################################
mutation_rates = [] #[mu]
prev_mut_count=0
for si,sweep in enumerate(gts[1:]):
    mut_count = np.sum(map(int, list(str(bin(sweep))[2:])))
    mutation_rates.append(mu**(mut_count-prev_mut_count))
    prev_mut_count = mut_count

data = np.asarray(gt_samples)
S=1
F=10
dynrange=20
L=data.shape[1]


######################################
# try to fit
######################################
ctl = ctl_fit.ctl_fit()
ctl.setup(N,mutation_rates,data,tp,sample_size*np.ones((len(tp),L)),F,S)
ctl.initial_guesses_one_parameter_ST()
uptolocus=L-1
ctl.multi_locus_fit(5)
ctl.likelihood_c(uptolocus, 1, 1)
fig = plt.figure(figsize = (7,7))
margins = [0.18, 0.14, 0.77, 0.8]
axes = plt.Axes(fig, margins) 
fig.add_axes(axes)
#plt.title(r'$N=10^{'+str(int(np.log10(N)))+'}$, $r='+str(r)+r'$, $\mu=10^{'+str(int(np.log10(mu)))+'}$')
model_traj=np.loadtxt('src/temp_test/traj.dat')
LH=np.loadtxt('src/temp_test/data.dat')
for locus in xrange(uptolocus+1):
    label_str = L*['-']
    for bi,b in enumerate(str(bin(gts[locus]))[2:]):
        if b=='1': label_str[bi]='+'
    plt.plot(gt[locus], ls='--', c=col[locus], lw=3)
    plt.plot(tp, 1.0*data[:,locus]/sample_size, ls='None',marker='o', c=col[locus], markersize=ms)    
    #if locus>0:  plt.plot(model_traj[:,locus], c=col[locus], lw=2,label = r"$"+"".join(label_str)+"$: "+r'$\epsilon_'+str(locus)+'='+str(round(ctl.fitness[locus-1],2))+r',\;\tau_'+str(locus)+'='+str(ctl.seedtimes[locus-1])+'$')
    if locus>0:  plt.plot(model_traj[:,locus], c=col[locus], lw=3,label = r"$"+"".join(label_str)+"$")
    else:  plt.plot(model_traj[:,locus], c=col[locus], lw=3, label='Founder')

plt.ylim([1.0/N, 2])
fig.text(0.03,0.93,'B', fontsize=44)
plt.xlabel(r'time ($\textrm{days}$)')
plt.ylabel('genotype frequencies')
plt.legend(loc=4)
axes.set_xlim([0, tp[-1]*1.05])
axes.set_ylim([0, 1.05])
plt.savefig('figures/sequential_valley_traj_fit.pdf')
ax=plt.gca()
ax.set_yscale('log')
axes.set_ylim([0.5/N, 2.0])
plt.yticks([1e-8, 1e-6, 1e-4, 1e-2, 1])
plt.savefig('figures/sequential_valley_traj_fit_log.pdf')

ctl.ctl_clean()


