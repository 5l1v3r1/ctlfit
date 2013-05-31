'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Tests the multilocus model against simulation data where many loci sweep slowly.
'''

import numpy as np
import pylab as plt
import FFPopSim as h
from scipy import optimize as opt

params = {'backend': 'ps',  
          'axes.labelsize': 20, 
          'text.fontsize': 20,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 18,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True}
plt.rcParams.update(params)


def allele_frequencies(f,N, mu,r,tmax):
    L=len(f)
    # Create population
    pop = h.haploid_lowd(L,rng_seed=0)
    pop.carrying_capacity=N
    pop.set_genotypes([0],[N])
    pop.recombination_model = h.CROSSOVERS
    pop.outcrossing_rate =  r
    pop.set_recombination_rates(0.3)
    pop.set_mutation_rates(mu)
    pop.set_fitness_additive(f*0.5)
    
    allele_frequencies = []
    while pop.generation<tmax:
        allele_frequencies.append(pop.get_allele_frequencies())
        pop.evolve()
    
    return np.asarray(allele_frequencies)

def logistic(x,t):
    return np.exp(x[0]*(t-x[1]))/(1.0+np.exp(x[0]*(t-x[1])))

def logistic_deviation(x,tp, data):
    return np.sum((logistic(x,tp)-data)**2)

def fit_logistic(tp,mut_count,sample_size):
    freq = 1.0*mut_count/sample_size
    if (np.sum(freq==1)+np.sum(freq==0) == len(freq)):
        freq[np.where(freq==0)[0][-1]]=1.0/sample_size
        freq[np.where(freq==1)[0][0]]=1-1.0/sample_size
    elif (np.sum(freq==1)+np.sum(freq==0) == len(freq)-1):
        nonzero_index = np.where((freq!=0)*(freq!=1))[0][0]
        if freq[nonzero_index]>0.5:
            freq[nonzero_index-1]=1.0/sample_size
        else:
            freq[nonzero_index+1]=1.0/sample_size

    x=np.zeros(2)
    x[0]=0.05
    try:
        x[1] = tp[np.where(freq>0.5)[0][0]]
    except:
        x[1] = tp[-1]
        
    x=opt.fmin(logistic_deviation, x, args=(tp, freq))
    return x
    

plt.ion()
mu=1e-5
f=np.asarray([0.5, 0.4, 0.4, 0.3, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05])
N=1e7
tmax=501
nrun=50
coinf=0.02
tp = [0,20,30,50,100,275,500]
sample_size = 10
fit_coeff = []

for ri in xrange(nrun):
    fit_coeff.append(np.zeros((len(f),3)))
    fit_coeff[-1][:,0]= f*(1+0.2*np.random.randn(f.shape[0]))
    AF = allele_frequencies( fit_coeff[-1][:,0],N,mu,coinf*0.5,tmax)
    
    for epi in xrange(AF.shape[1]):
        mut_count = np.random.binomial(sample_size, AF[tp,epi])
        fit_coeff[-1][epi,1:] = fit_logistic(tp, mut_count, sample_size)
        #plt.plot(AF[:,epi])
        #plt.plot(tp,1.0*mut_count/sample_size, ls='--', marker='x')

#ax=plt.gca()
#ax.set_yscale('log')
fit_coeff=np.asarray(fit_coeff)

plt.figure()
plt.xlabel(r'true $\epsilon$')
plt.ylabel(r'estimated $\epsilon$')
for ri in xrange(nrun):
    plt.scatter(fit_coeff[ri,:,0],fit_coeff[ri,:,1])

plt.plot(f*1.5,f*1.5, label = 'true coefficients')
plt.xlim([0,0.7])
plt.ylim([0,0.7])
plt.legend()
plt.savefig('true_vs_estimated.pdf')

plt.figure()
plt.xlabel(r'$t_{50}$')
plt.ylabel(r'estimated $\epsilon$')
for ri in xrange(nrun):
    plt.scatter(fit_coeff[ri,:,2], fit_coeff[ri,:,1])

ax=plt.gca()
ax.set_yscale('log')
plt.savefig('t50_vs_estimated.pdf')

