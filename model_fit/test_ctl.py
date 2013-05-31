import numpy as np
import pylab as plt
import ctl_fit
import os
plt.ion()

col=['k', 'r', 'b', 'g', 'c', 'm']

def test_data(N,mu,r,f,tp, sample_size):
    import FFPopSim as h
    L=len(f)

    # Create population
    pop = h.haploid_lowd(L)

    # Set genotypes, recombination and mutation rates, fitness landscape
    pop.set_genotypes([0],[N])
    pop.outcrossing_rate = r
    pop.set_mutation_rates(mu)
    pop.set_fitness_additive(f)

    # Evolve the population
    genotype_frequencies = []
    while pop.generation<np.max(tp)+1:
        genotype_frequencies.append(pop.get_genotype_frequencies())
        pop.evolve()
    genotype_frequencies = np.asarray(genotype_frequencies)

    gt = np.zeros((len(tp), L+1))
    gt_traj = np.zeros((np.max(tp)+1, L+1))
    for ti,t in enumerate(tp):
        remaining_samples = sample_size
        remaining_freq = 1.0
        for locus in xrange(L+1):
            if genotype_frequencies[t,2**locus-1]/remaining_freq<=1:
                gt[ti,locus] = np.random.binomial(remaining_samples, genotype_frequencies[t,2**locus-1]/remaining_freq)
            else:
                print genotype_frequencies[t,2**locus-1],remaining_freq
                gt[ti,locus] = remaining_samples
            remaining_samples-=gt[ti,locus]
            remaining_freq-=genotype_frequencies[t,2**locus-1]
            if remaining_samples<=0: break
    for locus in xrange(L+1):
        gt_traj[:,locus] = genotype_frequencies[:,2**locus-1]
    print genotype_frequencies.shape
    return gt, gt_traj

def add_initial_guess_to_plot(ctl_instance):
    temp_fit = np.array(ctl_instance.fitness)
    temp_st = np.array(ctl_instance.seedtimes)
    ctl_instance.fitness=ctl_instance.initial_fitness
    ctl_instance.seedtimes=ctl_instance.initial_seedtimes
    ctl_instance.likelihood_c(ctl_instance.L, 1)
    model_traj=np.loadtxt('traj.dat')
    for locus in xrange(ctl.L+1):
        plt.plot(model_traj[:,locus], c=col[locus],ls=':')
    
    ax=plt.gca()
    ax.set_yscale('log')
    plt.ylim([1.0/ctl_instance.N, 2])
    print ctl_instance.initial_fitness, ctl_instance.initial_seedtimes
    ctl_instance.fitness  = np.array(temp_fit)
    ctl_instance.seedtimes  = np.array(temp_st)
    
# Globals
N = 1e7                        # Population size
mu = 1e-5                       # Mutation rate
r = 0.0
dirname = '../N_'+str(N)+'_mu_'+str(mu)+'_r_'+str(r)+'_'
f = 0.5*np.array([0.4, 0.3, 0.15,  0.08, 0.04])       # Fitness (main/additive effects)
dirname+="_".join(map(str,2*f))
try:
    os.mkdir(dirname)
except:
    print "Can't create directory"
S=1
F=1
#tp = np.linspace(0,300,15) 
tp = [0,10,30,60,120,250,400]
sample_size=10
gt, gt_traj= test_data(N,mu,r,f,tp,sample_size)

ctl = ctl_fit.ctl_fit()
ctl.setup(N,mu,gt,tp,sample_size*np.ones_like(gt),F,S)
ctl.multi_locus_fit()

plt.figure()
model_traj=np.loadtxt('traj.dat')
LH=np.loadtxt('data.dat')
for locus in xrange(ctl.L+1):
    plt.plot(gt_traj[:,locus], ls='--', c=col[locus])
    plt.plot(tp, 1.0*gt[:,locus]/sample_size, ls='None',marker='o', c=col[locus])    
    if locus>0:  plt.plot(model_traj[:,locus], c=col[locus], label = r'$s='+str(round(ctl.fitness[locus-1],3))+'$')
    else:  plt.plot(model_traj[:,locus], c=col[locus])

ax=plt.gca()
ax.set_yscale('log')
plt.ylim([1.0/N, 2])
print LH[-1,2:-1]
print LH[:-1,:]
add_initial_guess_to_plot(ctl)
plt.xlabel('time')
plt.legend(loc=4)

if True:
    estimates_ML = {}
    estimates_1p = {}
    estimates_2p = {}
    nrun=20
    SList = [0,1] #,1,5,10,100]
    FList=[1,30,100]
    ssList = [10, 30,100]
    for F in FList:
        for sample_size in ssList:
            estimates_1p[(F,sample_size)] = []
            estimates_2p[(F,sample_size)] = []
            for S  in SList:
                estimates_ML[(F,S,sample_size)] = []
                for ri in xrange(nrun):
                    gt, gt_traj= test_data(N,mu,r,f,tp,sample_size)
                    ctl = ctl_fit.ctl_fit()
                    ctl.setup(N,mu,gt,tp,sample_size*np.ones_like(gt),F,S)
                    ctl.initial_guesses_one_parameter_ST()
                    estimates_1p[F,sample_size].append(np.array(ctl.initial_fitness))
                    ctl.initial_guesses_two_parameter()
                    estimates_2p[F,sample_size].append(np.array(ctl.initial_fitness))
                    ctl.multi_locus_fit()
                    estimates_ML[(F,S,sample_size)].append(np.array(ctl.fitness))
    
    ML_mean = np.zeros((len(FList),len(SList), len(ssList), len(f)))
    ML_std = np.zeros((len(FList),len(SList),len(ssList), len(f)))
    p1_mean = np.zeros((len(FList),len(ssList), len(f)))
    p1_std = np.zeros((len(FList),len(ssList), len(f)))
    p2_mean = np.zeros((len(FList), len(ssList), len(f)))
    p2_std = np.zeros((len(FList), len(ssList), len(f)))
    for fi,F in enumerate(FList):
        for si,sample_size in enumerate(ssList):
            estimates_1p[(F,sample_size)] = np.array(estimates_1p[(F,sample_size)])
            estimates_2p[(F, sample_size)] = np.array(estimates_2p[(F,sample_size)])
            p1_mean[fi,si,:]  = np.median(estimates_1p[(F,sample_size)], axis=0)
            p1_std[fi,si,:] = np.percentile(estimates_1p[(F,sample_size)], 75, axis=0)-np.percentile(estimates_1p[(F,sample_size)], 25, axis=0)
            p2_mean[fi,si,:]  = np.median(estimates_2p[(F,sample_size)], axis=0)
            p2_std[fi,si,:] = np.percentile(estimates_2p[(F,sample_size)], 75, axis=0)-np.percentile(estimates_2p[(F,sample_size)], 25, axis=0)
    for Si,S  in enumerate(SList):
        for fi,F in enumerate(FList):
            for si,sample_size in enumerate(ssList):
                estimates_ML[(F,S,sample_size)] = np.array(estimates_ML[(F,S,sample_size)])
                ML_mean[fi,Si,si,:]  = np.median(estimates_ML[(F,S,sample_size)], axis=0)
                ML_std[fi,Si,si,:] = np.percentile(estimates_ML[(F,S,sample_size)], 75, axis=0)-np.percentile(estimates_ML[(F,S,sample_size)], 25, axis=0)
        
    
    
    for fi,F in enumerate(FList):
        plt.figure()
        plt.title("1 parameter, prior F: "+str(F))
        for locus in xrange(ctl.L):
            plt.errorbar(np.array(ssList)+locus, p1_mean[fi,:,locus]/f[locus]/2, p1_std[fi,:,locus]/f[locus]/2,ls='--', c=col[locus])
        plt.ylim([0,2])
        plt.savefig(dirname+'/prior_'+str(F)+'_single_locus_normed_1p.pdf')
        plt.figure()
        plt.title("2 parameter, prior F: "+str(F))
        for locus in xrange(ctl.L):
            plt.errorbar(np.array(ssList)+locus, p2_mean[fi,:,locus]/f[locus]/2, p2_std[fi,:,locus]/f[locus]/2, ls=':', c=col[locus])
        plt.ylim([0,2])
        plt.savefig(dirname+'/prior_'+str(F)+'_single_locus_normed_2p.pdf')

    for Si,S in enumerate(SList):
        for fi,F in enumerate(FList):
            plt.figure()
            plt.title("prior F: "+str(F)+',S: '+str(S))
            for locus in xrange(ctl.L):
                plt.errorbar(np.array(ssList)+locus, ML_mean[fi,Si,:,locus]/f[locus]/2,ML_std[fi,Si,:,locus]/f[locus]/2, ls='-', c=col[locus])
            plt.ylim([0,2])
            plt.savefig(dirname+'/prior_'+str(F)+'_S_'+str(S)+'_multi_locus_normed.pdf')


if True:
    F = 1
    S = 1
    bins = np.linspace(0,1,51)
    for si,sample_size in enumerate(ssList):
        plt.figure()
        plt.title("prior F: "+str(F)+',S: '+str(S)+', sample size '+str(sample_size))
        for locus in xrange(ctl.L):
            y,x=np.histogram(estimates_ML[(F,S,sample_size)][:,locus], bins, normed='True')
            xbins = 0.5*(x[1:]+x[:-1])
            plt.plot(xbins,y, ls='-', c=col[locus])
            y,x=np.histogram(estimates_1p[(F,sample_size)][:,locus], bins, normed='True')
            plt.plot(xbins,y, ls='--', c=col[locus])
            y,x=np.histogram(estimates_2p[(F,sample_size)][:,locus], bins, normed='True')
            plt.plot(xbins,y, ls=':', c=col[locus])
            plt.plot([2*f[locus],2*f[locus]], [0,5], lw=2, c=col[locus])
    
ctl.ctl_clean()