'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Function for generating test data used for fitting.
          from test_data import test_data is needed.
          
'''

import numpy as np
import sys
import FFPopSim as h

col=['k', 'r', 'b', 'g', 'c', 'm']

params = {'f': np.array([0.5, 0.4, 0.25,  0.15 , 0.08]),
          'N': 1e7,
          'mu': 1e-5,
          'F': 10.0,
          'S': 1.0,
          'r':0.0,
          'tp':[0,20,40,60,120,250],
          'sample_size':20}


def test_data(N,mu,r,f,tp, sample_size, all_gts=False, epistatic=False):
    if (epistatic):
        L=f[2]
    else:
        L=len(f)
        
    # Create population
    pop = h.haploid_lowd(L,rng_seed=0)

    # Set genotypes, recombination and mutation rates, fitness landscape
    pop.set_genotypes([0],[N])
    pop.outcrossing_rate = r
    pop.set_mutation_rates(mu)
    if (epistatic):
        pop.set_fitness_coefficients(f[0],0.5*f[1])        
    else:
        pop.set_fitness_additive(0.5*f)

    # Evolve the population
    genotype_frequencies = []
    while pop.generation<np.max(tp)+1:
        genotype_frequencies.append(pop.get_genotype_frequencies())
        pop.evolve()
    genotype_frequencies = np.asarray(genotype_frequencies)
    print genotype_frequencies
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
    
    #option: either return all genotypes or just the dominant ones.
    if all_gts:
        return gt, gt_traj, pop, genotype_frequencies
    else:
        return gt, gt_traj, pop 




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
    
