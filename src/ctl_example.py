import numpy as np
import pylab as plt
import FFPopSim as h


f= np.array([0.5, 0.4, 0.0, 0.0])
L=len(f)
N=1e7
coinf =0.01
mu=1e-5

pop = h.haploid_lowd(L,rng_seed=0)
pop.carrying_capacity=N
pop.set_genotypes([0],[N])
pop.recombination_model = h.CROSSOVERS
pop.outcrossing_rate =  0.5*coinf
pop.set_recombination_rates(0.2)
pop.set_mutation_rates(mu)
pop.set_fitness_additive(f*0.5)

MF = []
for g in xrange(500):
    pop.evolve()
    MF.append(pop.get_allele_frequencies())
    stat=pop.get_fitness_statistics()
    pop.carrying_capacity = N*np.exp(stat.mean)
    if (g==100):
        f= np.array([0.5, 0.4, 0.3, 0.3])
        pop.set_fitness_additive(f*0.5)

    print pop.N
    
    


