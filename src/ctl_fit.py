import numpy as np
from scipy import optimize as op
import cfit
import os
import shutil
import errno
import pdb
import matplotlib.pyplot as plt

class ctl_fit:
    
    def setup(self,N,mu,genotype_counts, time_points, sample_sizes, F=5, S=1, runname='test'):
        '''
        called to set up the instance. inputs
        -N population size
        -mu mutation rate
        -genotype_counts: the number of times a particular genotype is observed (including wild type)
        -time_points: the time points at which the samples were taken
        -sample_sizes: the number of sequences sampled at each of the time points. this is expected to be a 
        matrix with dimensions #timepoints x #genotypes
        -runno: the run number (default to zero; needed only if several fitting runs are active at once)
        -F: the prior used to regularize fitness values
        '''
        self.af_min=1e-30
        self.af_max = 1-1e-2
        self.N = N
        #pdb.set_trace()
        self.L = np.array(genotype_counts).shape[1]-1
        self.ntp = len(time_points)
        if type(mu) in [ type(1), type(1.0)]:
            self.mu = mu*np.ones(self.L)
        else:
            self.mu = np.array(mu)       
        
        self.gt_count = np.array(genotype_counts)
        self.sample_sizes = np.array(sample_sizes)
        self.tp = np.array(time_points,dtype='int')
        self.a_freq = np.ones((len(self.tp), self.L))
        for locus in xrange(self.L):
            #self.a_freq[:,locus] = 1.0*np.sum(self.gt_count[:,(locus+1):], axis=1)/np.sum(self.gt_count[:,:], axis=1)
            self.a_freq[:,locus] = 1.0*np.sum(self.gt_count[:,(locus+1):], axis=1)/self.sample_sizes[:,locus]

        self.identify_degenerate_loci()
        
        self.F=F
        self.S=S
        self.runname = runname
        self.initial_fitness = np.zeros(self.L)
        self.initial_seedtimes = np.zeros(self.L)
        self.fitness = np.zeros(self.L)
        self.seedtimes = np.zeros(self.L)
    
    
    def identify_degenerate_loci(self):
        '''
        loop over loci and spot identical allele frequency trajectories
        these loci are degenerate in the sense that there is no information 
        to fit them separately. 
        '''
        self.sweep_index=np.zeros(self.L, dtype='int')
        for locus in xrange(1,self.L):
            if (self.a_freq[:,locus] == self.a_freq[:,locus-1]).all():
                self.sweep_index[locus]=self.sweep_index[locus-1]
            else:
                self.sweep_index[locus]=self.sweep_index[locus-1]+1
                
    
    def initial_guesses_two_parameter(self):
        '''
        produce initial fitness guesses by optimizing two unconstrained parameters
        in the current implementation those are the fitness coefficient and the seedtime
        '''
        for locus in xrange(self.L):
            xopt = op.fmin_powell(self.single_locus_LH_2p_fit, [0.2, 0], args=(locus,), direc=np.array([[0.01, 0],[0, 1.0]]))
            self.initial_fitness[locus]=xopt[0]
            self.initial_seedtimes[locus]=xopt[1]
        self.fitness = np.array(self.initial_fitness)
        self.seedtimes = np.array(self.initial_seedtimes)
        

    def initial_guesses_one_parameter_ST(self):
        '''
        produce initial fitness guesses by optimizing only the selection 
        coefficient. the seedtime is set to the value that maximizes the prior
        '''
        self.initial_seedtimes[0]=0
        t=np.arange(np.max(self.tp))
        for locus in xrange(self.L):
            if locus>0: xopt = op.fmin(self.single_locus_LH, [0.2], args=(self.initial_seedtimes[locus],locus,self.logistic))
            else: xopt = op.fmin(self.single_locus_LH, [0.2], args=(self.initial_seedtimes[locus],locus,self.mutsel))
            self.initial_fitness[locus]=xopt[0]
            if (locus<self.L-1): 
                freq= (t>self.initial_seedtimes[locus])*self.mutsel(t, xopt[0], self.initial_seedtimes[locus])
                self.initial_seedtimes[locus+1]=np.argmax(freq*np.exp(-self.N*self.mu[locus+1]*np.cumsum(freq)))
                #plt.plot(t, xopt[0]*freq*np.exp(-self.N*self.mu*np.cumsum(freq)))
        self.fitness = np.array(self.initial_fitness)
        self.seedtimes = np.array(self.initial_seedtimes)

    def initial_guesses_one_parameter(self):
        '''
        produce initial fitness guesses by optimizing only the selection 
        coefficient. the seedtime is set to 20*locus
        '''
        self.initial_seedtimes[0]=0
        t=np.arange(np.max(self.tp))
        for locus in xrange(self.L):
            xopt = op.fmin(self.single_locus_LH, [0.2], args=(0,locus,self.logistic))
            self.initial_fitness[locus]=xopt[0]
            if (locus<self.L-1): 
                freq= (t>self.initial_seedtimes[locus])*self.mutsel(t, xopt[0], self.initial_seedtimes[locus])
                self.initial_seedtimes[locus+1]=np.argmax(freq*np.exp(-self.N*self.mu[locus]*np.cumsum(freq)))
        self.fitness = np.array(self.initial_fitness)
        self.seedtimes = np.array(self.initial_seedtimes)



    def single_locus_LH_2p_fit(self,X,locus):
        '''
        wrapper for the single locus likelihood that accepts a length-2 vector 
        as argument that contains the selection coefficient and the seedtime
        '''
        return self.single_locus_LH(X[0], X[1], locus, self.logistic)

    def single_locus_LH(self,s,tau,locus,func):
        '''
        returns the likelihood of the a single locus assuming no interaction with other loci given s and tau
        -s, the selection coefficient.
        -tau time at which the initial condition is imposed
        '''
        a_freq = 1.0*self.a_freq[:,locus]
        af_model = func(self.tp,s,tau)
        af_model[np.where(af_model<self.af_min)]=self.af_min
        af_model[np.where(af_model>self.af_max)]=self.af_max
        LH = np.sum(self.sample_sizes[:,locus] * (a_freq * np.log(af_model) + (1.0-a_freq)*np.log(1.0-af_model)))
        return -LH+self.F*s

    def logistic(self,T,s,tau):
        '''
        computes the single locus trajectories based on the analytical solution to v' = sv(1-v)
        with initial condition v(tau)=1/Ns
        inputs:
        -T array containing the time points at whcih the function is to be evaluated
        -s, the selection coefficient.
        -tau time at which the initial condition is imposed
        '''
        return np.exp(s*(T-tau))/(np.exp(s*(T-tau)) + self.N*s)

    def mutsel(self, T,s,tau, locus=0):
        '''
        computes the single locus trajectories based on the analytical solution to v' = sv(1-v) + m.
        with initial condition v(tau)=0
        inputs:
        -T array containing the time points at whcih the function is to be evaluated
        -s, the selection coefficient.
        -tau time at which the initial condition is imposed
        '''
        root = np.sqrt(s**2+4.0*self.mu[locus]**2)
        alpha = np.arctanh((2.0*self.mu[locus]-s)/root)*2/root    
        freq = 1.0/(2.0*s)*(s-2*self.mu[locus]+np.tanh((alpha+T-tau)/2.0*root)*root)
        return freq
        

    def likelihood_c(self, uptolocus, printout=0, dataout=0):
        '''
        wrapper for cfit.cpp.
        inputs: 
        -uptolocus specifies how many genotypes are to be included into the fit;
        during the first round (where loci are successively added), this will be the current locus.
        it will be L in the later rounds, where single loci are optimized keeping all others fixed
        '''  
        
        #copy all fitness and data into arrays that are passed on the to the c function that integrates
        # and evaluates the likelihood
        fitness = np.array(self.fitness[:uptolocus])
        seedtimes= np.array(self.seedtimes[:uptolocus])

        temp_data = np.zeros([self.ntp, uptolocus+1])
        temp_data[:,:-1] = self.gt_count[:,:uptolocus]
        temp_data[:,-1] = np.sum(self.gt_count[:,uptolocus:], axis=1)
        temp_samplesizes = np.array(self.sample_sizes[:,:uptolocus+1])
        temp_tp = list(self.tp)
        
        try:
            os.makedirs('src/temp_' + str(self.runname) + '/')
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        
        #names of files into which the c-function dumps the trajectories and LH data
        trajpath = 'src/temp_' + str(self.runname) + '/traj.dat'
        datpath = 'src/temp_' + str(self.runname) + '/data.dat'
        seedpath = 'src/temp_' + str(self.runname) + '/seed.dat'
        LHpath = 'src/temp_' + str(self.runname) + '/LH.dat'
        
        #cLH = cfit.LH(fitness, seedtimes, self.N, self.mu,self.F,self.S, temp_tp, temp_data.tolist(), temp_samplesizes.tolist(), int(printout), trajpath, datpath, seedpath, LHpath)
        cLH = cfit.LH(fitness, seedtimes, self.N, self.mu,self.F,self.S, temp_tp, temp_data.tolist(), temp_samplesizes.tolist(), int(printout), int(dataout), trajpath, datpath, seedpath, LHpath)
        
        self.curr_LH = cLH
        return cLH

    def optimize_locus_range(self,lrange, uptolocus, printout=0, dataout=0):
        '''
        optimize the parameters of a range of loci by a greedy down hill search
        -lrange is the a set of loci, all assumed to have identical selection coefficients
        -uptolocus is the total number of loci to be included into the fit
        '''
        self.fitness[lrange]=np.mean(self.fitness[lrange])
        
        #step sizes used to optimize
        ds_list = [0.02, 0.01, 0.001]
        dtau_list = [1]

        #consecutively lower the step sizes 
        for ds in ds_list:
            for dtau in dtau_list:
                old_s=self.fitness[lrange[0]]
                oldLH = self.likelihood_c(uptolocus, 0, dataout)
                if len(lrange)>1:
                    #loop until no better nearest neighbor solution is found (breaks otherwise)
                    while True:
                        #make a list of new LH values, change selection coefficient first
                        new_LH=[]
                        self.fitness[lrange]+=ds
                        new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                        self.fitness[lrange]=old_s
                        if np.min(self.fitness[lrange])>ds:
                            self.fitness[lrange]-=ds
                            new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                        else:
                            new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
    
                        #reassign the old value
                        self.fitness[lrange]=old_s
                        
                        #loop over loci in locus range and attempt to change their seedtimes
                        for locus in lrange:
                            if locus<self.L-1 and locus>0:
                                if self.seedtimes[locus+1] - self.seedtimes[locus]>dtau:
                                    self.seedtimes[locus]+=dtau
                                    new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                                    self.seedtimes[locus]-=dtau   #undo
                                else:
                                    new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
                            elif locus==self.L-1:
                                self.seedtimes[locus]+=dtau
                                new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                                self.seedtimes[locus]-=dtau   #undo
                            else:
                                new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
                                
                            
                            if locus>0:
                                if self.seedtimes[locus] - self.seedtimes[locus-1]>dtau:
                                    self.seedtimes[locus]-=dtau
                                    new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                                    self.seedtimes[locus]+=dtau #undo
                                else:
                                    new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
                            else:
                                new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
    
                        #loop over loci in locus range and attempt to change their seedtimes
                        for locus in lrange:
                            if locus>0:
                                self.seedtimes[locus:]+=dtau
                                new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                                self.seedtimes[locus:]-=dtau   #undo
                            else:
                                new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
                            if locus>0:
                                if self.seedtimes[locus] - self.seedtimes[locus-1]>dtau:
                                    self.seedtimes[locus:]-=dtau
                                    new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                                    self.seedtimes[locus:]+=dtau #undo
                                else:
                                    new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
                            else:
                                new_LH.append(self.likelihood_c(uptolocus, printout, dataout)+1)
    
    
                        #produce a bunch of random moves
                        nrandom=20
                        random_moves = []
                        for ii in xrange(nrandom):
                            oldfitness = np.array(self.fitness)
                            oldseedtimes = np.array(self.seedtimes)
                            self.seedtimes[:uptolocus] += (1-2*(np.random.rand(uptolocus)>0.5))*dtau
                            self.seedtimes[0]=0
                            self.fitness[:uptolocus] +=(1-2*(np.random.rand(uptolocus)>0.5))*ds
                            for swi in xrange(np.max(self.sweep_index)):
                                loci= np.where(self.sweep_index == swi)[0]
                                self.fitness[loci]=np.mean(self.fitness[loci])
                            
                            if (self.seedtimes[:uptolocus] == sorted(self.seedtimes[:uptolocus])).all():
                                new_LH.append(self.likelihood_c(uptolocus, printout, dataout))
                            else:
                                new_LH.append(oldLH+1)
                            random_moves.append([np.array(self.fitness), np.array(self.seedtimes)])
                            self.fitness=oldfitness
                            self.seedtimes=oldseedtimes
                            
                            
                        #determine the change the reduces the LH most
                        best_move = np.argmin(new_LH)
                        minLH=new_LH[best_move]
                        #accept this best move if it increases the LH compared to the previous best LH
                        if (oldLH>minLH):
                            sign = 1-2*(best_move%2)
                            locus=0
                            #print self.seedtimes,
                            if (best_move<2):
                                self.fitness[lrange]+=sign*ds
                            elif (best_move<2*(1+len(lrange))):
                                locus =lrange[0]+ best_move/2-1
                                self.seedtimes[locus]+=sign*dtau
                            elif (best_move<2*(1+2*len(lrange))):
                                locus =lrange[0]+ best_move/2-1-len(lrange)
                                self.seedtimes[locus:]+=sign*dtau
                            else:
                                print "random move"
                                ri = best_move - 2*(1+2*len(lrange))
                                self.fitness=np.array(random_moves[ri][0])                            
                                self.seedtimes=np.array(random_moves[ri][1])
                                
                            old_s=self.fitness[lrange[0]]
                            oldLH = self.likelihood_c(uptolocus, 0, 0)
                            #print oldLH, lrange, ds, dtau, best_move, sign,locus, self.fitness,self.seedtimes
                        else:
                            #break otherwise and the lower step size
                            break
                else:
                    print "LH surface optimization"
                    locus=lrange[0]
                    frange = np.linspace(self.fitness[locus]-5*ds, self.fitness[locus]+5*ds, 11)
                    strange = np.linspace(self.seedtimes[locus]-5*dtau, self.seedtimes[locus]+5*dtau, 11)
                    LH,minLH,s,t = self.LH_slice(lrange[0], frange, strange, uptolocus, False)
                    if (minLH<oldLH) and not np.isnan(minLH):
                        print 'LH surface: ', minLH,s,t
                        if (locus==0):
                            si = np.where(strange==0)[0]
                            self.fitness[locus]=frange[np.argmin(LH[:,si])]
                            self.seedtimes[locus]=0
                            oldLH = np.min(LH[:,si])
                        else:
                            self.fitness[locus]=s
                            self.seedtimes[locus]=t
                            oldLH=minLH
                    else:
                        break

    def MCMC(self, T, n, dn):
        '''
        optimize the parameters of a range of loci by a greedy down hill search
        -lrange is the a set of loci, all assumed to have identical selection coefficients
        -uptolocus is the total number of loci to be included into the fit
        '''
        
        #step sizes used to optimize
        ds=0.01
        dtau=1
        self.MCMC_fitness = []
        self.MCMC_seedtimes = []
         
        #consecutively lower the step sizes 
        old_LH = self.likelihood_c(self.L, 0, 0)
        uptolocus=self.L
        cycle=0
        while cycle<n:
            cycle+=1
            oldfitness = np.array(self.fitness)
            oldseedtimes = np.array(self.seedtimes)
            self.seedtimes[:uptolocus] += (1-2*(np.random.rand(uptolocus)>0.5))*dtau
            self.seedtimes[0]=0
            self.fitness[:uptolocus] +=(1-2*(np.random.rand(uptolocus)>0.5))*ds
            for sweep in range(np.max(self.sweep_index)+1):
                lrange= np.where(self.sweep_index==sweep)[0]
                self.fitness[lrange]=np.mean(self.fitness[lrange])
            
            if (self.seedtimes[:uptolocus] == sorted(self.seedtimes[:uptolocus])).all():
                new_LH = self.likelihood_c(uptolocus, 0, 0)
            else:
                new_LH = self.likelihood_c(uptolocus, 0, 0)+1e5
            if np.random.rand()>np.exp((old_LH-new_LH)/T):
                self.fitness=oldfitness
                self.seedtimes=oldseedtimes
            else:
                old_LH=new_LH
            
            if cycle%dn==0:
                self.MCMC_fitness.append(np.array(self.fitness))
                self.MCMC_seedtimes.append(np.array(self.seedtimes))
                print cycle, self.fitness, self.seedtimes, old_LH
        
    def multi_locus_fit(self,n_iter=5, printout=0, dataout=0):
        '''
        fit by optimizing the multi-locus LH
        -n_iter are the number of times each locus is optimized
        '''

        self.initial_guesses_one_parameter_ST()
        if np.sum(self.fitness==sorted(self.fitness,reverse=True))!=self.L:
            print "Funny fitnesses:", self.initial_fitness
            print "Falling back onto simple estimates: ",
            self.initial_guesses_one_parameter()
            print self.initial_fitness, self.initial_seedtimes
        
        #do sequential fitting first and successively increase the number of loci to be included
        for locus in xrange(np.max(self.sweep_index)+1):
            #determine the set of loci that have the same sweep pattern, those are fitted together with the same coefficient
            locus_range = np.where(self.sweep_index==locus)[0]
            if len(locus_range)==1:
                locus=locus_range[0]
                if locus>0:
                    ds=self.initial_fitness[locus]/10
                    dtau=max(int(self.initial_seedtimes[locus]/15),1)
                    frange = np.linspace(self.initial_fitness[locus]-10*ds, self.initial_fitness[locus]+10*ds, 21)
                    strange = np.arange(max(0,self.initial_seedtimes[locus]-15*dtau), self.initial_seedtimes[locus]+15*dtau, dtau)
                    LH,minLH,s,t = self.LH_slice(locus, frange, strange, locus+1, False)
                    self.fitness[locus]=s
                    self.seedtimes[locus]=t
#                    if (locus==0):
#                        si = np.where(strange==0)[0]
#                        self.fitness[locus]=frange[np.argmin(LH[:,si])]
#                        self.seedtimes[locus]=0
            else:
                self.optimize_locus_range(locus_range, np.max(locus_range)+1)
            print "preiteration locus:", locus_range, self.likelihood_c(locus+1, 1, 1)
        
        #loop over all loci again n_iter times to refine the estimate.
        for ii in xrange(n_iter):
            for locus in xrange(np.max(self.sweep_index)+1):
                locus_range = np.where(self.sweep_index==locus)[0]
                self.optimize_locus_range(locus_range, self.L)
                print "iteration",ii,"locus:", self.likelihood_c(self.L, 1, 1)
        
        print self.fitness, self.seedtimes

    def LH_slice(self, locus, fit_range, seed_range, uptolocus, makefigure=False):
        oldfitness = np.array(self.fitness)
        oldseedtimes = np.array(self.seedtimes)
        
        LH = np.zeros((len(fit_range), len(seed_range)))
        for fi, f in enumerate(fit_range):
            for si, st in enumerate(seed_range):
                self.fitness[locus]=f
                self.seedtimes[locus]=st
                if (self.seedtimes[:uptolocus] == sorted(self.seedtimes[:uptolocus])).all():
                    LH[fi,si] = self.likelihood_c(uptolocus, 0, 0)
                else:
                    LH[fi,si]=np.nan

        LH=np.ma.masked_invalid(LH)
        self.fitness=oldfitness
        self.seedtimes=oldseedtimes
        minLH = np.min(LH.flatten())
        try:
            fi,si = np.where(LH==minLH)
        except:
            print "messed up LH surface", LH
            return LH, np.nan, 0,0
        print "Minimal neg  LH", minLH,"at s=",np.round(fit_range[fi[0]],3), "and t=",seed_range[si[0]]
        if makefigure:
            plt.figure()
            plt.title("Minimal neg  LH "+str(np.round(minLH,3))+" at s="+str(np.round(fit_range[fi[0]],3))+" and t="+str(seed_range[si[0]]))
            plt.imshow(LH, interpolation='nearest')
            ii=int(len(fit_range)/10)+1
            plt.yticks(np.arange(0,len(fit_range),ii), map(str, np.round(fit_range[::ii],3)))
            ii=int(len(seed_range)/10)+1
            plt.xticks(np.arange(0,len(seed_range),ii), map(str, seed_range[::ii]))
            plt.colorbar()
        return LH, minLH, fit_range[fi[0]], seed_range[si[0]]

    def ctl_clean(self):
        if os.path.isdir('src/temp_' + str(self.runname) + '/'):
            shutil.rmtree('src/temp_' + str(self.runname) + '/')
