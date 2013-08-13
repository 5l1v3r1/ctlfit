#include "cfit.h"
using namespace std;
#include <fstream>
#include <string>

#define DEBUG 1

#define NOTHING 1e-10

double LH(int nlocifit, double *fitness, int nlociseed, double *seedtimes, double N, int nmut, double *mu,
		  double fitprior,double stprior, int ntp, int *tp, int ntp1, int ngt1, int *data, int ntp2, int ngt2,
		  int *samplesizes, int printflag, int trajflag, char *trajfile, char *datfile, char *seedfile, char *LHfile){

	ofstream outfilestream;
	ofstream datfilestream;
	ofstream seedfilestream;
	ofstream LHfilestream;
	if (printflag){
		outfilestream.open(trajfile);
		datfilestream.open(datfile);
		seedfilestream.open(seedfile);
		LHfilestream.open(LHfile);
	}

	double cumfreq[nlocifit+1];		//cumulative frequency of genotypes, need for the seedtime LH
	double freq[nlocifit+1];		//frequencies of genotypes
	double dfreq[nlocifit+1];		//frequency increments
	double perlocus_LH[nlocifit+1];	//array holding the sampling likelihood of different genotypes
	double seed_LH[nlocifit+1];		//array holding the seed time likelihood of different genoytpes
	double data_LH;
	double missing_LH;
	double meanfit;
	double freqsum;		//sum of the frequencies of all tracked genotypes
	double newmeanfit;
	double delta = 1;	//integration time step
	for (int locus=0; locus<nlocifit+1; locus++){cumfreq[locus]=0;freq[locus]=0;}
	freq[0] = 1.0;
	double LH=0;
	double offset = 0;
	double s=0;
	
	double gtfitness[nlocifit+1];
	double gtdummy = 0;
	//we need to compute the genotype fitnesses.
	for (int locus=0; locus<nlocifit+1; locus++){
		if (locus > 0){
			gtdummy += fitness[locus-1];
			LH+=fitprior*fitness[locus-1];
		}
		gtfitness[locus] = gtdummy;
	}
	int tpi = 0;
	for (int gen = 0; gen < tp[ntp-1]+1; gen++){
		meanfit = 0.;
		freqsum = 0.;
		//track the mean fitness. Mean fitness is calculated using only the observed genotypes
		for (int locus=0; locus<nlocifit+1; locus++){
			meanfit+=exp(gtfitness[locus])*freq[locus];
			freqsum+=freq[locus];
		}
		meanfit/=freqsum; //normalize

		/******************************************
		 * evalation of the likelihood of the data
		 * whenever a time point with data is reached
		 ******************************************/
		if (gen == tp[tpi]){ //if we've hit a data point
			int unaccounted_obs = samplesizes[(nlocifit+1)*tpi];
			// of this could be done in the LH_time_point routine, but we want to record the contribution of different loci
			for (int locus=0; locus<nlocifit+1; locus++){	// note that locus really means genotype
				unaccounted_obs-=data[locus+(nlocifit+1)*tpi];
				data_LH = data_LH_perlocus(N, locus, nlocifit, tpi, freq, data, samplesizes); //compute the likelihood contribution from each locus
				perlocus_LH[locus]+=data_LH; 	//save the per-locus likelihood contribution
				if (printflag == 1){			//print stuff to a file for read in from python
					datfilestream << data_LH << '\t';
				}
				LH+=data_LH;					//add it to the total for this data point
			}
			if (unaccounted_obs>0){				//there could be genotypes that are not simulated. those are lumped together as unobserved
				missing_LH =-unaccounted_obs*(log(1.0-freqsum+1e-2)-log(1.0*unaccounted_obs/samplesizes[(nlocifit+1)*tpi]));
			}else{missing_LH=0;}
			if (printflag == 1){
				datfilestream << missing_LH << '\t'; //save the per-locus likelihood contribution
			}
			LH+=missing_LH;
			tpi++;
			if (printflag == 1){
				datfilestream << '\n';
			}
		}


		/******************************************
		 * propagation of the system of differential equations
		 ******************************************/
		dfreq[0] = freq[0]*selection(gtfitness[0], meanfit); //selection on 0th locus
		dfreq[0]-= mu[0]*freq[0]; //mutation away from wild-type

		for (int locus=1; locus<nlocifit+1; locus++){
			if (gen>seedtimes[locus-1]+0.5){ //if we've hit the seed time of this locus
			  //selection on current locus
			  dfreq[locus]=freq[locus]*selection(gtfitness[locus], meanfit);
			  //mutation terms
			  if (locus<nlocifit){ //anything but the last locus gains and looses by mutation
				dfreq[locus]+=mu[locus-1]*(freq[locus-1]-freq[locus]);
			  }
			  else{//the last locus gains only
				dfreq[locus]+=mu[locus-1]*freq[locus-1];
			  }
			}
			else if (abs(gen-seedtimes[locus-1])<0.5 and locus>1){ //if we've just now passed the seed time of a locus
				seed_LH[locus] = stprior*seedtime_LH(N,mu[locus-1], gtfitness[locus], meanfit, freq[locus-1], 
													 cumfreq[locus-1]); //save the likelihood for that locus
				LH+=seed_LH[locus]; //add it to the total likelihood
				if (printflag == 1){
					seedfilestream << seed_LH[locus] << '\n'; //for debugging purposes, output the seed time LH contribution
				}
				dfreq[locus] = (1+selection(gtfitness[locus], meanfit))/N; //no establishment factor. 
				                      //this should be ok since the initial dynamics is mutation fed
			}else{dfreq[locus]=0;} //if we haven't hit this seedtime yet, forget about this locus for now
		}
		for (int locus=0; locus<nlocifit+1; locus++){
			freq[locus]+=delta*dfreq[locus]; //add the accumulated changes to the locus' frequency
			cumfreq[locus]+=delta*freq[locus]; //accumulate the integral of the frequency of each genotype
			//sanity checks on integration to catch instabilities due to numerical inaccuracy
			if (freq[locus] < 0){
				freq[locus] = 0; //nothing should ever be negative
			}
			else if (freq[locus] > 1){
				freq[locus] = 1; //nor greater than one
			}
			if (trajflag == 1){
				outfilestream <<freq[locus]<<'\t'; //save the trajectory, for debugging purposes
			}
		}

		if (printflag == 1){
			outfilestream << '\n';
		}
	}

	// output likelihood details whenever the printflag is set
	if (printflag==1){
		for (int locus=0; locus<nlocifit+1; locus++){
			datfilestream << seed_LH[locus] << '\t'; //save the per-locus likelihood contribution
		}
		datfilestream << 0 << '\n';
	}
	outfilestream.close();
	datfilestream.close();
	seedfilestream.close();
	return LH;
}

/*
 * return growth rate of genotypes
 */
double selection(double gtfitness, double meanfit){
	double increment;
	increment = exp(gtfitness)/meanfit-1.0;
	return increment;
}

/*
 * return the log likelihood of seeding in the current generation (up to a constant)
 * needs as argument the cumulative frequency of the source populations
 * establishment factors in the exponent are ignored, as of now.
 */
double seedtime_LH(double N, double mu, double gtfitness, double meanfit, double freq, double intfreq){
	double ST_LH;
	ST_LH = -log(freq+0*mu)+N*mu*intfreq;
	return ST_LH;
}

/*
 * returns the sampling likelihood of a particular genotype (locus) at a particular time point (tpi)
 */
double data_LH_perlocus(double N, int locus, int nlocifit, int tpi, double *freq, int *data, int *samplesizes){
	double dnu = 1e-2;
	double data_LH = 0.;
	int index;
	double temp_freq;
	index = locus+tpi*(nlocifit+1);

	if (freq[locus] > 1-dnu){	//note the asymmetric treatment of high and low frequencies. We think that alleles rarely go to fixation due to latency etc. The logistic model goes to 1 exponentially resulting in absurdly low LH if founder sequences are sampled late
		temp_freq = 1.0-dnu;
	}
	else{	//to avoid log(0), add a small number. Doesn't matter since it makes negligible contributions to the entropy
		temp_freq = freq[locus]+1e-15;
	}
	data_LH -= data[index]*(log(temp_freq) - log(1.0*(data[index]+1e-8)/samplesizes[index])); // the unaccounted observations are added later in the main routine
	return data_LH;
}

/*
 * function that returns the fixation probability of a single genotype with fitness gtfitness given 
 * a mean fitness of meanfit
 */
double pfix(double gtfitness,double meanfit){
  return 2*selection(gtfitness, meanfit)/(1+2*selection(gtfitness, meanfit));
}
