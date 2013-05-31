#include <cmath>
#include <iostream>


//double data_LH_timepoint(double N, int nlocifit, int tpi, double *freq, int *data, int *samplesizes);
double LH(int nlocifit, double *fitness, int nlociseed, double *seedtimes, double N, int nmut,double *mu,  double Fprior,double stprior, int ntp, int *tp, int ntp1, int ngt1, int *data, int ntp2, int ngt2, int *samplesizes, int printflag, int trajflag, char *trajfile, char *datfile, char *seedfile, char *LHfile);
double selection(double gtfitness, double meanfit);
double seedtime_LH(double N, double mu, double gtfitness, double meanfit, double freq, double intfreq);
double data_LH_perlocus(double N, int locus, int nlocifit, int tpi, double *freq, int *data, int *samplesizes);
double pfix(double gtfitness,double meanfit);
