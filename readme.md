### Inferring HIV Escape Rates from Multi-Locus Genotype Data
#### [Front. Immunol., 03 September 2013](http://journal.frontiersin.org/article/10.3389/fimmu.2013.00252/abstract)
#### Taylor A. Kessinger, Alan S. Perelson, Richard A. Neher.

#### Abstract:
Cytotoxic T-lymphocytes (CTLs) recognize viral protein fragments displayed by major histocompatibility complex molecules on the surface of virally infected cells and generate an anti-viral response that can kill the infected cells. Virus variants whose protein fragments are not efficiently presented on infected cells or whose fragments are presented but not recognized by CTLs therefore have a competitive advantage and spread rapidly through the population. We present a method that allows a more robust estimation of these escape rates from serially sampled sequence data. The proposed method accounts for competition between multiple escapes by explicitly modeling the accumulation of escape mutations and the stochastic effects of rare multiple mutants. Applying our method to serially sampled HIV sequence data, we estimate rates of HIV escape that are substantially larger than those previously reported. The method can be extended to complex escapes that require compensatory mutations. We expect our method to be applicable in other contexts such as cancer evolution where time series data is also available.

#### Contents:

##### src/
Contains the class ctl_fit (in ctl_fit.py). This is where most of the fitting is implemented.
The computationally expensive calculation of the likelihood for a given set of estimates is done in cfit.cpp.
Make sure src/ is in sys.path for whichever scripts are run.

	cfit.py
	Implements a forward Euler simulation for our system of ODEs.
	Temporary files for each simulation are output to src/temp_(run number).
	Compile cfit.cpp using make, then SWIGify it.
	
	ctl_fit.py
	The "meat" of the fitting scheme.
	setup must be run every time, as must one of the initial_guesses functions.
	multi_locus_fit is the main script for refining the estimates.
	LH_slice is needed to generate likelihood surfaces.
	MCMC is used for posterior sampling.
	finally, likelihood_c must be run with both optional parameters set to 1 (or true) to generate trajectories.
	The remaining functions are usually called by one of these main wrapping scripts.
	

##### model_fit/
Generates figures for the testing/evaluation of our model, compared against one-parameter and two-parameter logistic fitting.
Most model testing is done with simulation data from the population genetics package FFPopSim:
[get it here](http://neherlab.github.io/ffpopsim/)

	test_data.py
	Contains the function test_data, which generates dummy data from FFPopSim.
	Most of the other scripts call this at some point.

	illustrate_dominant_genotypes.py
	Creates sample plots illustrating the dominance of our selected genotypes over time;
	saved as
	figures/dominant_genotypes_demonstration.png
	figures/dominant_genotypes_demonstration_log.png
	
	illustrate_sequential_fitting.py (make pdfs instead of pngs)
	Shows likelihood surfaces for initial and final estimates, log scale and linear scale of refined fits.
	saved as
	figures/sequential_LH_initial.png
	figures/sequential_traj.png
	figures/sequential_traj_refined.png
	figures/sequential_LH_final_refined.png
	figures/sequential_traj_refined_linear_scale.png (here with larger axis labels)

	sampling_frequency_and_priors_ctl.py
	Use to produce pickle data for varying sample frequencies.
	Remake plots with make_sample_freq_figures.py

	sample_size_and_priors_ctl.py
	Use to produce pickle data for varying sample sizes.
	Remake plots with make_samplesize_figures.py
	
	model_vary_analyze.py
	Use to produce plots based on pickled data from varying model mu and N, as well as simulation r.
	Requires that submit/submit_jobs_model_variation has already been run.
	This obsoletes model_variation_ctl.py.
	
	many_slow.py, varied_CTL_throughout.py, const_CTL_staggered.py, const_CTL_throughout.py
	These perform fits on simulation data under various assumptions, respectively:
	-many loci under weak selection.
	-loci whose selection coefficients vary with time.
	-loci for whom selection turns on at random times.
	-loci under constant selection throughout (i.e., no noise).
	
	test_ctl.py
	Generates a single fit from simple "patient-like" data.
	
##### patient_fit/
Generates figures acquired from the fitting of patient data using our model.

	fit_all_patients_LH.py
	A submit script that calls fit_patients.py with all three patient values.
	Generates trajectories and likelihood plots.

	fit_patients.py
	Contains the function fit_patient_LH, which performs trajectory fitting and LH surface for a single patient.
	
	fit_patients_posterior_only.py
	Produces posterior sampling of escape rates (selection coefficients) for a given patient and parameter set.
	
	plot_all_patient_posteriors.py
	Calls the function plot_patient_posterior (contained in plot_patient_posterior.py) for a range of patients and parameter values.
	Output: posterior distribution plots (based on output from fit_parients_posterior_only.py).

##### submit/
Contains submit scripts appropriate for implementation on a computing cluster.
These generally output pickle files containing aggregated data from many simulation runs.

	simple_wrap.py, submit_script*
	Basic wrappers called by the submit_jobs scripts.
	
	submit_jobs_model_variation.py
	Batch script for generating data based on varying S and F (our prior weights).
	
	submit_jobs_vary_params.py
	Batch script for generating data based on varying N, mu, and r.
	
	fit_all_patients_posterior.py
	Performs posterior sampling for all patients on the cluster.

##### figures/
The output directory for figures.

##### gt_data/
Contains the data used for patient genotypes.
All genotype counts are reconstructed by eyeballing (since the data are half genomes and we don't always know linkage).
Also contains results for each patient given different F priors, tau, and population size.
Each pickle file contains the posterior distribution of the estimated coefficients.
