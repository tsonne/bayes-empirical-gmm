# Bayesian Data Analysis Book Chapter

## 4.2. Bayesian Inference of Empirical GMMs

* This documentation and Matlab code by Tim Sonnemann (tsonne@hi.is)
* last update: 2020-12-09

### Abbreviations and Terminology
* PGA = peak ground acceleration
* PSA = pseudo-spectral acceleration, calculated from ground acceleration for single-degree-of-freedom oscillator at various frequencies, maximum value is relevant data
* RJB = Joyner-Boore distance, i.e. horizontal distance from receiver to earthquake source fault
* MW = moment magnitude of earthquake
* SITE = site type, i.e. type of rock under receiver
* DEPTH = hypocentral depth, i.e. depth to earthquake source
* Event = earthquake
* station = acceleration measurement site/device
* GMM = ground motion model, i.e. mathematical model to estimate ground motion caused by earthquake
* MCMC = Markov chain Monte Carlo, sampling method
* DRAM = Delayed-Rejection Adaptive Metropolis sampler

### Structure
* 4.2.1. Icelandic strong-motion data
    * [x] data: (PGA,PSA)(RJB,MW,SITE,DEPTH)(Events,Stations)
		* `Data_PSA_values.csv`
			* PSA values at multiple periods (frequencies) for 83 station-event pairs
			* data were used in Kowsari et al. 2020
			* first row: periods
			* all other rows: each one station-event PSA record
		* `Explanatory.csv`
			* explanatory variables for all station-event pairs
			* columns: station code, eventID, MW, RJB, SITE, DEPTH
    * [x] f: read dataset table
		* `load_42_input.m`
    * [x] f: plot dataset MW vs RJB
		* `plot_mag_dist.m`
* 4.2.2. Bayesian random effects
    * [x] data: GMM functional form, prior probability values
		* original published coefficient values
			* `coeff_Am05.m`
			* `coeff_LL08.m`
		* GMM functions, prior settings in config files
			* `config_Y2.m` for model Y2 (Am05)
			* `config_Y3.m` for model Y3 (LL08)
    * [x] f: read GMM prior data
		* done by `setup_BayesInf.m`
    * [x] f: GMM sim
		* done through config evaluation in `GMM_BayesInf.m`
    * [x] f: global, local optimization (mode search)
		* `simulannealbnd()` from Global Optimization Toolbox
		* `fminsearch` from Optimization Toolbox
		* could also look for public implementations
    * [x] f: proposal covariance matrix
		* based on Hessian at mode from `my_hessian.m`
		* actually using lower/upper parameter value setting to derive suitable proposal standard deviation, Hessian not always stable or reliable
    * [x] f: log-likelihood
		* `lnmvnpdf.m`: multivariate Normal, plain without random effect, i.e. covariance is unity matrix times sigma squared
		* `ln_re_mvnpdf.m`: multivariate Normal, random effect, i.e. intra- and inter-event covariance estimated
    * [x] f: log prior functions
		* `lnPriorNorm.m`: Normal distribution for all
		* `lnPriorUnif.m`: Uniform distribution for all
		* `lnPriorNU.m`: Mixed, some Normal, some Uniform
    * [x] f: MCMC loop, all parallel, must set type in `config_M1.m`
        * [x] f: basic
			* "Roberts1997", fixed covariance, no adaptation
			* can produce garbage if proposal not quite good
			* built in `GMM_BayesInf.m`
        * [x] f: staged, tuned
			* "MSAM", multi-stage adaptive Metropolis
			* primitive tuning by attempting to adapt by acceptance rate, multiple stages with chain combination after each stage adapt more
			* can fail due to numerical issues with covariance matrix
			* `ts_MSAMv2.m`, `wcov.m`, `adapt_tuning.m`
        * [x] f: DRAM
			* "DRAM", delayed rejection adaptive Metropolis
			* not strictly Markovian, but shown to converge to posterior accurately
			* highly self-adapting and efficient
			* `tsdramrun.m`, `covupd.m`
		* [x] f: DRAM + basic
			* "BiDRAM", first use DRAM MCMC for burn-in only, then run basic MCMC with adapted values from DRAM
			* can get same result as in DRAM only, but takes much longer, but is strictly Markovian after burn-in
			* built in `GMM_BayesInf.m`, `tsdramrun.m`, `covupd.m`, `adapt_tuning.m`
* 4.2.3. Results
    * [x] f: combine chains, estimate statistics
		* wrapped in `analyse_MCMC_results.m`
		* Gelman-Rubin statistic in `gpar.m`
		* chain percentiles, mean, median, variance, autocorrelation
		* `random_effect_etai`: estimate event terms
    * [x] f: plot chains
		* plot full chains and autocorrelations up to lag 50
		* plot chain histograms (marginal posteriors)
		* built in `GMM_BayesInf.m`
    * [x] f: plot 2d histograms
		* `kde2d.m`: bivariate (2D) marginal posterior density plots
		* plot parameter correlation coefficient matrix
    * [x] f: plot residuals, bias and event terms
		* built in `GMM_BayesInf.m`
    * [x] f: table of posterior statistics
		* `tsLatexTable.m` for convenient LaTeX table from matrix
    * [x] f: summary plot of residuals and bias through all frequencies
		* `plot_bias_vs_FRM.m`
		* `plot_residuals_vs_MR.m`
	* [x] f: summary tables of posterior parameter statistics
		* built in `GMM_BayesInf.m`, `tsLatexTable`
	* [x] f: write output ASCII table of inferred coefficient values
		* to `dat/(MODEL_NAME)/(MODEL_NAME)_coeff.txt`
	* [x] f: save all chains and statistics to `.mat` files
		* located in `dat/(MODEL_NAME)/`
		* data for each period zipped in archive file

### How to run
* Requires: Matlab, Parallel Computing Toolbox, Global Optimization Toolbox, Optimization Toolbox.
* run `full_ch42` to load data, do inference, analyze, plot all
* `GMM_BayesInf` does almost everything and ties the functions together, sets up LaTeX file, sorts out input configurations, does inference, plots and tables
* output should be generated in local directory
	* new subdirectories are created: `aux dat fig pdf tex`
	* each model will have separate subdirectories in `dat` and files in the other directories
	* LaTeX script will be generated and compilation to PDF file will be attempted (if failed, a warning text should come up)

### Acknowledgments
* Instruction in Bayesian statistics by Professor Birgir Hrafnkelsson
	* Professor of Statistics, University of Iceland
* Borrowed or modified code is listed here:
* Matlab base code and global optimization toolbox
	* required basics, require license
* DRAM functions `tsdramrun.m` and `covupd.m`
	* original author: Marko Laine
	* strongly modified here to allow more advanced form
* `format_ticks.m` to nicely format plot ticks
	* modified from original
	* original from Matlab Central, ref.:
	* Alexander Hayes (2016). Format Tick Labels (https://www.mathworks.com/matlabcentral/fileexchange/15986-format-tick-labels), MATLAB Central File Exchange. Retrieved June 23, 2016.
* `gpar.m` for statistical analysis
	* required to get Gelman-Rubin statistic and estimated effective sample size
	* original code written by Andrew Gelman
	* slight modification to fix minor issue
	* not sure where to find online...
* `ts_lsq.m` for least squares regression with significance test and upper/lower confidence limits of regression ordinate
	* original written by Benedikt Halldorsson
	* slight modification to avoid text output
* `tight_subplot.m` to better control subplot spacing
	* original written by Pekka Kumpulainen
	* Matlab Central ref.:
	* Pekka Kumpulainen (2010). tight_subplot(Nh, Nw, gap, marg_h, marg_w) (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w), MATLAB Central File Exchange. Retrieved June 24, 2010.
* `kde2d` 2D kernel density estimator for 2D probability density plots
	* Matlab Central ref.:
	* Zdravko Botev (2015). kernel density estimation (https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation), MATLAB Central File Exchange. Retrieved December 2, 2015. 
* Kowsari, M., Sonnemann, T., Halldorsson, B., Hrafnkelsson, B., Snæbjörnsson, J. Þ., and Jónsson, S. (2020). Bayesian inference of empirical ground motion models to pseudo-spectral accelerations of south Iceland seismic zone earthquakes based on informative priors. Soil Dynamics and Earthquake Engineering, 132, 106075.

