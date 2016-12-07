# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Configuration

Set these variables:

```
~/.bashrc.ext
export PATH=$PATH:/project/projectdirs/m1187/concorde/TSP/
export WFIRST=/project/projectdirs/m1187/wfirst/
```

* Dependencies

You'll need:

* module load python/2.7-anaconda
* pip install sncosmo --user, sep, anything else...

* Database configuration
* How to run tests
* Deployment instructions



Here’s a typical script for getting a single core on Cori (on the shared queue):

```
#!/bin/bash -l

#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --job-name=unity
#SBATCH —mem=5000

module load python/2.7-anaconda
cd /project/projectdirs/m1187/wfirst/scripts/stan_cosmo/hybrid_6_shared/
srun -n 1 -c 1 python ../STEP1_simulate_survey.py ../paramfiles/paramfile_hybrid.csv pickle_hybrid.txt > log1.txt
```

You can plot the results with:

```
python ../plot_LCs.py pickle_quick.txt
```


maximum_trigger_fraction			SN rates are scaled by this before simulating the survey. Represents high-extinction/near-core that we can’t use/find

adjust_each_SN_exp_time			Adjust exposure time to give same rest-frame S/N, irrespective of extinction?

normalization_wavelength_range		I’ve been running with 5000-6000 rest-frame (~ V band), Charlie has tried 3900-4900 (~ B band)

shallow_SNR						Integrated S/N at maximum light in the above rest-frame range for shallow/medium/deep/reference spectra

medium_SNR

deep_SNR

reference_SNR

spectra_depth_and_phase			depth and rest-frame phase for each spectrum. The code will pick the visit after discovery but closest to this phase for each observation. The reference is always 1 year (observer-frame) after maximum light to reproduce the observation conditions.

targeted_parallels 					Maximize the number of SNe in parallel imaging as a way to choose which SNe to follow

parallel_filters						Which filters should we use for parallel observations?



### Running FoM ###

The following is a typical script for running the FoM. It should run in about 10-14 hours and use 8 cores (for ~200 CPU hours with a change factor of 2).

```
#!/bin/bash -l
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --job-name=unity
#SBATCH --mem=15000
module load python/2.7-anaconda
cd /project/projectdirs/m1187/wfirst/scripts/stan_cosmo/survey_tests_new/survey_00178/FoM_IPC=0.02_nredcoeff=2_fundcal=0.005_crnl=0.0025_include_sys=1_read_noise_floor=4_graydisp=0.08_dark_current=0.01_TTel=282
srun -n 1 -c 8 python ../../../STEP2_SimultUNITY.py -p ../pickle*t -nrestlamb 20 -neigen 1  -gray 0.08 -nredcoeff 2 -IFCIPC 0.02 -IFCdark 0.01 -crnl 0.0025 -fund 0.005 -IFCRNfloor 4 -TTel 282 -mwnorm 0.05 -mwZP 0.005 -mwRV 0.2 -IGext 0.25 > log2.txt
```

* -p tells the code which survey pickle file to use
* -nrestlamb specifies the number of rest-frame wavelengths to model. The nominal number is 100, but 20 (a resolution of ~10) seems to capture most of the cosmological information, and it’s 5x faster.
* -neigen tells the code how many eigenvectors to use. Nominally, we use 5 or 10, but the resolution is too low with -nrestlamb 20. I just set this to 1.
* -gray is the gray dispersion (not counting R_V changes, lensing, or peculiar velocity uncertainty)
* -nredcoeff specifies the number of redshift coefficients for population/R_V drift. 2 allows for an approximately linear drift with scale factor. We’ve also tried 1 (constant with redshift) and 3 (which allows for quite a bit of flexibility).
* -crnl is the size of the count-rate nonlinearity (which is correlated in count-rate). It effectively also sets the size of the calibration difference between nearby and distant SNe.
* -fund is the fundamental color calibration uncertainty (energy vs wavelength). It is correlated in wavelength.
* -mwnorm is the uncertainty on the MW extinction normalization (fractional uncertainty)
* -mwZP is the uncertainty on the MW extinction zeropoint (magnitudes E(B-V))
* -mwRV is the uncertainty on MW R_V
* -IGext is the fractional uncertainty on Brice Menard’s intergalactic extinction model.
* -IFCmaxwave sets the maximum IFC wavelength. That’s useful for colder telescopes.
* -IFCminwave sets the minimum IFC wavelength. That’s 4200A right now, but if we need to make the case for 4200, running survey variants with other values would be useful.


The following are override parameters for the survey (they are all specified when the survey is run). They won’t change the survey (as that has already been run), but they will change the S/N of the existing observations. I’m not sure I recommend setting any of these options to a value different than the survey was originally run with; they’re really there for testing purposes.

* -IFCRNfloor is the read noise floor for the IFC.
* -IFCIPC is the inter-pixel capacitance
* -IFCdark is the dark current
* -TTel is the telescope temperature



### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact