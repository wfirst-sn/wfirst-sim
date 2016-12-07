from numpy import *
import cPickle as pickle
import sys

import os
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/pixel-level/")
sys.path.append(wfirst_path + "/scripts/cosmo/")
import astro_functions





def get_fom(picke_file):
    print "Loading pickle file...", pickle_file

    
    (fit_params, stan_data, other_inputs) = pickle.load(open(pickle_file, 'rb'))

    for key in fit_params:
        print "fit_params:", key
    for key in other_inputs:
        print "other_inputs:", key

    if not other_inputs.has_key("jacobian_names"):
        other_inputs["jacobian_names"] = " "*stan_data["nsys"]
        other_inputs["redshifts"] = arange(stan_data["nred"])*0.1 + 0.05
        print "Hack, but fixed in later versions!"

    cmat = cov(transpose(fit_params["mus"]))

    z_set = concatenate(([0.05], arange(len(cmat) - 1, dtype=float64)*0.05 + 0.125))
    print "z_set", z_set
    FoM = astro_functions.get_FoM(cmat, z_set, zp = 0.3, addcmb = 1, verbose = True, shift_constraint = 0.002)[0]
    print "Final_FoM", FoM

    return 


########################################## Starts Here ##########################################

print "Usage: python STEP3_get_fom.py pickle_hybrid.txt"
redo_old_results = 0


pickle_file = sys.argv[1]
get_fom(pickle_file)


print "Done!"

