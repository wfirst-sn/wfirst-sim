from numpy import *
import pyfits
import os
import sys
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/cosmo/")
from astro_functions import get_FoM

for item in sys.argv[1:]:
    f = pyfits.open(item)
    cmat = f[0].data
    f.close()

    cmat -= cmat.min()

    """
    for i in range(len(cmat)):
        for j in range(len(cmat)):
            if i != j:
                cmat[i,j] = 0
    """

    z_list = concatenate(([0.05], arange(len(cmat) - 1)*0.05 + 0.125))
    print z_list

    FoM = get_FoM(cmat, z_list)[0]

    print "FoM", item, FoM
    print "Diagonal Only ", get_FoM(diag(diag(cmat)), z_list)
    print "Diagonal Only (Median Removed) ", get_FoM(diag(diag(cmat - median(cmat))), z_list)
    print "Half the number of SNe ", get_FoM(cmat*2., z_list)



