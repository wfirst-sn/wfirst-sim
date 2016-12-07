import commands
from numpy import *
import pyfits
from scipy.interpolate import interp1d
from DavidsNM import save_img

norm = interp1d([1000., 6000., 10500, 12500, 15500, 20000, 40000],
                [0.88, 0.88, 0.852273685192, 0.837328058564, 0.807873205993, 0.75, 0.75], kind = 'linear')

print "HACK!!!!"*100
norm = lambda x:1.

commands.getoutput("rm -f *.fits")

waves = exp(arange(log(3000.), log(30001.), log(10.)/23.))

for wave in waves:
    dat = zeros([401, 401], dtype=float64)

    dat[200, 200] = 1.0001*norm(wave)/3.
    dat[200, 200 - 30] = norm(wave)/3.
    dat[200, 200 + 30] = norm(wave)/3.
    
    save_img(dat, "%05i_5mas.fits" % wave)

