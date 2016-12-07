import commands
from numpy import *
import pyfits
from scipy.interpolate import RectBivariateSpline
from scipy.special import jv
from DavidsNM import save_img

print "This is a dumb hack to get PSFS while we wait for better ones."
f = open("tiny.par")
orig_lines = f.read()
f.close()
commands.getoutput("rm -f *.fits")

waves = exp(arange(log(3000.), log(30001.), log(10.)/23.))

for wave in waves:
    xs, ys = meshgrid(arange(-200., 201.)*0.005, arange(-200., 201.)*0.005)


    rs = sqrt(xs**2. + ys**2.)
    rs = clip(rs, 1e-4, 1e10)

    psf = (
        3*jv(1, rs*1e5*3.14159**2./(9.*wave)) - 10.*jv(1, 1.e6 * 3.14159**2. * rs/(27.*wave))
        )**2./(91.*3.14159*rs**2.)

    psf *= 0.005**2.

    print sum(psf)

    save_img(psf, "%05i_5mas.fits" % wave)
    
    

