import commands
from numpy import *
import pyfits
from scipy.interpolate import interp1d

def make_file(wave, psf_size, root = ""):
    psfroot = "%05i_5mas%s" % (wave, root)
    lines = orig_lines.replace("OOOOO", psfroot)
    lines = lines.replace("DDDDD", "%.2f" % psf_size)
    lines = lines.replace("MMMMM", "%.4f" % (wave/10000.))
    lines = lines.replace("PPPPP", "%i" % (psf_size*200))
    
    f = open("tmp.par", 'w')
    f.write(lines)
    f.close()

    commands.getoutput("/Applications/tinytim-7.5/tiny2 tmp.par")
    return psfroot + "00_psf.fits"

f = open("tiny.par")
orig_lines = f.read()
f.close()
commands.getoutput("rm -f *.fits")


waves = exp(arange(log(3000.), log(30002.), log(30001./3000.)/4.))
psf_size = 30.
central_600_norms = []

for wave in waves:
    print wave
    flname = make_file(wave, psf_size, root = "_big")

    f = pyfits.open(flname)
    dat = f[0].data
    f.close()

    ind_i, ind_j = where(dat == dat.max())
    ind_i = ind_i[0]
    ind_j = ind_j[0]
    
    print "dat.sum()", dat.sum()
    dat /= dat.sum()

    norm = dat[ind_i - 300: ind_i+300, ind_j - 300: ind_j+300].sum()

    print wave, norm

    central_600_norms.append(norm)

normfn = interp1d(waves, central_600_norms, kind = 'linear')


commands.getoutput("rm -f *_big*fits")



waves = exp(arange(log(3000.), log(30001.), log(10.)/23.))

psf_size = 3.
for wave in waves:
    print wave
    flname = make_file(wave, psf_size)

    f = pyfits.open(flname, 'update')
    f[0].data *= normfn(wave)/f[0].data.sum() 
    print wave, f[0].data.sum()
    f.flush()
    f.close()

    f = pyfits.open(flname)
    assert abs(f[0].data.sum() - normfn(wave)) < 1e-3, "Normalization error!"
    f.close()

commands.getoutput("rm -f *.tt3")
