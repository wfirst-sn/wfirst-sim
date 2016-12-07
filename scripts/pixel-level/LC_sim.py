from numpy import *
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pixel_level_ETC import initialize_PSFs, get_imaging_SN, solve_for_exptime
import sys
from FileRead import readcol, writecol
import commands
import argparse
import glob
import sncosmo
from astropy import cosmology
from DavidsNM import save_img


def do_n_point_dither(redshift, exp_time, filt, phase, args):
    total_StoN2 = 0.
    
    for offset_i, offset_j in [(0, 0),
                               (0.5*22, 0.25*22),
                               (0, 0.5*22),
                               (0.5*22, 0.75*22)]:

        total_StoN2 += get_imaging_SN(redshift = redshift, exp_time = exp_time/4., effective_meters2_fl = filt, phase = phase,
                                      offset_i = offset_i, offset_j = offset_j, **args)["PSF_phot_S/N"]**2.
    
    return sqrt(total_StoN2)


def solve_for_exptime_imaging(SN_target, redshift, filt, phase, args):
    exp_times = [1, 10, 100, 1000, 100000]

    SN_vals = [do_n_point_dither(redshift = redshift, exp_time = exp_time, filt = filt, phase = phase, args = args) for exp_time in exp_times]

    last_time = 10
    this_time = 100
    while abs(log(last_time/this_time)) > 0.005:
        last_time = this_time
        this_time = interp1d(SN_vals, exp_times, kind = 'linear')(SN_target)
        exp_times.append(this_time)
        SN_vals.append(do_four_point_dither(redshift = redshift, exp_time = this_time, filt = filt, phase = phase, args = args))

    print exp_times, SN_vals
    return this_time

def realize_SN_model(redshift, x1, c, MV, beta = 3.1, alpha = 0.13):
    sncosmo_model = sncosmo.Model(source="salt2-extended")
    cosmo = cosmology.FlatLambdaCDM(Om0 = 0.3, H0 = 70.)
    ten_pc_z = 2.33494867e-9
    assert abs(cosmo.distmod(z=ten_pc_z).value) < 1.e-3, "Distance modulus zeropoint wrong!"

    ampl = 1.
    for i in range(2):
        sncosmo_model.set(z=ten_pc_z, t0=0., x0=ampl, x1 = x1, c = c)

        mag = sncosmo_model.bandmag('bessellv', 'ab', 0.)
        print "mag, ampl ", mag, ampl
        ampl *= 10.**(0.4*(mag - (MV - alpha*x1 + beta*c)))
        
    mu = cosmo.distmod(redshift).value
    print "mu ", mu

    sncosmo_model.set(z=redshift, t0=0., x0 = ampl*10**(-0.4*mu))
    return sncosmo_model


def make_plot():
    redshift = float(sys.argv[1])

    colors = {"Z087.txt": 'm', "Y106.txt": 'b', "J129.txt": 'g', "H158.txt": 'orange'} # "F184.txt": 'red'

    sncosmo_model = realize_SN_model(redshift = redshift, x1 = 0.0, c = 0.0, MV = -19.08)
    filt_SN_squared = {}
    for key in colors:
        filt_SN_squared[key] = 0.

    for date in arange(-15*(1. + redshift), 45*(1 + redshift), 5):
        phase = date/(1. + redshift)
        f_lamb_SN = sncosmo_model.flux(date, obs_waves)
        args["mdl"] = f_lamb_SN

        for filt in colors:

            #S_to_N = do_four_point_dither(redshift = redshift, exp_time = 150., filt = filt, phase = phase, args = args)
            total_SN = 0
            for dither in range(1):
                ETC_result = get_imaging_SN(redshift = redshift, exp_time = 28.25, effective_meters2_fl = filt, phase = phase,
                                            offset_i = random.random()*22, offset_j = random.random()*22, **args)

                total_SN += ETC_result["PSF_phot_S/N"]**2.
            total_SN = sqrt(total_SN)
            fluxerr = ETC_result["SN_photons/s"]/total_SN
            plt.errorbar(phase + random.random()*0.05, ETC_result["SN_photons/s"] + random.normal()*fluxerr, yerr = fluxerr, fmt = '.', color = colors[filt])
            filt_SN_squared[filt] += total_SN**2.

    plt.title("redshift " + str(redshift) + "  " + "  ".join(["%s=%.1f" % (filt, sqrt(filt_SN_squared[filt])) for filt in colors]))
    plt.show()

"""
PSFs = initialize_PSFs(scales = [22])
obs_waves = arange(6000., 22000., 100.)
args = {"waves": obs_waves, "PSFs": PSFs, "gal_flamb": lambda x: 0}

make_plot()
stop_here
"""

def Kim13_PCs():
    [bands, pcs, phases, mags] = readcol("input/LC_PCs.txt", 'aiff')
    bands = array(bands)

    interpfns = {}
    for i, band in enumerate('griz'):
        plt.subplot(2,2,i+1)
        plt.title(band)
        for pc in range(4):
            inds = where((bands == band)*(pcs == pc))
            phase = phases[inds]
            mag = mags[inds]
            
            phase = concatenate(([-100.], phase, [100.]))
            mag = concatenate(([mag[0]], mag, [mag[-1]]))
            
            interpfns[(band, pc)] = interp1d(phase, mag, kind = 'linear')
            plt.plot(arange(-10., 36.), interpfns[(band, pc)](arange(-10., 36.)), label = str(pc))
        plt.legend(loc = 'best')
    plt.savefig("pc_interps.pdf")
    plt.close()

    redshift = 1.2

    sncosmo_model = realize_SN_model(redshift = redshift, x1 = 0.0, c = 0.0, MV = -19.08 + 2.1*0.0)
    if redshift == 1.2:
        colors = {"Z087": 'm', "Y106": 'b', "J129": 'g', "H158": 'orange'} # "F184": 'red'
        obs_to_pc_filt = {"Z087": 'g', "Y106": 'r', "J129": 'i', "H158": 'z'}
    else:
        colors = {"Y106": 'b', "J129": 'g', "H158": 'orange', "F184": 'red'}
        obs_to_pc_filt = {"Y106": 'g', "J129": 'r', "H158": 'i', "F184": 'z'}

    dates = arange(-10*(1. + redshift), 35*(1 + redshift), 5)
    phases = dates/(1. + redshift)
    jacobian = zeros([len(phases)*4, 6], dtype=float64) # parameters are daymax, mag, pc0123
    jacobian[:,1] = 1.

    weight_matrix = zeros([len(phases)*4]*2, dtype=float64)
    total_SNR_all = 0.

    for i, date in enumerate(dates):
        phase = date/(1. + redshift)
        for j, filt in enumerate(colors):

            total_SN = 0
            f_lamb_SN = sncosmo_model.flux(date, obs_waves)
            args["mdl"] = f_lamb_SN


            for dither in range(1):
                ETC_result = get_imaging_SN(redshift = redshift, exp_time = 100.*(1. + (filt == "F184")), effective_meters2_fl = filt + ".txt", phase = phase,
                                            offset_i = random.random()*22, offset_j = random.random()*22, **args)
                
                total_SN += ETC_result["PSF_phot_S/N"]**2.
            total_SNR_all += total_SN
            total_SN = sqrt(total_SN)

            dmag = (2.5/log(10.))/total_SN
            weight_matrix[i + j*len(dates), i + j*len(dates)] = 1./(dmag**2.)

            f_lamb_SN = sncosmo_model.flux(date + 0.2*(1. + redshift), obs_waves)
            args["mdl"] = f_lamb_SN
            ETC_result_dphase = get_imaging_SN(redshift = redshift, exp_time = 100., effective_meters2_fl = filt + ".txt", phase = phase + 0.2,
                                               offset_i = random.random()*22, offset_j = random.random()*22, **args)

            jacobian[i + j*len(dates), 0] = (ETC_result_dphase["AB_mag"] - ETC_result["AB_mag"])/0.2


            for k in range(4):
                jacobian[i + j*len(dates), k+2] = interpfns[(obs_to_pc_filt[filt], k)](phase)

    total_SNR_all = sqrt(total_SNR_all)
    print weight_matrix
    save_img(jacobian, "jacobian.fits")
    save_img(weight_matrix, "weight_matrix.fits")
    param_wmat = dot(transpose(jacobian), dot(weight_matrix, jacobian))
    param_cmat = linalg.inv(param_wmat)
    print param_cmat
    print total_SNR_all, sqrt(diag(param_cmat))
    
PSFs = initialize_PSFs(scales = [22])
obs_waves = arange(6000., 22000., 100.)
args = {"waves": obs_waves, "PSFs": PSFs, "gal_flamb": lambda x: 0}

Kim13_PCs()
ffklds
"""
"""



def realize_SN(redshift):

    x1 = random.normal()
    c = random.normal()*0.1 + random.exponential()*0.1 - 0.1*log(2.)
    MV = -19.08 + random.normal()*0.12 + random.normal()*0.055*redshift - 0.13*x1 + 2.1*c # beta is 2.1, as this is MV
 
    sncosmo_model = realize_SN_model(redshift = redshift,
                                     x1 = x1,
                                     c = c,
                                     MV = MV)

    filt_SN_squared = {"Y106.txt": 0., "J129.txt": 0., "H158.txt": 0., "F184.txt": 0.}
    ordered_filts = ["Y106.txt", "J129.txt", "H158.txt", "F184.txt"]
    
    for date in arange(-15*(1. + redshift), 45*(1 + redshift), 5):
        phase = date/(1. + redshift)
        f_lamb_SN = sncosmo_model.flux(date, obs_waves)
        args["mdl"] = f_lamb_SN

        for filt in filt_SN_squared:
            #S_to_N = do_four_point_dither(redshift = redshift, exp_time = 150., filt = filt, phase = phase, args = args)

            for dither in range(2):
                            
                """
                import cProfile, pstats, StringIO
                pr = cProfile.Profile()
                pr.enable()
                """

                offset_i = random.random()*22
                offset_j = random.random()*22

                ETC_result = get_imaging_SN(redshift = redshift, exp_time = 99.3, effective_meters2_fl = filt, phase = phase,
                                            offset_i = offset_i, offset_j = offset_j, **args)

                """
                pr.disable()
                s = StringIO.StringIO()
                sortby = 'cumulative'
                ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
                ps.print_stats()
                print s.getvalue()
                fdlkafjdslf
                """

                filt_SN_squared[filt] += ETC_result["PSF_phot_S/N"]**2.
                fStoN.write("%.1f %.1f %s %.4f %.4f\n" % (offset_i, offset_j, filt, ETC_result["AB_mag"], ETC_result["PSF_phot_S/N"]))

    S_to_Ns = array([sqrt(filt_SN_squared[filt]) for filt in ordered_filts])
    print S_to_Ns

    return MV, x1, c, S_to_Ns


PSFs = initialize_PSFs(scales = [22])
obs_waves = arange(6000., 25000., 100.)
args = {"waves": obs_waves, "PSFs": PSFs, "gal_flamb": lambda x: 6e-19}



args = {"waves": obs_waves, "PSFs": PSFs, "gal_flamb": lambda x: 0.0, "mdl": 10**(-0.4*float(sys.argv[1]))*0.10884806248/obs_waves**2., "TTel": 282, "zodi_fl": "aldering.txt"}
dithers = 1
exp_time = float(sys.argv[2])

SNs = []
for i in range(100/dithers):
    summed_SN = 0.
    for j in range(dithers):
        ETC_result = get_imaging_SN(redshift = 0.0, exp_time = exp_time/dithers, effective_meters2_fl = sys.argv[3] + ".txt", phase = 0.0,
                                    offset_i = random.random()*22, offset_j = random.random()*22, verbose = True, **args)

        print ETC_result
        summed_SN += ETC_result["PSF_phot_S/N"]**2.
    SNs.append(sqrt(summed_SN))

print SNs
print "dithers ", dithers, median(SNs), mean(SNs), std(SNs, ddof=1), sys.argv[3], exp_time, sys.argv[1]

fdasfdsafdsfdd
"""




f = open("output/all_SNe.txt", 'a')
fStoN = open("output/S_to_N.txt", 'a')

total_found = 0
while total_found < 100000:
    redshift = random.randint(40)*0.05 + 0.5
    redshift = random.random()**0.5
    redshift = redshift*2.05 + 0.475
    redshift = around(redshift/0.05)*0.05

    MV, x1, c, S_to_Ns = realize_SN(redshift = redshift)

    found = sum(S_to_Ns >= 20.) >= 3.
    f.write("%.3f  %.4f  %.4f  %.4f %.3f %.3f %.3f %.3f\n" % (redshift, MV, x1, c, S_to_Ns[0], S_to_Ns[1], S_to_Ns[2], S_to_Ns[3]))
    f.flush()
    #plt.plot(c, MV, 'o', color = found*'b' + (1 - found)*'r')
    total_found += found
    print "total_found ", total_found
f.close()
fStoN.close()
#make_plot()

"""
