from numpy import *
import pystan
import matplotlib.pyplot as plt
import triangle
import commands
import pyfits

def save_img(dat, imname):

    commands.getoutput("rm -f " + imname)
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()


stan_code = """
data {
    int nobs;
    int nsne; // Also equal to the number of galaxies
    matrix [nsne, nobs] flux;
    matrix [nsne, nobs] fluxerr;
    int inds [nsne];
}

parameters {
    real <lower = 0, upper = 2> power; // SN are distributed proportional to trueflux^power
    real <lower = -6, upper = -2> log10_lower_flux;
    real <lower = -1, upper = 1> log10_upper_flux;
    real <lower = 0, upper = 8> flux_power;

    matrix <lower = log10_lower_flux, upper = log10_upper_flux> [nsne, nobs] log10_trueflux;
}


model {
    matrix [nsne, nobs] unnorm_probs;

    // Generate probabilities from truefluxes
    unnorm_probs <- exp(power*log10_trueflux*log(10.)); // unnorm_probs = trueflux^power

    for (i in 1:nsne) {
    // Probability of finding a SN here is equal to the prob
        increment_log_prob(log(unnorm_probs[i, inds[i] + 1]/sum(unnorm_probs[i])   ));
    }

    for (i in 1:nsne) {
        flux[i] ~ normal(exp(log10_trueflux[i]*log(10.)), fluxerr[i]);
    }

    increment_log_prob(flux_power*log(log10_upper_flux - log10_trueflux) - (1. + flux_power)*log(log10_upper_flux - log10_lower_flux) + log(1. + flux_power));


}
"""

nobs = 900 # pixels per galaxy
nsne = 200 # SNe (one per galaxy)
iterations = 400

sn_hosts = glob.glob("/Users/Miles/research/simulated_data/host_with_noise*")
sn_subs = glob.glob("/Users/Miles/research/simulated_data/sub*")

flux = []
inds = []

n = 0
while n < len(sn_hosts):
 
    sn_host = pyfits.open(sn_hosts[n])
    sn_host = sn_host[0].data

    sn_sub = pyfits.open(sn_subs[n])
    sn_sub = sn_sub[0].data
    
    sn_host_1d = sn_host.flatten()
    sn_sub_1d = sn_sub.flatten()
    
    ind = argmax(sn_sub_1d)

    flux.append(sn_host_1d)
    inds.append(ind)
    n += 1



fluxerr = ones((nsne, nobs), dtype=float64)*0.03
obsflux = flux + random.normal(size = [nsne, nobs])*fluxerr

save_img(obsflux, "obsflux.fits")
save_img(flux, "trueflux.fits")


# Run the sampling
fit = pystan.stan(model_code = stan_code, data={"nobs": nobs, "nsne": nsne, "flux": obsflux, "fluxerr": fluxerr, "inds": inds},
                  iter=iterations, chains=4, n_jobs = 4, refresh = 10)
fit_params = fit.extract(permuted = True)

print fit


parameters = ["power", "log10_lower_flux", "log10_upper_flux", "flux_power"]

samples = []
for parameter in parameters:
    samples.append(fit_params[parameter])


figure = triangle.corner(transpose(array(samples)), labels = parameters)

figure.savefig("power.pdf")
