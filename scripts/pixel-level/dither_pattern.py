from DavidsNM import miniNM_new
from numpy import *
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import sys

def clone(P):
    cloned_x = []
    cloned_y = []

    pad = 1
    for dx in arange(-pad, pad+1):
        for dy in arange(-pad,pad+1):
            cloned_x = concatenate((P[:npts] + dx, cloned_x))
            cloned_y = concatenate((P[npts:] + dy, cloned_y))
    return cloned_x, cloned_y


"""
def chi2fn(P, NA):
    P = P%1.
    cloned_x, cloned_y = clone(P)

    energy = 0.
    for i in range(len(cloned_x)):
        for j in range(i+1, len(cloned_x)):
            energy += 1./(
                (cloned_x[i] - cloned_x[j])**2. + (cloned_y[i] - cloned_y[j])**2.
            )

            #energy += 1./sqrt(
            #    (cloned_x[i] - cloned_x[j])**2. + (cloned_y[i] - cloned_y[j])**2.
            #)
    return energy
"""
def chi2fn(P, NA):
    P = P%1.
    cloned_x, cloned_y = clone(P)

    vor = Voronoi(zip(cloned_x, cloned_y))

    # print vor.vertices
    max_distance_squared = 0.

    for i in range(len(vor.vertices)):
        vx = vor.vertices[i][0]
        vy = vor.vertices[i][1]

        if vx >= 0 and vx < 1 and vy >= 0 and vy < 1:
            this_min_squared = min((vx - cloned_x)**2. + (vy - cloned_y)**2.)
            max_distance_squared = max(max_distance_squared, this_min_squared)

    return max_distance_squared*1000.



npts = int(sys.argv[1])

best_start = None
best_chi2 = 1e10

for i in range(10000):
    ministart = random.random(size = npts*2)
    ministart[0] = 0
    ministart[npts] = 0

    chi2 = chi2fn(ministart, None)
    if chi2 < best_chi2:
        best_chi2 = chi2
        best_start = ministart
        print i, best_chi2

miniscale = [0.1]*(npts*2)
miniscale[0] = 0
miniscale[npts] = 0

bestP, NA, NA = miniNM_new(best_start, miniscale, chi2fn = chi2fn, passdata = None, verbose = False, compute_Cmat = False, maxiter = 500)
cloned_x, cloned_y = clone(bestP)
plt.plot(cloned_x, cloned_y, 'o')
plt.plot([0,1,1,0,0], [0,0,1,1,0], color = 'c')
plt.ylim(-0.5, 1.5)
plt.xlim(-0.5, 1.5)

for i in range(npts):
    print "%.2f, %.2f" % (bestP[i], bestP[npts+i])
plt.axes().set_aspect('equal', 'datalim')

plt.savefig("dithers_%i.pdf" % npts)

