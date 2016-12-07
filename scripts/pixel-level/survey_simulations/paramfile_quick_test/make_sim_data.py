from numpy import *
import cPickle as pickle


z_array = []
for z in arange(0.05, 0.96, 0.05):
    z_array += [z]*20

ncomp = 3
exp_array = [500.]*len(z_array)

C_array = []
P_array = []
jacobian_array = []

uncor = 1

for i in range(len(z_array)):
    true_mag = random.normal()*0.1 + random.normal()*0.055*z_array[i] + random.normal()*2.1714724*0.001/z_array[i]
    true_EBV = 100.
    while true_EBV > 2.:
        true_EBV = random.exponential()
    true_EBV -= log(2.)
    true_EBV *= 0.1 + 0.05*z_array[i]
    true_projs = random.normal(size = ncomp)


    if uncor:
        errs = [0.03, 0.02] + [0.3]*ncomp
        P_array.append(concatenate(([true_mag, true_EBV], true_projs)) + random.normal(size = ncomp+2)*errs)
        C_array.append(diag(errs)**2.)
    else:
        rho_mat = diag(ones(ncomp + 2.))
        rho_mat[0,1] = -0.9
        rho_mat[1,0] = -0.9

        rho_mat[0,2:] = 0.2
        rho_mat[2:,0] = 0.2
        rho_mat[1,2:] = -0.2
        rho_mat[2:,1] = -0.2

        errs = [0.025, 0.1] + [3.]*ncomp
        
        Cmat = rho_mat*0.
        for m in range(len(Cmat)):
            for n in range(len(Cmat)):
                Cmat[m,n] = rho_mat[m,n]*errs[m]*errs[n]
        Lmat = linalg.cholesky(Cmat)
        P_array.append(concatenate(([true_mag, true_EBV], true_projs)) + dot(Lmat, random.normal(size = ncomp+2)))
        C_array.append(Cmat)
    jacobian_array.append(zeros([ncomp+2., 2]))

pickle.dump([P_array, C_array, z_array, exp_array, jacobian_array], open("sim_data.pickle", 'wb'))
