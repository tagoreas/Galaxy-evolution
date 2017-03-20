import numpy as np
import emcee
from mpi4py import MPI
import os, sys
from ctypes import *
import time


# numpy seed
np.random.seed(11)
runcount = 0

sampler = None

# link to C++ library
thisdir    = os.path.abspath(os.path.dirname(__file__))
libgalevol = cdll.LoadLibrary(os.path.join(thisdir, 'galevol/libgalevol.so'))

galevol_init1 = libgalevol.galevol_init1
galevol_init1.argtypes   = [c_void_p]
galevol_init1.restype    = None

galevol_init = libgalevol.galevol_init
galevol_init.argtypes   = [c_void_p]
galevol_init.restype    = c_void_p

galevol_sample = libgalevol.galevol_sample
galevol_sample.argtypes = [c_void_p, c_void_p, np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')]
galevol_sample.restype  = c_double

galevol_killme = libgalevol.galevol_killme
galevol_killme.argtypes = None
galevol_killme.restype  = None


# parallelization

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

galevolstruct = galevol_init1(MPI._addressof(comm))
comm.barrier()
galevolstruct = galevol_init(MPI._addressof(comm))
comm.barrier()

masterranks = np.loadtxt('master-ranks.dat')
if len(masterranks.shape)==0:
    masterranks = np.array([masterranks])

# only run one rank per node
# pthreads will make use of all CPU cores
if not rank in masterranks:
    print ('Rank',rank,'is sleeping.')
    if 0!=rank:
        while 1:
            time.sleep(10000)

print ('rank,size =',rank,size)
sys.stdout.flush()
if 0 != rank:
    status = MPI.Status()
    while True:
        # block and receive something
        task = comm.recv (source=0, tag=MPI.ANY_TAG, status=status)
        if task[0]==0:
            break
        # evaluate parameters
        if task[0]==1:
            arr = task[1::].astype (np.float64)
            ans = galevol_sample (MPI._addressof(comm), galevolstruct, arr)
            comm.isend (np.array([ans]), dest=0, tag=rank)
    sys.exit (0)


# function and enum decl.
def mcmcfunc(x):

    global runcount, galevolstruct, comm

    # signal to all threads to start computing
    arr = np.append (np.array([1]), np.array(x))
    
    requests = []
    for i in masterranks: #range (1,size):
        r = comm.isend (arr, dest=i, tag=i)
        requests.append (r)
    MPI.Request.waitall (requests)

    # get answers from all threads
    ans = 0;
    for i in masterranks:
        ans += comm.recv (source=i, tag=i)

    if np.isnan(ans): ans = -1e50

    if True or runcount%100==0:
            print (time.ctime(),'--',runcount,ans)
            sys.stdout.flush()


    # save chains along the way
    # this is very hackey. replace with proper code
    if (runcount+1)%nwalkers==0:        
        np.savetxt('mcmc.dat',np.hstack((sampler.flatchain,-2.0*np.transpose(np.array([sampler.flatlnprobability])))))
        halved=sampler.chain[::,sampler.chain.shape[1]/2:,::]
        np.savetxt('mcmc-half.dat',halved.reshape((halved.shape[0]*halved.shape[1],halved.shape[2])))
        halved=sampler.chain[::,sampler.chain.shape[1]*4/5:,::]
        np.savetxt('mcmc-fifth.dat',halved.reshape((halved.shape[0]*halved.shape[1],halved.shape[2])))


    runcount +=1
    return ans
    

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


# lens galaxy parameters
xiind      = enum ('z_i', 'm_star_i', 'r_e_i', 'num_xi')
etaind     = enum ('m_dm5_i', 'gamma_dm_i', 'alpha_imf_i', 'num_eta')

# galaxy population model parameters
psiind     = enum ('zeta_star', 'mu_star_0', 'sigma_star', 'zeta_r', \
                   'beta_r', 'mu_r_0', 'sigma_r', 'num_psi')
thetaind   = enum ('zeta_dm', 'beta_dm', 'xi_dm', 'm_dm_0', \
                   'sigma_dm', 'gamma_0', 'sigma_gamma', 'zeta_imf', \
                   'beta_imf', 'xi_imf', 'alpha_imf_0', \
                   'sigma_imf', 'llambda', 'num_theta')

# selection function
llambdaind = enum ('r_sel', 'sigma_sel', 'num_llambda')

# lens catalogs stuff
catind     = enum ('slacs', 'sl2s', 'num_cat')
catfuncind = enum ('mu_star', 'mu_r', 'num_catfunc')

# data format (observables)
dataind    = enum ('r_ein_i', 'sigma_i', 'm_star_i', 'r_eff_i', 'z_i', 'zsrc_i', 'Dod', 'Sigma_cr', 'num_observ')

guess_psi_slacs     = np.zeros (7)
guess_theta         = np.zeros (12)
guess_llambda_slacs = np.zeros (2)

guess_psi_slacs[psiind.mu_star_0]   = 11.5
guess_psi_slacs[psiind.zeta_star]   = 1.5
guess_psi_slacs[psiind.sigma_star]  = -0.75
guess_psi_slacs[psiind.mu_r_0]      = -0.3
guess_psi_slacs[psiind.zeta_r]      = -0.79
guess_psi_slacs[psiind.beta_r]      = 0.22
guess_psi_slacs[psiind.sigma_r]     = -0.75

guess_theta[thetaind.gamma_0]      = 0.80
guess_theta[thetaind.sigma_gamma]  = 0.34
guess_theta[thetaind.zeta_dm]      = 0
guess_theta[thetaind.beta_dm]      = 0
guess_theta[thetaind.xi_dm]        = 0
guess_theta[thetaind.m_dm_0]       = 11
guess_theta[thetaind.sigma_dm]     = -0.5
guess_theta[thetaind.zeta_imf]     = 0
guess_theta[thetaind.beta_imf]     = 0
guess_theta[thetaind.xi_imf]       = 0
guess_theta[thetaind.alpha_imf_0]  = 0
guess_theta[thetaind.sigma_imf]    = -1

guess_llambda_slacs[llambdaind.r_sel]     = np.log10(1.5)
guess_llambda_slacs[llambdaind.sigma_sel] = 0


ndim, nwalkers = 19, 100
p00 = np.array([ \
        guess_psi_slacs[psiind.mu_star_0], \
            guess_psi_slacs[psiind.zeta_star], \
            guess_psi_slacs[psiind.sigma_star], \
            guess_psi_slacs[psiind.mu_r_0], \
            guess_psi_slacs[psiind.zeta_r], \
            guess_psi_slacs[psiind.beta_r], \
            guess_psi_slacs[psiind.sigma_r], \
            guess_theta[thetaind.zeta_dm], \
            guess_theta[thetaind.beta_dm], \
            guess_theta[thetaind.xi_dm], \
            guess_theta[thetaind.m_dm_0], \
            guess_theta[thetaind.sigma_dm], \
            guess_theta[thetaind.zeta_imf], \
            guess_theta[thetaind.beta_imf], \
            guess_theta[thetaind.xi_imf], \
            guess_theta[thetaind.alpha_imf_0], \
            guess_theta[thetaind.sigma_imf], \
            guess_llambda_slacs[llambdaind.r_sel], \
            guess_llambda_slacs[llambdaind.sigma_sel] \
            ])


mcmcfunc(p00)

print ("Startup model:")
print (p00)
print ('Starting MCMC')
sys.stdout.flush()

iszero = (p00==0).astype('int')
randss = np.array([1, 5.5, 0.75, 2.5, 5.75, 7.25, 0.75, 6.25, 4.75, 2, 1, 1, 2.25, 3, 3, 1, 1, 1, 1])
p00 = np.array([11.5, 4.5, -0.75, -1.5, -0.75, -2.75, -0.75, 3.75, -1.25, 0, 11, -0.5, -0.25, 1, 2, 0, -1, 0, 0])
p0 = [ p00 + (np.random.rand(ndim)-0.5)*2*randss for i in range(nwalkers)]

nruns=2000

sampler = emcee.EnsembleSampler(nwalkers, ndim, mcmcfunc)

print ('start time: ',time.ctime())
sys.stdout.flush()


sampler.run_mcmc(p0, nruns)

np.savetxt('mcmc.dat',np.hstack((sampler.flatchain,-2.0*np.transpose(np.array([sampler.flatlnprobability])))))
halved=sampler.chain[::,sampler.chain.shape[1]/2:,::]
np.savetxt('mcmc-half.dat',halved.reshape((halved.shape[0]*halved.shape[1],halved.shape[2])))
halved=sampler.chain[::,sampler.chain.shape[1]*4/5:,::]
np.savetxt('mcmc-fifth.dat',halved.reshape((halved.shape[0]*halved.shape[1],halved.shape[2])))


print ('end time: ',time.ctime())
sys.stdout.flush()

# close all threads
for i in masterranks: #range(1,size):
    comm.send (np.array([0]), dest=i)
    
time.sleep(10)
galevol_killme()
