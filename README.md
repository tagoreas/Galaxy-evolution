# Galaxy-evolution
Constrain galaxy evolution using strong lenses

**What does galevol do?**

galevol takes observations of strong lenses and uses them to constrain the evolution of early-type galaxies.
Specifically, it constrains the evolution of the stellar mass, dark matter mass, and the IMF function.
In a hierarchical Bayesian framework, each individual lens is modeled but marginalized over (i.e., N integrals over the individual lens parameters are calculated, where N is the number of lenses).
The parameters of interest are sampled using emcee (an affine-invariant MCMC sampler).

**Requirements**

- python-emcee
- libopenblas.so

The MCMC sampling is done in Python. The emcee module is required.
The calculations are done in C++, and MPI compilers are required (mpic++).
Additionally, an optimized BLAS implementation is required. We use openblas by default.

**Installation**
To compile the C++ code, enter the galevol directory. There is a compile.sh script which should work if the above requirements are met.

For the first run, you will need to tabulate Einstein radii for your data set. Set the tabulate_rein_bool variables in galevol/galevol_main.cpp and alevol/galevol_calcs.cpp to 1. Compile and run.
Reset the  variables to 0, compile, and run. The newly tabulated Einstein radii will be used.
In the future, this step will hopefully be automated.

**Input Data**

An example of input data is included in the data directory. 
There is one lens per line.
Each line has 12 entries, corresponsding to 
1.  Lens redshift
2.  Lens redshift (1-sigma uncertainty)
3.  Source redshift
4.  Source redshift (1-sigma uncertainty)
5.  Lens effective radius
6.  Lens effective radius (1-sigma uncertainty)
7.  Lens Einstein radius
8.  Lens Einstein radius (1-sigma uncertainty)
9.  Lens velocity dispersion
10. Lens velocity dispersion (1-sigma uncertainty)
11. Lens stellar mass
12. Lens stellar mass (1-sigma uncertainty)

**Possible Issues**

galevol relies on two binary files, which contain pre-tabulated velocity dispersions and cosmological distance measurements.
These are loaded directly into C++ double arrays and assume 64-bit, Little Endian byte order.
