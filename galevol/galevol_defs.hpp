#ifndef GALEVOL_DEFS_HPP_
#define GALEVOL_DEFS_HPP_

#include <cmath>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/****************************************************/
/******************** CONSTANTS *********************/
/****************************************************/

#define GEDEFS galevol_defs::
class galevol_defs
{
public:
    
    const static double pi, sq2, sq2pi, twopi, fourpi,
	rad2arc, arc2rad, c2by4pig, G, pc2m,
	qags_ftol, mstarlim[2], mdm5lim[2], huge_penalty;
};


/****************************************************/
/************** GAL STRUCTS/ TYPEDEFS ***************/
/****************************************************/

// galaxy population parms
typedef struct
{
    int fatalerror;
    int *numlenses;
  //  double ***data, ***dataerr, ***dataorig, ***datatablower, ***datatabupper;
  double ***dataerr, ***dataorig;
    double **psi, *theta;
    void   ***catfunc;
    double **lambda;
  double *vdlut, *erlut, erlutints[5];
  double vdlutlims[9], erlutlims[15];
  size_t erlutsize[5],erlutsizemul[4];
  int vdlutsize[3];
    double *ADD;
  int numhosts;
  int *hostranks;
  char *hostnames;
  //  gsl_spline2d ****erlut_spline;
  //  gsl_interp_accel ****erlut_xacc;
} galevolstruct;

// individual galaxy parms
typedef struct
{
  //  double ***data, ***dataerr, ***dataorig, ***datatablower, ***datatabupper;
  double ***data, ***datatablower, ***datatabupper;
  int galind, catind;
  double integral_norms[4];
  double *eta, *xi;
  double *msnorms;
} galmodelstruct;

// galaxy parms
typedef struct
{
  galevolstruct *ges;
  galmodelstruct *gms;
  //double *rein_tab;
  //void *norm_spline, *norm_acc;
} galstruct;

// function pointer to catalog-specific functions
typedef double (*ge_catfunc_ptr) (void*);

// struct for multithreading
struct threadstruct
{
    int catind, startlens, numlenses, numthreads;
    double integral;
    galevolstruct *gestruct;
};

/****************************************************/
/******************** INDEXING **********************/
/****************************************************/

// lens galaxy parameters
enum xiind    {xiind_z_i, xiind_m_star_i, xiind_r_e_i, xiind_num_xi};
enum etaind   {etaind_m_dm5_i, etaind_gamma_dm_i, etaind_alpha_imf_i, etaind_num_eta};

// galaxy population model parameters
enum psiind   {psiind_zeta_star, psiind_mu_star_0, psiind_sigma_star,
	       psiind_zeta_r, psiind_beta_r, psiind_mu_r_0,
	       psiind_sigma_r, psiind_num_psi};
enum thetaind {thetaind_zeta_dm, thetaind_beta_dm, thetaind_xi_dm,
	       thetaind_m_dm_0, thetaind_sigma_dm,
	       thetaind_gamma_0, thetaind_sigma_gamma, thetaind_zeta_imf,
	       thetaind_beta_imf,
	       thetaind_xi_imf, thetaind_alpha_imf_0,
	       thetaind_sigma_imf, thetaind_llambda, thetaind_num_theta};

// selection function
enum lambdaind {lambdaind_r_sel, lambdaind_sigma_sel, lambdaind_num_lambda};

// lens catalogs
//enum catind     {catind_slacs, catind_sl2s, catind_euclid, catind_num_cat};
enum catind     {catind_slacs, catind_num_cat};
enum catfuncind {catfuncind_mu_star, catfuncind_mu_r, catfuncind_num_catfunc};

// data format (observables + some cosmology stuff)
enum dataind {dataind_z_i, dataind_zsrc_i, dataind_r_eff_i,
	      dataind_r_ein_i, dataind_sigma_i, dataind_m_star_i,
	      dataind_Dod, dataind_Sigma_cr, dataind_log10_r_eff_i, 
	      dataind_log10_r_eff_physical_i, dataind_num_observ};


/****************************************************/
/**************** GSL MINIMIZATION ******************/
/****************************************************/

struct gsl_mdm
{
    gsl_mdm() : T(gsl_multimin_fminimizer_nmsimplex2rand),
		s(NULL), ss(NULL), iter(0) {}

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_vector *ss;
    gsl_multimin_function minex_func;
    size_t iter;
    int status;
    double size;
};

#endif
