#include <cstring>
#include <algorithm>
#include <gsl/gsl_multifit.h>
#include <sstream>
#include <ctime>
#include <sys/time.h>
#include <typeinfo>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include "galevol_calcs.hpp"
#include "galevol_funcs_templates.cpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "fastonebigheader.h"

#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX  1
#include <mpi.h>

//const static int num_zzr         = 1;
//const static int num_mm          = 1001;
const static int num_zzr         = 7;
const static int num_mm          = 51;
const static size_t calls        = 10000;
const static size_t calls_xi_int = 10000;
const static int debugmode       = 0;

const static int num_zzr_m1      = num_zzr-1;
const static int num_zzr_m2      = num_zzr-2;
const static int num_mm_m1       = num_mm-1;
const static int num_mm_m2       = num_mm-2;

// global variables for tabulating Einstein radii (only needs to happen once)
// tabulation expects 1 node with 16 cores.
// One core is the master node and does not calculate R_Ein, so there are 15 working nodes
int tabulate_rein_bool = 0; // tabulate or don't
int rein_block_index=0; // tracks which parameters a particular node loops over
pthread_mutex_t rein_mutex; // for assigning rein_block_index

// BLAS matrix vector multiplication
extern "C" void dgemv_ ( const char *transa, int *m, int *n, double *alpha,
                                             const double *a, int *lda, double*x, int *incx,
                                             double *beta, double *y, int *incy);


double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double mdm5_norm_interp (void *parms_)
{
    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *datatabupper= gms->datatabupper[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    double *theta   = ges->theta;
    void **catfunc  = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    double *erlut   = ges->erlut;
    //
    // END setup variables
    //
    
    
    
    double intnorm;
    double *msnorms = gms->msnorms;
    double reints[4] = { (datatabupper[dataind_z_i]    -datatablower[dataind_z_i])    /num_zzr_m1, 
			 (datatabupper[dataind_zsrc_i] -datatablower[dataind_zsrc_i]) /num_zzr_m1, 
			 (datatabupper[dataind_r_eff_i]-datatablower[dataind_r_eff_i])/num_zzr_m1,
			 2.0/num_mm_m1 };
    double refracs[4] = { (data[dataind_z_i]     -datatablower[dataind_z_i])     /reints[0], 
			  (data[dataind_zsrc_i]  -datatablower[dataind_zsrc_i])  /reints[1], 
			  (data[dataind_r_eff_i] -datatablower[dataind_r_eff_i]) /reints[2],
			  (mstar-GEDEFS mstarlim[0])/reints[3] };
    int reinds[4] = {(int)refracs[0],(int)refracs[1],(int)refracs[2],(int)refracs[3]};
    double wt1t = refracs[0]-reinds[0];
    double wt2t = refracs[1]-reinds[1];
    double wt3t = refracs[2]-reinds[2];
    double wt4t = refracs[3]-reinds[3];
    double wt1[2] = {1.0-wt1t,wt1t};
    double wt2[2] = {1.0-wt2t,wt2t};
    double wt3[2] = {1.0-wt3t,wt3t};
    double wt4[2] = {1.0-wt4t,wt4t};
    double avg = 0;
    int msnormind, inds[3];
    double weights[3];
    if (1==num_zzr)
      {
	for (int g=0; g<2; ++g)
	  {
	    msnormind = (reinds[3]+g);
	    double val = msnorms[msnormind];
	    if (val<std::log(1e-300*2)+1e-5) return std::log(1e-300);
	    else avg += wt4[g] *val;
	  }
      }
    else
      {
	  if (0)
	  {
	      for (int m=0; m<2; ++m)
	      {
		  inds[0] = (reinds[0]+m)*(num_zzr*num_zzr*num_mm);
		  weights[0] = wt1[m];
		  for (int k=0; k<2; ++k)
		  {
		      inds[1] = inds[0] + (reinds[1]+k)*(num_zzr*num_mm);
		      weights[1] = weights[0] *wt2[k];
		      for (int f=0; f<2; ++f)
		      {
			  inds[2] = inds[1] + (reinds[2]+f)*(num_mm);
			  weights[2] = weights[1] *wt3[f];
			  for (int g=0; g<2; ++g)
			  {
			      msnormind = inds[2] + (reinds[3]+g);
			      double val = msnorms[msnormind];
			      if (val<std::log(1e-300*2)+1e-5) return std::log(1e-300);
			      else avg += weights[2] *wt4[g] *val;
			  }
		      }
		  }
	      }
	  }    
	  else
	  {
	      // This interpolation uses a higher-dimensional implentation
	      // of the thin-plate spline
	      // It uses the bounding 2^4 points plus 4 points close to the corner of the n-dimensional box
	      double valsarr[21+4], distarr[21+4];
	      const static double *kinv=0;
	      int valind=0;	  
	      int cornerind=0;
	      int cinds[4] = {0,0,0,0};
	      int signs[4] = {-1,-1,-1,-1};
	      if ( ( (1-wt1[0])>0.5 && reinds[0]!=num_zzr_m2 ) || 0==reinds[0]) {cornerind += 8; cinds[0]=1; signs[0]=1;}
	      if ( ( (1-wt2[0])>0.5 && reinds[1]!=num_zzr_m2 ) || 0==reinds[1]) {cornerind += 4; cinds[1]=1; signs[1]=1;}
	      if ( ( (1-wt3[0])>0.5 && reinds[2]!=num_zzr_m2 ) || 0==reinds[2]) {cornerind += 2; cinds[2]=1; signs[2]=1;}
	      if ( ( (1-wt4[0])>0.5 && reinds[3]!=num_mm_m2  ) || 0==reinds[3]) {cornerind += 1; cinds[3]=1; signs[3]=1;}
#include "kinv_mdm5_1.cpp"
#include "kinv_mdm5_2.cpp"
#include "kinv_mdm5_3.cpp"
#include "kinv_mdm5_4.cpp"
#include "kinv_mdm5_5.cpp"
#include "kinv_mdm5_6.cpp"
#include "kinv_mdm5_7.cpp"
#include "kinv_mdm5_8.cpp"
#include "kinv_mdm5_9.cpp"
#include "kinv_mdm5_10.cpp"
#include "kinv_mdm5_11.cpp"
#include "kinv_mdm5_12.cpp"
#include "kinv_mdm5_13.cpp"
#include "kinv_mdm5_14.cpp"
#include "kinv_mdm5_15.cpp"
#include "kinv_mdm5_16.cpp"
	      switch (cornerind)
	      {
	      case 0:  kinv = kinv1;  break;
	      case 1:  kinv = kinv2;  break;
	      case 2:  kinv = kinv3;  break;
	      case 3:  kinv = kinv4;  break;
	      case 4:  kinv = kinv5;  break;
	      case 5:  kinv = kinv6;  break;
	      case 6:  kinv = kinv7;  break;
	      case 7:  kinv = kinv8;  break;
	      case 8:  kinv = kinv9;  break;
	      case 9:  kinv = kinv10; break;
	      case 10: kinv = kinv11; break;
	      case 11: kinv = kinv12; break;
	      case 12: kinv = kinv13; break;
	      case 13: kinv = kinv14; break;
	      case 14: kinv = kinv15; break;
	      case 15: kinv = kinv16; break;
	      default: std::cerr << "No Way!!" << std::endl; exit(0);
	      }
	      
	      for (int m=0; m<2; ++m)
	      {
		  inds[0] = (reinds[0]+m)*(num_zzr*num_zzr*num_mm);
		  weights[0] = (1-wt1[m])*(1-wt1[m]);
		  for (int k=0; k<2; ++k)
		  {
		      inds[1] = inds[0] + (reinds[1]+k)*(num_zzr*num_mm);
		      weights[1] = weights[0] + (1-wt2[k])*(1-wt2[k]);
		      for (int f=0; f<2; ++f)
		      {
			  inds[2] = inds[1] + (reinds[2]+f)*(num_mm);
			  weights[2] = weights[1] + (1-wt3[f])*(1-wt3[f]);
			  for (int g=0; g<2; ++g)
			  {
			      msnormind = inds[2] + (reinds[3]+g);
			      double val = msnorms[msnormind];
			      if (val<std::log(1e-300*2)+1e-5) return std::log(1e-300);
			      valsarr[valind] = val;
			      distarr[valind] = weights[2] + (1-wt4[g])*(1-wt4[g]);
			      if (cornerind==valind)
			      {
				  int msnormindjj[4] = { (reinds[0]+m+signs[0])*(num_zzr*num_zzr*num_mm) + (reinds[1]+k)*(num_zzr*num_mm) + (reinds[2]+f)*(num_mm) + (reinds[3]+g),
							    (reinds[0]+m)*(num_zzr*num_zzr*num_mm) + (reinds[1]+k+signs[1])*(num_zzr*num_mm) + (reinds[2]+f)*(num_mm) + (reinds[3]+g),
							    (reinds[0]+m)*(num_zzr*num_zzr*num_mm) + (reinds[1]+k)*(num_zzr*num_mm) + (reinds[2]+f+signs[2])*(num_mm) + (reinds[3]+g),
							    (reinds[0]+m)*(num_zzr*num_zzr*num_mm) + (reinds[1]+k)*(num_zzr*num_mm) + (reinds[2]+f)*(num_mm) + (reinds[3]+g+signs[3]) };
				  distarr[16+0] = distarr[valind] - (1-wt1[cinds[0]])*(1-wt1[cinds[0]]) + (1-wt1[cinds[0]]+signs[0])*(1-wt1[cinds[0]]+signs[0]);
				  distarr[16+1] = distarr[valind] - (1-wt2[cinds[1]])*(1-wt2[cinds[1]]) + (1-wt2[cinds[1]]+signs[1])*(1-wt2[cinds[1]]+signs[1]);
				  distarr[16+2] = distarr[valind] - (1-wt3[cinds[2]])*(1-wt3[cinds[2]]) + (1-wt3[cinds[2]]+signs[2])*(1-wt3[cinds[2]]+signs[2]);
				  distarr[16+3] = distarr[valind] - (1-wt4[cinds[3]])*(1-wt4[cinds[3]]) + (1-wt4[cinds[3]]+signs[3])*(1-wt4[cinds[3]]+signs[3]);
				  for (int jj=0; jj<4; ++jj)
				  {
				      val = msnorms[msnormindjj[jj]];
				      if (val<std::log(1e-300*2)+1e-5) return std::log(1e-300);
				      valsarr[16+jj] = val;
				  }
			      }
			      ++valind;
			  }
		      }
		  }
	      }
	      	      
	      valsarr[16+4] = valsarr[17+4] = valsarr[18+4] = valsarr[19+4] = valsarr[20+4] = 0;
	      char tr = 'N';	  
	      double one = 1.0;	  
	      double zer = 0.0;	  
	      int onei = 1;	  
	      int tnrow = 21+4;
	      int tncol = 21+4;
	      double weightsarr[21+4];
	      dgemv_ (&tr, &tnrow, &tncol, &one, kinv, &tnrow,
		      valsarr, &onei, &zer, weightsarr, &onei);	      
	      intnorm = weightsarr[16+4] + weightsarr[17+4]*wt1t + 
		  weightsarr[18+4]*wt2t +weightsarr[19+4]*wt3t +weightsarr[20+4]*wt4t;
	      for (int jj=0; jj<16+4; ++jj)
		  intnorm += weightsarr[jj]*distarr[jj];
	  }
      }
	      
    intnorm = avg;
    return intnorm;
}

double GECALCS vel_dispersion (void *parms_)
{
    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *datatabupper= gms->datatabupper[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    double *theta   = ges->theta;
    void **catfunc  = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    double *erlut   = ges->erlut;
    double *vdlutlims = ges->vdlutlims;
    int *vdlutsize = ges->vdlutsize;
    double *erlutlims = ges->erlutlims;
    //
    // END setup variables
    //

    double reff = data[dataind_log10_r_eff_physical_i] + std::log10(5.0);
    
    double reints[3] = { (vdlutlims[3]-vdlutlims[0])/(vdlutlims[6]-1), 
			 (vdlutlims[4]-vdlutlims[1])/(vdlutlims[7]-1),
			 (vdlutlims[5]-vdlutlims[2])/(vdlutlims[8]-1) };
    double refracs[3] = { (reff  -vdlutlims[0]) /reints[0],
			  (mstar -vdlutlims[1]) /reints[1],
			  (mdm5  -vdlutlims[2]) /reints[2] };
    int reinds[3] = {(int)refracs[0],(int)refracs[1],(int)refracs[2]};
    double wt1t = refracs[0]-reinds[0];
    double wt2t = refracs[1]-reinds[1];
    double wt3t = refracs[2]-reinds[2];
    double wt1[2] = {1.0-wt1t,wt1t};
    double wt2[2] = {1.0-wt2t,wt2t};
    double wt3[2] = {1.0-wt3t,wt3t};
    double avg = 0;
    int inds[2];
    double weights[2];
    if (0)
    {
	for (int m=0; m<2; ++m)
	{
	    inds[0] = (reinds[0]+m)*vdlutsize[1]*vdlutsize[2];
	    weights[0] = wt1[m];
	    for (int n=0; n<2; ++n)
	    {
		inds[1] = inds[0] + (reinds[1]+n)*vdlutsize[2];
		weights[1] = weights[0] *wt2[n];
		for (int k=0; k<2; ++k)
		{
		    avg += weights[1] *wt3[k] *vdlut[inds[1]+(reinds[2]+k)];
		}
	    }
	}
    }
    else
    {
	double valsarr[12], distarr[12];
	int valind=0;
	for (int m=0; m<2; ++m)
	{
	    inds[0] = (reinds[0]+m)*vdlutsize[1]*vdlutsize[2];
	    weights[0] = (1-wt1[m])*(1-wt1[m]);
	    for (int n=0; n<2; ++n)
	    {
		inds[1] = inds[0] + (reinds[1]+n)*vdlutsize[2];
		weights[1] = weights[0] + (1-wt2[n])*(1-wt2[n]);
		for (int k=0; k<2; ++k)
		{
		    distarr[valind] = weights[1] + (1-wt3[k])*(1-wt3[k]);
		    valsarr[valind] = vdlut[inds[1]+(reinds[2]+k)];
		    ++valind;
		}
	    }
	}
#include "kinv_veldisp.cpp"
	valsarr[8] = valsarr[9] = valsarr[10] = valsarr[11] = 0;
	char tr = 'N';
	double one = 1.0;	
	double zer = 0.0;	
	int onei = 1;	
	int tnrow = 12;	
	int tncol = 12;	
	double weightsarr[12];	
	dgemv_ (&tr, &tnrow, &tncol, &one, kinv, &tnrow,
		valsarr, &onei, &zer, weightsarr, &onei);	
	avg = weightsarr[8] + weightsarr[9]*wt1t +
	    weightsarr[10]*wt2t +weightsarr[11]*wt3t;
	for (int jj=0; jj<8; ++jj)
	    avg += weightsarr[jj]*distarr[jj];
    }
    
    return avg;
}

double ein_rad_interp (void *parms_)
{
    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *datatabupper= gms->datatabupper[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    double *theta   = ges->theta;
    void **catfunc  = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    double *erlut   = ges->erlut;
    double *vdlutlims = ges->vdlutlims;
    double *erlutlims = ges->erlutlims;
    size_t *erlutsize = ges->erlutsize;
    size_t *erlutsizemul = ges->erlutsizemul;
    double *erlutints = ges->erlutints;
    //
    // END setup variables
    //
    
    if (debugmode) return data[dataind_r_eff_i];

    double refracs[5] = { (data[dataind_z_i]     -erlutlims[0]) /erlutints[0],
			  (data[dataind_zsrc_i]  -erlutlims[1]) /erlutints[1],
			  (data[dataind_r_eff_i] -erlutlims[2]) /erlutints[2],
			  (mstar                 -erlutlims[3]) /erlutints[3],
			  (mdm5                  -erlutlims[4]) /erlutints[4]  };
    size_t reinds[5] = {(size_t)refracs[0],(size_t)refracs[1],(size_t)refracs[2],
			(size_t)refracs[3],(size_t)refracs[4]};
    double wt1t = refracs[0]-reinds[0];
    double wt2t = refracs[1]-reinds[1];
    double wt3t = refracs[2]-reinds[2];
    double wt4t = refracs[3]-reinds[3];
    double wt5t = refracs[4]-reinds[4];
    double wt1[2] = {1.0-wt1t,wt1t};
    double wt2[2] = {1.0-wt2t,wt2t};
    double wt3[2] = {1.0-wt3t,wt3t};
    double wt4[2] = {1.0-wt4t,wt4t};
    double wt5[2] = {1.0-wt5t,wt5t};
    double avg = 0;
    size_t inds[4];
    double weights[4];
    for (size_t m=0; m<2; ++m)
      {
	inds[0] = (reinds[0]+m)*erlutsizemul[0];
	weights[0] = wt1[m];
      for (size_t n=0; n<2; ++n)
	{
	  inds[1] = inds[0] + (reinds[1]+n)*erlutsizemul[1];
	  weights[1] = weights[0] *wt2[n];
	for (size_t k=0; k<2; ++k)
	  {
	    inds[2] = inds[1] + (reinds[2]+k)*erlutsizemul[2];
	    weights[2] = weights[1] *wt3[k];
	  for (size_t f=0; f<2; ++f)
	    {
	      inds[3] = inds[2] + (reinds[3]+f)*erlutsizemul[3];
	      weights[3] = weights[2] *wt4[f];
	    for (size_t g=0; g<2; ++g)
	      {
	      avg += weights[3] *wt5[g]	*erlut[inds[3]+(reinds[4]+g)];
	      }
	    }
	  }
	}
      }

    return avg;
}

double ein_rad_interp_4mdm5norm (void *parms_, size_t reinds[5], double wt1[2], double wt2[2], double wt3[2], double wt4[2], double wt5[2])
{
    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *datatabupper= gms->datatabupper[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    double *theta   = ges->theta;
    void **catfunc  = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    double *erlut   = ges->erlut;
    double *vdlutlims = ges->vdlutlims;
    double *erlutlims = ges->erlutlims;
    size_t *erlutsize = ges->erlutsize;
    size_t *erlutsizemul = ges->erlutsizemul;
    //const static size_t erlutsizemul[5] = {(size_t)76*(size_t)76*(size_t)76*(size_t)76,(size_t)76*(size_t)76*(size_t)76,(size_t)76*(size_t)76,(size_t)76};
    double *erlutints = ges->erlutints;
    //
    // END setup variables
    //

    if (debugmode) return data[dataind_r_eff_i];

    double refrac5 =  (mdm5-erlutlims[4]) /erlutints[4];
    size_t reind5  = (size_t)refrac5;
    double wt5t    = refrac5-reind5;
    wt5[0] = 1.0-wt5t;
    wt5[1] = wt5t;
    double avg = 0;
    size_t inds[4];
    double weights[4];
    for (size_t m=0; m<2; ++m)
      {
	inds[0] = (reinds[0]+m)*erlutsizemul[0];
	weights[0] = wt1[m];
      for (size_t n=0; n<2; ++n)
	{
	  inds[1] = inds[0] + (reinds[1]+n)*erlutsizemul[1];
	  weights[1] = weights[0]*wt2[n];
	for (size_t k=0; k<2; ++k)
	  {
	    inds[2] = inds[1] + (reinds[2]+k)*erlutsizemul[2];
	    weights[2] = weights[1]*wt3[k];
	  for (size_t f=0; f<2; ++f)
	    {
	      inds[3] = inds[2] + (reinds[3]+f)*erlutsizemul[3];
	      weights[3] = weights[2]*wt4[f];
	    for (size_t g=0; g<2; ++g)
	      {
		avg += weights[3] *wt5[g] *erlut[inds[3]+(reind5+g)];
	      }
	    }
	  }
	}
      }
    return avg;
}

double GECALCS ein_radius_min_func (const gsl_vector *v, void *parms_)
{
    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *theta   = ges->theta;
    void **catfunc = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    //
    // END setup variables
    //

    double re, reff, Dod, sigmacr, xr3, xr5;

    // exponentiate log masses
    mstar = std::pow (10.0, mstar);
    mdm5  = std::pow (10.0, mdm5);
    
    re      = gsl_vector_get (v, 0);
    // reject small and negative angles a priori, and also for numerical stability
    if (re<1e-4)
        return 1e10*(1+re*re);
    Dod     = data[dataind_Dod];
    reff    = data[dataind_r_eff_i];
    sigmacr = data[dataind_Sigma_cr];

    // de Vaucouleurs
    double k, kappa0, zfac, pararg, remodel;
    k       = 7.66925001;
    kappa0  = std::pow (k,8.0) /(40320.0 *GEDEFS pi *reff*reff) *mstar;
    kappa0  = kappa0 /sigmacr;
    zfac    = k *std::pow (re/reff, 0.25);
    pararg  = 1.0;
    for (int t=7; t>=1; --t)
        pararg = pararg * zfac/t+1.0;
    pararg  = 1.0 - pararg *std::exp (-zfac);
    remodel = kappa0 *40320.0 *std::pow (k,-8.0) *reff*reff /re *pararg;
        
    // DM halo
    double rs, xr, xr2, lim, rs2, lim2, fac1, fac2, thisre;
    rs      = 10.0 *reff;
    xr      = re/rs;
    xr2     = xr*xr;
    lim     = 5.0 /Dod *GEDEFS rad2arc;
    rs2     = rs*rs;
    lim2    = lim*lim;
    fac1    = std::sqrt (std::abs (rs2-lim2));
    fac2    = std::sqrt (std::abs (1.0-xr2));

    // get mass normalization
    if (lim<rs)
        kappa0  = mdm5 /
	    (GEDEFS twopi *rs2 *(
		2.0*rs*atanh(fac1/rs)/fac1 + std::log(lim2/(4.0*rs2))
		));
    else
        kappa0 = mdm5 /
	    (GEDEFS fourpi *rs2 *(
		rs*atan(fac1/rs)/fac1 + std::log(lim/(2.0*rs))
		));
    kappa0  = kappa0 /sigmacr;

    // get deflection
    if (std::abs(xr-1.0)<1e-3)
    {
	// treat near unity xr as unity
	thisre = 4.0 *kappa0 *rs *(std::log(0.5*1.0)+1.0) /1.0;
    }
    else if (xr<0.1)
    {
	// use series expansion for small r (more numerically stable)
	xr3 = xr2*xr;
	xr5 = xr3*xr2;
	thisre = 4.0 *kappa0 *rs
	    *(0.25 *(-1.0+2.0 *std::log(2.0)-2.0 *std::log(xr))*xr +
	      0.03125 *(-7.0+12.0*std::log(2.0)-12.0*std::log(xr)) *xr3+
	      0.0052083333333333333333 *(-37.0+60.0 *std::log(2.0)-60.0
					 *std::log(xr)) *xr5);
    }
    else if (xr<1.0)
    {
	// use exact formula for other small xr values
	thisre = 4.0 *kappa0 *rs *(std::log(0.5*xr)+atanh(fac2)/fac2) /xr;
    }
    else
    {
	// xr>1
        thisre = 4.0 *kappa0 *rs *(std::log(0.5*xr)+atan(fac2)/fac2) /xr;
    }
    remodel += thisre;

    return (remodel-re)*(remodel-re);
}

double GECALCS ein_radius (void *parms_)
{
    double einradius;
    int maxiter = 1000;
    double tolerance = 0.00001;
    gsl_vector *x  = gsl_vector_alloc (1);
    gsl_vector *ss = gsl_vector_alloc (1);
    gsl_vector_set (x, 0, 10);
    gsl_vector_set (ss, 0, 2);

    struct gsl_mdm mdm;
    mdm.minex_func.n      = 1;
    mdm.minex_func.f      = &GECALCS ein_radius_min_func;
    mdm.minex_func.params = parms_;
    mdm.s                 = gsl_multimin_fminimizer_alloc (mdm.T, mdm.minex_func.n);
    gsl_multimin_fminimizer_set (mdm.s, &mdm.minex_func, x, ss);
    do
    {
	if ((int)mdm.iter>maxiter)
	{
	    std::cerr << "Einstein radius minimization failed. Too many iterations." << std::endl;
	    
	    galstruct *parms    = (galstruct*)parms_;
	    galevolstruct *ges  = parms->ges;
	    galmodelstruct *gms = parms->gms;
	    int catind      = gms->catind;
	    int galind      = gms->galind;
	    double mstar    = gms->xi[xiind_m_star_i];
	    double mdm5     = gms->eta[etaind_m_dm5_i];
	    double gammadm  = gms->eta[etaind_gamma_dm_i];
	    double *data    = gms->data[catind][galind];
	    double *dataerr = ges->dataerr[catind][galind];
	    double *theta   = ges->theta;
	    void **catfunc = ges->catfunc[catind];
	    double *psi     = ges->psi[catind];
	    double *lambda  = ges->lambda[catind];
	    double *vdlut   = ges->vdlut;
	    double re, reff, Dod, sigmacr, xr3, xr5;
	    Dod     = data[dataind_Dod];
	    reff    = data[dataind_r_eff_i];
	    sigmacr = data[dataind_Sigma_cr];
	    std::cerr << mstar << " " << mdm5 << " " << reff << " " << Dod << " " << sigmacr << std::endl;
	    
	    break;
	}

	++mdm.iter;
	mdm.status = gsl_multimin_fminimizer_iterate (mdm.s);

	if (mdm.status)
	    break;

	mdm.size   = gsl_multimin_fminimizer_size (mdm.s);
	mdm.status = gsl_multimin_test_size (mdm.size, tolerance);
    }
    while (mdm.status == GSL_CONTINUE);

    // copy best value
    einradius = gsl_vector_get (gsl_multimin_fminimizer_x (mdm.s), 0);

    // cleanup
    gsl_multimin_fminimizer_free (mdm.s);
    gsl_vector_free (x);
    gsl_vector_free (ss);

    return einradius;
}

typedef struct 
{
  void *parms_;
  size_t reinds[5];
  double wt1[2];
  double wt2[2];
  double wt3[2];
  double wt4[2];
  double wt5[2];
  int numeval;
} mdm5intstruct;
double GECALCS integral_mdm5_norm (double var, void *parms__)
{
  void *parms_   = ((mdm5intstruct*)parms__)->parms_;
  size_t *reinds = ((mdm5intstruct*)parms__)->reinds;
  double *wt1   = ((mdm5intstruct*)parms__)->wt1;
  double *wt2   = ((mdm5intstruct*)parms__)->wt2;
  double *wt3   = ((mdm5intstruct*)parms__)->wt3;
  double *wt4   = ((mdm5intstruct*)parms__)->wt4;
  double *wt5   = ((mdm5intstruct*)parms__)->wt5;
  //((mdm5intstruct*)parms__)->numeval++;

    // set mdm5
    ((galstruct*)parms_)->gms->eta[etaind_m_dm5_i] = var;
    
    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *theta   = ges->theta;
    void **catfunc  = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    double *erlut   = ges->erlut;
    //
    // END setup variables
    //

    double mudm, re, int2, int3, Dod;
    double facs[2], args[2];
    double product;

    Dod     = data[dataind_Dod];
    mudm = theta[thetaind_zeta_dm] *(data[dataind_z_i]-0.7)  +
      theta[thetaind_beta_dm] *(mstar-11.5) +
      theta[thetaind_xi_dm]   *(mstar-11.5 -
				2.0 *data[dataind_log10_r_eff_physical_i]) +
      theta[thetaind_m_dm_0];
    
    re = ein_rad_interp_4mdm5norm (parms_, reinds, wt1, wt2, wt3, wt4, wt5);
    
    if (0)
      {
	// logmdm5 term
	facs[0] = std::log (1.0 /(GEDEFS sq2pi *theta[thetaind_sigma_dm]));
	args[0] = (-(mdm5-mudm)*(mdm5-mudm)
		   /(2.0*theta[thetaind_sigma_dm]*theta[thetaind_sigma_dm]));
	
	// selection function
	facs[1] = std::log (re*re);
	args[1] = (-(re-lambda[lambdaind_r_sel])*(re-lambda[lambdaind_r_sel])
		   /(2.0*lambda[lambdaind_sigma_sel]*lambda[lambdaind_sigma_sel]));
	
	product = std::exp (facs[0]+facs[1]+args[0]+args[1]);
	product = product<1e-300 ? 1e-300 : product;
      }
    else
      {    
	// faster evaluation
	facs[0] = re*re /(GEDEFS sq2pi *theta[thetaind_sigma_dm]);
	args[0] = (-(mdm5-mudm)*(mdm5-mudm)/(2.0*theta[thetaind_sigma_dm]*theta[thetaind_sigma_dm])) + 
	  (-(re-lambda[lambdaind_r_sel])*(re-lambda[lambdaind_r_sel])
	   /(2.0*lambda[lambdaind_sigma_sel]*lambda[lambdaind_sigma_sel]));
	
	product = 0;
	if (args[0]>-87.33654480436835770919)  // if fast float exp works
	  product = facs[0]*fastexp(args[0]);
	else if (args[0]>-745)                 // if double exp works
	  product = facs[0]*std::exp(args[0]);
	else if (facs[0]>1)                    // if factor will increase 'product'
	  {
	    double newarg = fastlog(facs[0]) + args[0];
	    if (newarg>-87.33654480436835770919) // if fast float exp works
	      product = fastexp(args[0]);
	    else if (newarg>-745)                // if double exp works
	      product = std::exp(args[0]);
	  }
	product = product<1e-300 ? 1e-300 : product;
      }
    
    return product;
}

double xiintnorm (double *vars, size_t dimensions, void *parms_)
{
  {
    // set vars
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;
    int catind      = gms->catind;
    int galind      = gms->galind;
    double *data    = gms->data[catind][galind];

    gms->data[catind][galind][dataind_z_i]     = vars[0];
    gms->data[catind][galind][dataind_r_eff_i] = vars[1];
    gms->xi[xiind_m_star_i]    = vars[2];

    double Dod = ges->ADD[(int)(data[dataind_z_i]/6.5*99999)];
    double Dos = ges->ADD[(int)(data[dataind_zsrc_i]/6.5*99999)];
    double Dds = Dos - (1+data[dataind_z_i]) /(1+data[dataind_zsrc_i]) *Dod;
    data[dataind_Dod] = Dod;
    data[dataind_Sigma_cr] = 39084.4403*Dos*Dod/Dds;
    data[dataind_log10_r_eff_i] = std::log10(data[dataind_r_eff_i]);
    data[dataind_log10_r_eff_physical_i] = data[dataind_log10_r_eff_i] + std::log10(GEDEFS arc2rad*Dod/5.0);
  }

    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *theta   = ges->theta;
    void **catfunc  = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    double *vdlutlims   = ges->vdlutlims;
    double *erlut   = ges->erlut;
    double *datatabupper= gms->datatabupper[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    //
    // END setup variables
    //

    double mur, mustar, fac, arg;
    mustar = psi[psiind_zeta_star] *(data[dataind_z_i]-0.7) +psi[psiind_mu_star_0];
    mur    = psi[psiind_zeta_r] *(data[dataind_z_i]-0.7) + psi[psiind_beta_r]
      *(mstar-11.5) + psi[psiind_mu_r_0];
    fac    = std::log (1.0 /(GEDEFS twopi *psi[psiind_sigma_star] *psi[psiind_sigma_r]));
    arg    = ((-(mstar-mustar)*(mstar-mustar)
	       /(2*psi[psiind_sigma_star]*psi[psiind_sigma_star])
	       -(data[dataind_log10_r_eff_i]-mur)
	       *(data[dataind_log10_r_eff_i]-mur)
	       /(2*psi[psiind_sigma_r]*psi[psiind_sigma_r])));
    return std::exp(fac+arg);    
}

double int2dfunc (double *vars, size_t dimensions, void *parms_)
{

  {
    // set vars
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;
    int catind      = gms->catind;
    int galind      = gms->galind;

    if (ges->fatalerror)
      return 0.00;

    if (1==num_zzr)
      {
	gms->xi[xiind_m_star_i]    = vars[0];
	gms->eta[etaind_m_dm5_i]   = vars[1];
      }
    else
      {
	gms->data[catind][galind][dataind_z_i]     = vars[0];
	gms->data[catind][galind][dataind_zsrc_i]  = vars[1];
	gms->data[catind][galind][dataind_log10_r_eff_i] = vars[2];
	gms->xi[xiind_m_star_i]    = vars[3];
	gms->eta[etaind_m_dm5_i]   = vars[4];
      }
  }
    
    //
    // START setup variables
    //
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;

    int catind      = gms->catind;
    int galind      = gms->galind;
    double mstar    = gms->xi[xiind_m_star_i];
    double mdm5     = gms->eta[etaind_m_dm5_i];
    double gammadm  = gms->eta[etaind_gamma_dm_i];
    double *data    = gms->data[catind][galind];
    double *integral_norms = gms->integral_norms;
    double *dataorig= ges->dataorig[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *theta   = ges->theta;
    void **catfunc  = ges->catfunc[catind];
    double *psi     = ges->psi[catind];
    double *lambda  = ges->lambda[catind];
    double *vdlut   = ges->vdlut;
    double *vdlutlims   = ges->vdlutlims;
    double *erlut   = ges->erlut;
    double *datatabupper= gms->datatabupper[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    //
    // END setup variables
    //

    if (num_zzr!=1)
    {
      double Dod = ges->ADD[(int)(data[dataind_z_i]/6.5*99999)];
      double Dos = ges->ADD[(int)(data[dataind_zsrc_i]/6.5*99999)];
      double Dds = Dos - (1+data[dataind_z_i]) /(1+data[dataind_zsrc_i]) *Dod;
      data[dataind_Dod] = Dod;
      data[dataind_Sigma_cr] = 39084.4403*Dos*Dod/Dds;
      data[dataind_r_eff_i] = std::pow(10.0,data[dataind_log10_r_eff_i]);
      data[dataind_log10_r_eff_physical_i] = data[dataind_log10_r_eff_i] + std::log10(GEDEFS arc2rad*Dod/5.0);
      
      if (data[dataind_log10_r_eff_physical_i]+(std::log10(5.0)-1e-8)<vdlutlims[0]) return 1e-300;
      if (data[dataind_log10_r_eff_physical_i]+(std::log10(5.0)+1e-8)>vdlutlims[3]) return 1e-300;
    }

    // get mdm5 norm
    double intnorm = mdm5_norm_interp (parms_);
    if (intnorm</*std::log(1e-300*2)+1e-5*/std::log(1e-290)) return 1e-300;

    // get mstar norm
    double mustar = psi[psiind_zeta_star] *(data[dataind_z_i]-0.7) +psi[psiind_mu_star_0];

    double mur, mudm, lmstarerr, simf, muimf, re, vd, Dod;

    Dod     = data[dataind_Dod];
    mur    = psi[psiind_zeta_r] *(data[dataind_z_i]-0.7) + psi[psiind_beta_r]
      *(mstar-11.5) + psi[psiind_mu_r_0];
    mudm = theta[thetaind_zeta_dm] *(data[dataind_z_i]-0.7)  +
	theta[thetaind_beta_dm] *(mstar-11.5) +
	theta[thetaind_xi_dm]   *(mstar-11.5 -
				  2.0 *data[dataind_log10_r_eff_physical_i]) +
	theta[thetaind_m_dm_0];

    lmstarerr = dataerr[dataind_m_star_i];
    simf      = theta[thetaind_sigma_imf];
    muimf     = theta[thetaind_zeta_imf] *(data[dataind_z_i]-0.7)  +
	theta[thetaind_beta_imf] *(mstar-11.5) +
	theta[thetaind_xi_imf]   *(mstar-11.5 -
				   2.0 *data[dataind_log10_r_eff_physical_i]) +
	theta[thetaind_alpha_imf_0];
    
    
    re = ein_rad_interp (parms_);
    vd = GECALCS vel_dispersion (parms_);
    if (ges->fatalerror)
	return 0.00;
    
    double reff_norm, zs_norm, zd_norm, mstar_norm;
    zd_norm    = integral_norms[0];
    zs_norm    = integral_norms[1];
    reff_norm  = integral_norms[2];
    mstar_norm = integral_norms[3];
    
    double facs[5], args[5], product;
   
    if (0)
      {
	facs[0] = std::log (1.0 /(GEDEFS twopi *psi[psiind_sigma_star] *psi[psiind_sigma_r]) *mstar_norm);
	args[0] = (-(mstar-mustar)*(mstar-mustar)
		   /(2*psi[psiind_sigma_star]*psi[psiind_sigma_star])
		   -(data[dataind_log10_r_eff_i]-mur)
		   *(data[dataind_log10_r_eff_i]-mur)
		   /(2*psi[psiind_sigma_r]*psi[psiind_sigma_r]));
	
	
	facs[1] = std::log (1.0 /(GEDEFS sq2pi *theta[thetaind_sigma_dm]));
	args[1] = (-(mdm5-mudm)*(mdm5-mudm)
		   /(2*theta[thetaind_sigma_dm]*theta[thetaind_sigma_dm]));
	
	// selection function
	facs[2] = std::log (re*re);
	args[2] = (-(re-lambda[lambdaind_r_sel])*(re-lambda[lambdaind_r_sel])
		   /(2.0*lambda[lambdaind_sigma_sel]*lambda[lambdaind_sigma_sel]));
	
	// observational uncertainties
	if (1==num_zzr)
	{
	    facs[3] = std::log (1.0 /(GEDEFS twopi *dataerr[dataind_r_ein_i] *dataerr[dataind_sigma_i]));
	    args[3] = (-(data[dataind_r_ein_i]-re)*(data[dataind_r_ein_i]-re)
		   /(2*dataerr[dataind_r_ein_i]*dataerr[dataind_r_ein_i])
		   -(data[dataind_sigma_i]-vd)*(data[dataind_sigma_i]-vd)
		       /(2*dataerr[dataind_sigma_i]*dataerr[dataind_sigma_i]));
	}	
	else
	{
	    facs[3] = std::log (1.0 /(GEDEFS twopi *dataerr[dataind_r_ein_i] *dataerr[dataind_sigma_i]))
	  + std::log (1.0 /(GEDEFS twopi *dataerr[dataind_r_eff_i] *dataerr[dataind_z_i])*reff_norm*zd_norm)
	  + std::log (1.0 /(GEDEFS sq2pi *dataerr[dataind_zsrc_i])*zs_norm);
	args[3] = (-(data[dataind_r_ein_i]-re)*(data[dataind_r_ein_i]-re)
		   /(2*dataerr[dataind_r_ein_i]*dataerr[dataind_r_ein_i])
		   -(data[dataind_sigma_i]-vd)*(data[dataind_sigma_i]-vd)
		   /(2*dataerr[dataind_sigma_i]*dataerr[dataind_sigma_i])
		   -(data[dataind_r_eff_i]-dataorig[dataind_r_eff_i])*(data[dataind_r_eff_i]-dataorig[dataind_r_eff_i])
		   /(2*dataerr[dataind_r_eff_i]*dataerr[dataind_r_eff_i])
		   -(data[dataind_z_i]-dataorig[dataind_z_i])*(data[dataind_z_i]-dataorig[dataind_z_i])
		   /(2*dataerr[dataind_z_i]*dataerr[dataind_z_i])
		   -(data[dataind_zsrc_i]-dataorig[dataind_zsrc_i])*(data[dataind_zsrc_i]-dataorig[dataind_zsrc_i])
		   /(2*dataerr[dataind_zsrc_i]*dataerr[dataind_zsrc_i]));
	}
	
	
	
	// log alpha integral
	facs[4] = std::log (1.0 /(GEDEFS sq2pi *std::sqrt (lmstarerr*lmstarerr + simf*simf)));
	args[4] = (-(data[dataind_m_star_i]-mstar+muimf)*(data[dataind_m_star_i]-mstar+muimf)
		   /(2 *(lmstarerr*lmstarerr + simf*simf)));
	
	product = std::exp (facs[0]+facs[1]+facs[2]+facs[3]+facs[4]+
			    args[0]+args[1]+args[2]+args[3]+args[4]-intnorm);
	if (intnorm<std::log(1e-300*2)+1e-5) product = 0;
	product = product<1e-300 ? 1e-300 : product;
      }
    else
      {
	  if (1==num_zzr)
	  {
	      facs[0] = std::log
	  (
	   (1.0 /(GEDEFS twopi *psi[psiind_sigma_star] *psi[psiind_sigma_r]) *mstar_norm) *
	   (1.0 /(GEDEFS sq2pi *theta[thetaind_sigma_dm])) *
	   (re*re) *
	   (1.0 /(GEDEFS twopi *dataerr[dataind_r_ein_i] *dataerr[dataind_sigma_i])) *
	   (1.0 /(GEDEFS sq2pi *std::sqrt (lmstarerr*lmstarerr + simf*simf)))
	   );
	args[0] = ((-(mstar-mustar)*(mstar-mustar)
		    /(2*psi[psiind_sigma_star]*psi[psiind_sigma_star])
		    -(data[dataind_log10_r_eff_i]-mur)
		    *(data[dataind_log10_r_eff_i]-mur)
		    /(2*psi[psiind_sigma_r]*psi[psiind_sigma_r]))) +
	  ((-(mdm5-mudm)*(mdm5-mudm)
	    /(2*theta[thetaind_sigma_dm]*theta[thetaind_sigma_dm]))) +
	  ((-(re-lambda[lambdaind_r_sel])*(re-lambda[lambdaind_r_sel])
	    /(2.0*lambda[lambdaind_sigma_sel]*lambda[lambdaind_sigma_sel]))) +
	  ((-(data[dataind_r_ein_i]-re)*(data[dataind_r_ein_i]-re)
	    /(2*dataerr[dataind_r_ein_i]*dataerr[dataind_r_ein_i])
	    -(data[dataind_sigma_i]-vd)*(data[dataind_sigma_i]-vd)
	    /(2*dataerr[dataind_sigma_i]*dataerr[dataind_sigma_i]))) +
	  ((-(data[dataind_m_star_i]-mstar+muimf)*(data[dataind_m_star_i]-mstar+muimf)
	    /(2 *(lmstarerr*lmstarerr + simf*simf))));
	  }
	  else
	  {
	      facs[0] = std::log
	  (
	   (1.0 /(GEDEFS twopi *psi[psiind_sigma_star] *psi[psiind_sigma_r]) *mstar_norm) *
	   (1.0 /(GEDEFS sq2pi *theta[thetaind_sigma_dm])) *
	   (re*re) *
	   (1.0 /(GEDEFS twopi *dataerr[dataind_r_ein_i] *dataerr[dataind_sigma_i])) *
	   (1.0 /(GEDEFS twopi *dataerr[dataind_r_eff_i] *dataerr[dataind_z_i])*reff_norm*zd_norm) *
	   (1.0 /(GEDEFS sq2pi *dataerr[dataind_zsrc_i])*zs_norm) *
	   (1.0 /(GEDEFS sq2pi *std::sqrt (lmstarerr*lmstarerr + simf*simf)))
	   );
	args[0] = ((-(mstar-mustar)*(mstar-mustar)
		    /(2*psi[psiind_sigma_star]*psi[psiind_sigma_star])
		    -(data[dataind_log10_r_eff_i]-mur)
		    *(data[dataind_log10_r_eff_i]-mur)
		    /(2*psi[psiind_sigma_r]*psi[psiind_sigma_r]))) +
	  ((-(mdm5-mudm)*(mdm5-mudm)
	    /(2*theta[thetaind_sigma_dm]*theta[thetaind_sigma_dm]))) +
	  ((-(re-lambda[lambdaind_r_sel])*(re-lambda[lambdaind_r_sel])
	    /(2.0*lambda[lambdaind_sigma_sel]*lambda[lambdaind_sigma_sel]))) +
	  ((-(data[dataind_r_ein_i]-re)*(data[dataind_r_ein_i]-re)
	    /(2*dataerr[dataind_r_ein_i]*dataerr[dataind_r_ein_i])
	    -(data[dataind_sigma_i]-vd)*(data[dataind_sigma_i]-vd)
	    /(2*dataerr[dataind_sigma_i]*dataerr[dataind_sigma_i])
	    -(data[dataind_r_eff_i]-dataorig[dataind_r_eff_i])*(data[dataind_r_eff_i]-dataorig[dataind_r_eff_i])
	    /(2*dataerr[dataind_r_eff_i]*dataerr[dataind_r_eff_i])
	    -(data[dataind_z_i]-dataorig[dataind_z_i])*(data[dataind_z_i]-dataorig[dataind_z_i])
	    /(2*dataerr[dataind_z_i]*dataerr[dataind_z_i])
	    -(data[dataind_zsrc_i]-dataorig[dataind_zsrc_i])*(data[dataind_zsrc_i]-dataorig[dataind_zsrc_i])
	    /(2*dataerr[dataind_zsrc_i]*dataerr[dataind_zsrc_i]))) +
	  ((-(data[dataind_m_star_i]-mstar+muimf)*(data[dataind_m_star_i]-mstar+muimf)
	    /(2 *(lmstarerr*lmstarerr + simf*simf))));
	  }
	  
	product = std::exp (facs[0]+args[0]-intnorm);
	product = product<1e-300 ? 1e-300 : product;
      }
    
    
    //std::cout << "parms " << mstar << " " << mdm5 << " " << gammadm << " " << std::endl;
    //std::cout << "facargs ";
    //for(int hh=0; hh<5; ++hh)std::cout << facs[hh] << " " << args[hh] << " ";
    //std::cout << intnorm << " " << product << std::endl;
    //exit(0);  
    
    return product;
}

void tabulate_r_ein (void *parms_)
{
  int myreinind;
  pthread_mutex_lock (&rein_mutex);
  myreinind = rein_block_index++; 
  if (myreinind>=15) return;
  std::cout << "found index " << myreinind << std::endl;
  pthread_mutex_unlock (&rein_mutex);
  
    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;
    int catind      = gms->catind;
    int galind      = gms->galind;
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    double *datatabupper= gms->datatabupper[catind][galind];
    double zd = data[dataind_z_i];
    double zs = data[dataind_zsrc_i];
    double reff = data[dataind_r_eff_i];
    double zderr = dataerr[dataind_z_i];
    double zserr = dataerr[dataind_zsrc_i];
    double refferr = dataerr[dataind_r_eff_i];
   
    int status;
    double integral, intnorm, res, err;

    // setup integration limits
    double xl[5],xu[5];
    xl[0] = 0.05;       xu[0] = 1.65;
    xl[1] = 0.30;       xu[1] = 6.45;
    xl[2] = 0.05;       xu[2] = 6.15;
    xl[3] = GEDEFS mstarlim[0];  xu[3] = GEDEFS mstarlim[1];
    xl[4] = GEDEFS mdm5lim[0];   xu[4] = GEDEFS mdm5lim[1];
    datatablower[dataind_z_i]     = xl[0];
    datatablower[dataind_zsrc_i]  = xl[1];
    datatablower[dataind_r_eff_i] = xl[2];
    datatabupper[dataind_z_i]     = xu[0];
    datatabupper[dataind_zsrc_i]  = xu[1];
    datatabupper[dataind_r_eff_i] = xu[2];

    // set variables from outer integral
    data[dataind_z_i]     = dataorig[dataind_z_i];
    data[dataind_zsrc_i]  = dataorig[dataind_zsrc_i];
    data[dataind_r_eff_i] = dataorig[dataind_r_eff_i];
    zd   = dataorig[dataind_z_i];
    zs   = dataorig[dataind_zsrc_i];
    reff = dataorig[dataind_r_eff_i];
    double Dod = ges->ADD[(int)(data[dataind_z_i]/6.5*99999)];
    double Dos = ges->ADD[(int)(data[dataind_zsrc_i]/6.5*99999)];
    double Dds = Dos - (1+data[dataind_z_i]) /(1+data[dataind_zsrc_i]) *Dod;
    data[dataind_Dod] = Dod;
    data[dataind_Sigma_cr] = 39084.4403*Dos*Dod/Dds;
    data[dataind_log10_r_eff_i] = std::log10(data[dataind_r_eff_i]);
    data[dataind_log10_r_eff_physical_i] = std::log10(data[dataind_r_eff_i]*GEDEFS arc2rad*Dod/5.0);

    // tabulating Einstein radius
    int mystartind = 5*(0+myreinind);
    int mystopind = 14==myreinind ? 76 : 5*(1+myreinind);
    double *reintab;
    size_t totalsize = mystopind-mystartind;
    for (int j=0; j<4; ++j) totalsize *= (size_t)76;
    GEFUNCS ge_malloc (&reintab,totalsize);
    std::fill (reintab,reintab+totalsize,0);
    {
    double thisms, thismd, thiszs, thiszd, thisreff;
    double inter1=(xu[0]-xl[0])/75., inter2=(xu[1]-xl[1])/75., inter3=(xu[2]-xl[2])/75., inter4=2./75., inter5=2./75.;

    size_t erindex = 0;
    for (int xx=0; xx<76; ++xx) 
      {
	for (int yy=0; yy<76; ++yy) 
	  {
	    for (int aa=0; aa<76; ++aa) 
	      {
		for (int bb=0; bb<76; ++bb) 
		  {
		    for (int cc=mystartind; cc<mystopind; ++cc) 
		      {
			thiszd   = xl[0] + xx*inter1;
			thiszs   = xl[1] + yy*inter2;
			thisreff = xl[2] + aa*inter3;
			thisms   = GEDEFS mstarlim[0] + bb*inter4;
			thismd   = GEDEFS mdm5lim[0]  + cc*inter5;
			if (75 ==xx) thiszd   = xu[0];
			if (75 ==yy) thiszs   = xu[1];
			if (75 ==aa) thisreff = xu[2];
			if (75 ==bb) thisms   = GEDEFS mstarlim[1];
			if (75 ==cc) thismd   = GEDEFS mdm5lim[1];
			
			gms->xi[xiind_m_star_i]  = thisms;
			gms->eta[etaind_m_dm5_i] = thismd;
			data[dataind_z_i]        = thiszd;
			data[dataind_zsrc_i]     = thiszs;
			data[dataind_r_eff_i]    = thisreff;
			double Dod = ges->ADD[(int)(data[dataind_z_i]/6.5*99999)];
			double Dos = ges->ADD[(int)(data[dataind_zsrc_i]/6.5*99999)];
			double Dds = Dos - (1+data[dataind_z_i]) /(1+data[dataind_zsrc_i]) *Dod;
			data[dataind_Dod] = Dod;
			data[dataind_Sigma_cr] = 39084.4403*Dos*Dod/Dds;
			data[dataind_log10_r_eff_i] = std::log10(data[dataind_r_eff_i]);
			data[dataind_log10_r_eff_physical_i] = std::log10(data[dataind_r_eff_i]*GEDEFS arc2rad*Dod/5.0);
			if (thiszd>thiszs-0.05)
			  reintab[erindex++] = -1;
			else
			  reintab[erindex++] = GECALCS ein_radius (parms_);
		      }
		  }
	      }
	  }
	string Result  = GEFUNCS ge_tostring(myreinind);
	string Result2 = GEFUNCS ge_tostring(xx);
	string mystr = "erlut/erlut."+Result+"."+Result2+".bin";
	std::cout << "writing " << mystr << std::endl;
	FILE *file = fopen (mystr.c_str(), "ab");
	fwrite ((void*)(reintab+(totalsize/76)*xx), sizeof(reintab[0]), totalsize/76, file);
	fclose (file);
      }
    }
}

double GECALCS integrate_one_lens_allerr (void *parms_)
{
    if (((galstruct*)parms_)->ges->fatalerror)
	return 0.00;

    galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;
    galmodelstruct *gms = parms->gms;
    int catind      = gms->catind;
    int galind      = gms->galind;
    double *data    = gms->data[catind][galind];
    double *dataerr = ges->dataerr[catind][galind];
    double *dataorig= ges->dataorig[catind][galind];
    double *datatablower= gms->datatablower[catind][galind];
    double *datatabupper= gms->datatabupper[catind][galind];
    double zd = data[dataind_z_i];
    double zs = data[dataind_zsrc_i];
    double reff = data[dataind_r_eff_i];
    double zderr = dataerr[dataind_z_i];
    double zserr = dataerr[dataind_zsrc_i];
    double refferr = dataerr[dataind_r_eff_i];
    double *vdlutlims = ges->vdlutlims;
    double *erlutlims = ges->erlutlims;
    size_t *erlutsize = ges->erlutsize;
    size_t *erlutsizemul = ges->erlutsizemul;
    double *erlutints = ges->erlutints;
    double *integral_norms = gms->integral_norms;
    double *psi     = ges->psi[catind];

    int status;
    double integral, intnorm, res, err;

    if (1==num_zzr)
      {
	double Dod = ges->ADD[(int)(data[dataind_z_i]/6.5*99999)];
	double Dos = ges->ADD[(int)(data[dataind_zsrc_i]/6.5*99999)];
	double Dds = Dos - (1+data[dataind_z_i]) /(1+data[dataind_zsrc_i]) *Dod;
	data[dataind_Dod] = Dod;
	data[dataind_Sigma_cr] = 39084.4403*Dos*Dod/Dds;
	data[dataind_log10_r_eff_i] = std::log10(data[dataind_r_eff_i]);
	data[dataind_log10_r_eff_physical_i] = data[dataind_log10_r_eff_i] + std::log10(GEDEFS arc2rad*Dod/5.0);
      }
    
    // setup integration limits
    double xl[5],xu[5];
    xl[0] = zd   -3*zderr;       xu[0] = zd   +3*zderr;
    xl[1] = zs   -3*zserr;       xu[1] = zs   +3*zserr;
    xl[2] = reff -3*refferr;     xu[2] = reff +3*refferr;
    xl[3] = GEDEFS mstarlim[0];  xu[3] = GEDEFS mstarlim[1];
    xl[4] = GEDEFS mdm5lim[0];   xu[4] = GEDEFS mdm5lim[1];
    if (xl[0]<0.05)  xl[0] = 0.05;
    if (xu[0]>xl[1]) {xu[0] = xl[1]; xl[1] +=0.1;}
    if (xu[0]<xl[0]) {std::cerr << "err: what?! 1" << std::endl;exit(1);}
    if (xu[1]<xl[1]) {std::cerr << "err: what?! 2" << std::endl;exit(1);}
    if (xl[2]<0.05)  xl[2] = 0.05;

    if (xl[0]<erlutlims[0]) xl[0] = erlutlims[0] + 1e-8;
    if (xl[1]<erlutlims[1]) xl[1] = erlutlims[1] + 1e-8;
    if (xl[2]<erlutlims[2]) xl[2] = erlutlims[2] + 1e-8;
    if (xu[0]>erlutlims[5]) xu[0] = erlutlims[5] - 1e-8;
    if (xu[1]>erlutlims[6]) xu[1] = erlutlims[6] - 1e-8;
    if (xu[2]>erlutlims[7]) xu[2] = erlutlims[7] - 1e-8;

    datatablower[dataind_z_i]     = xl[0];
    datatablower[dataind_zsrc_i]  = xl[1];
    datatablower[dataind_r_eff_i] = xl[2];
    datatabupper[dataind_z_i]     = xu[0];
    datatabupper[dataind_zsrc_i]  = xu[1];
    datatabupper[dataind_r_eff_i] = xu[2];


    // get integral norms
    double seconds0 = get_wall_time();
    {
      if (1!=num_zzr)
	{
	    /*
	  double xl2[3] = {xl[0],xl[2],xl[3]};
	  double xu2[3] = {xu[0],xu[2],xu[3]};
	  const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_monte_function G = { &xiintnorm, 3, parms_ };
	  //gsl_rng_env_setup ();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
	  int numiter = 0;
	  status = gsl_monte_vegas_integrate (&G, xl2, xu2, 3, calls_xi_int, r, s, &res, &err);
	  do
	    {
	      if (20==numiter) {status = 1010101; break;}
	      status += gsl_monte_vegas_integrate (&G, xl2, xu2, 3, calls_xi_int, r, s, &res, &err);
	      ++numiter;
	    }
	  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
	  gsl_monte_vegas_free (s);
	  gsl_rng_free (r);
	    */

	  double zd_norm,zs_norm,reff_norm;
	  zd_norm = 2.0/( gsl_sf_erf((datatabupper[dataind_z_i]-dataorig[dataind_z_i])
				     /(GEDEFS sq2*dataerr[dataind_z_i])) -
			  gsl_sf_erf((datatablower[dataind_z_i]-dataorig[dataind_z_i])
				     /(GEDEFS sq2*dataerr[dataind_z_i])) );
	  zs_norm = 2.0/( gsl_sf_erf((datatabupper[dataind_zsrc_i]-dataorig[dataind_zsrc_i])
				     /(GEDEFS sq2*dataerr[dataind_zsrc_i])) -
			  gsl_sf_erf((datatablower[dataind_zsrc_i]-dataorig[dataind_zsrc_i])
				     /(GEDEFS sq2*dataerr[dataind_zsrc_i])) );
	  reff_norm = 2.0/( gsl_sf_erf((datatabupper[dataind_r_eff_i]-dataorig[dataind_r_eff_i])
				       /(GEDEFS sq2*dataerr[dataind_r_eff_i])) -
			    gsl_sf_erf((datatablower[dataind_r_eff_i]-dataorig[dataind_r_eff_i])
				       /(GEDEFS sq2*dataerr[dataind_r_eff_i])) );
	  
	  double mustar = psi[psiind_zeta_star] *(data[dataind_z_i]-0.7) +psi[psiind_mu_star_0];
	  res = 2.0/( gsl_sf_erf((GEDEFS mstarlim[1]-11.5)/(GEDEFS sq2*psi[psiind_sigma_star])) -
		      gsl_sf_erf((GEDEFS mstarlim[0]-11.5)/(GEDEFS sq2*psi[psiind_sigma_star])) );

	  integral_norms[0] = zd_norm;
	  integral_norms[1] = zs_norm;
	  integral_norms[2] = reff_norm;
	  integral_norms[3] = res;
	}
      else
	{
	  double x2 = data[dataind_log10_r_eff_i];
	  double a = psi[psiind_zeta_r] *(data[dataind_z_i]-0.7) + psi[psiind_beta_r]*(-11.5) + psi[psiind_mu_r_0];
	  double b = psi[psiind_beta_r];
	  double s1 = psi[psiind_sigma_star];
	  double s2 = psi[psiind_sigma_r];
	  double x01 = psi[psiind_zeta_star] *(data[dataind_z_i]-0.7) +psi[psiind_mu_star_0];
	  double l1 = GEDEFS mstarlim[0];
	  double l2 = GEDEFS mstarlim[1];
	  double abs12 = a*b*s1*s1;
	  double b2l2s12 = b*b*l2*s1*s1;
	  double bs12x2 = b*s1*s1*x2;
	  double b2l1s12 = b*b*l1*s1*s1;
	  double b2s12 = b*b*s1*s1;
	  double s22 = s2*s2;
	  
	  /*
	  res = std::exp(-((a + b*x01 - x2)*(a + b*x01 - x2)/(2*(b2s12 + s22)))) *
	    (-gsl_sf_erf((abs12 + b2l1s12 - bs12x2 + s22*(l1 - x01))/(std::sqrt(2)*s1*s2*std::sqrt(b2s12 + s22))) +
	     gsl_sf_erf((abs12 + b2l2s12 - bs12x2 + s22*(l2 - x01))/(std::sqrt(2)*s1*s2*std::sqrt(b2s12 + s22))))/
	    (2*GEDEFS twopi*std::sqrt(b2s12 + s22));	  
	  */

	  double mustar = psi[psiind_zeta_star] *(data[dataind_z_i]-0.7) +psi[psiind_mu_star_0];
	  res = 2.0/( gsl_sf_erf((GEDEFS mstarlim[1]-11.5)/(GEDEFS sq2*psi[psiind_sigma_star])) -
		      gsl_sf_erf((GEDEFS mstarlim[0]-11.5)/(GEDEFS sq2*psi[psiind_sigma_star])) );
	  
	  integral_norms[0] = 1;
	  integral_norms[1] = 1;
	  integral_norms[2] = 1;
	  integral_norms[3] = res;
	}
    }
    

    // tabulate normalizations for mdm5 integrals
    double seconds1 = get_wall_time();
    {
      GEFUNCS ge_malloc (&(gms->msnorms), num_mm*num_zzr*num_zzr*num_zzr);
      double *normvals = gms->msnorms;
      mdm5intstruct mis;
      mis.parms_ = parms_;
      gsl_integration_workspace *w2 = gsl_integration_workspace_alloc (10000);
      double thisms, thismd, thiszs, thiszd, thisreff;
      double inter1=2./num_mm_m1, inter2=(xu[0]-xl[0])/num_zzr_m1, 
	inter3=(xu[1]-xl[1])/num_zzr_m1, inter4=(xu[2]-xl[2])/num_zzr_m1;
      size_t normind = 0;
      for (int aa=0; aa<num_zzr; ++aa) 
	{
	  thiszd   = xl[0] + aa*inter2;
	  if (num_zzr_m1 ==aa) thiszd   = xu[0]-1e-8;
	  if (1==num_zzr) thiszd = zd;
	  data[dataind_z_i]        = thiszd;
	  double thisrefrac=(thiszd-erlutlims[0])/erlutints[0]; 
	  size_t thisreind=(size_t)thisrefrac; double thiswtt=thisrefrac-thisreind;
	  mis.wt1[0]=1.0-thiswtt; mis.wt1[1]=thiswtt; mis.reinds[0]=thisreind;
	  double Dod = ges->ADD[(int)(thiszd/6.5*99999)];
	  for (int bb=0; bb<num_zzr; ++bb) 
	    {
	      thiszs   = xl[1] + bb*inter3;
	      if (num_zzr_m1 ==bb) thiszs   = xu[1]-1e-8;
	      if (1==num_zzr) thiszs = zs;
	      data[dataind_zsrc_i]     = thiszs;
	      double thisrefrac=(thiszs-erlutlims[1])/erlutints[1]; 
	      size_t thisreind=(size_t)thisrefrac; double thiswtt=thisrefrac-thisreind;
	      mis.wt2[0]=1.0-thiswtt; mis.wt2[1]=thiswtt; mis.reinds[1]=thisreind;
	      double Dos = ges->ADD[(int)(thiszs/6.5*99999)];
	      double Dds = Dos - (1+thiszd) /(1+thiszs) *Dod;
	      data[dataind_Dod] = Dod;
	      data[dataind_Sigma_cr] = 39084.4403*Dos*Dod/Dds;
	      for (int cc=0; cc<num_zzr; ++cc) 
		{
		  thisreff = xl[2] + cc*inter4;
		  if (num_zzr_m1 ==cc) thisreff = xu[2]-1e-8;
		  if (1==num_zzr) thisreff = reff;
		  data[dataind_r_eff_i]    = thisreff;
		  double thisrefrac=(thisreff-erlutlims[2])/erlutints[2]; 
		  size_t thisreind=(size_t)thisrefrac; double thiswtt=thisrefrac-thisreind;
		  mis.wt3[0]=1.0-thiswtt; mis.wt3[1]=thiswtt; mis.reinds[2]=thisreind;
		  data[dataind_log10_r_eff_i] = std::log10(data[dataind_r_eff_i]);
		  data[dataind_log10_r_eff_physical_i] = data[dataind_log10_r_eff_i] + 
		    std::log10(GEDEFS arc2rad*Dod/5.0);
		  for (int xx=0; xx<num_mm; ++xx) 
		    { 
		      thisms   = GEDEFS mstarlim[0] + xx*inter1;
		      if (num_mm_m1==xx)   thisms   = GEDEFS mstarlim[1]-1e-8;
		      gms->xi[xiind_m_star_i]  = thisms;
		      double thisrefrac=(thisms-erlutlims[3])/erlutints[3]; 
		      size_t thisreind=(size_t)thisrefrac; double thiswtt=thisrefrac-thisreind;
		      mis.wt4[0]=1.0-thiswtt; mis.wt4[1]=thiswtt; mis.reinds[3]=thisreind;
		      
		      mis.numeval=0;
		      gsl_function F2;
		      F2.function = &GECALCS integral_mdm5_norm;
		      F2.params = &mis;
		      status = gsl_integration_qags (&F2, GEDEFS mdm5lim[0], GEDEFS mdm5lim[1],
						     1e-280*0, 1e-3, 10000, w2, &intnorm, &err);

		      if (status || ((galstruct*)parms_)->ges->fatalerror)
			{
			  //std::cerr << "QAGS failed: mdm5 integral norm" << std::endl;
			  intnorm = 1e-300;
			}
		      normvals[normind++] = std::log(intnorm);
		    }
		}
	    }
	}
      gsl_integration_workspace_free (w2);
    }

    // do actual integration
    double seconds2=get_wall_time();
    {
      if (1==num_zzr)
	{
	  const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_monte_function G = { &int2dfunc, 2, parms_ };
	  //gsl_rng_env_setup ();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
	  int numiter = 0;
	  status = gsl_monte_vegas_integrate (&G, xl+3, xu+3, 2, calls, r, s, &res, &err);
	  do
	    {
	      if (20==numiter) {status = 1010101; break;}
	      status += gsl_monte_vegas_integrate (&G, xl+3, xu+3, 2, calls, r, s, &res, &err);
	      ++numiter;
	    }
	  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
	  gsl_monte_vegas_free (s);
	  gsl_rng_free (r);
	  integral = std::log (res);
	}
      else
	{
	    xl[2] = std::log10(xl[2]) + 1e-8;
	    xu[2] = std::log10(xu[2]) - 1e-8;
	  const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_monte_function G = { &int2dfunc, 5, parms_ };
	  //gsl_rng_env_setup ();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (5);
	  int numiter = 0;
	  status = gsl_monte_vegas_integrate (&G, xl, xu, 5, calls, r, s, &res, &err);
	  do
	    {
	      if (20==numiter) {status = 1010101; break;}
	      status += gsl_monte_vegas_integrate (&G, xl, xu, 5, calls, r, s, &res, &err);
	      ++numiter;
	    }
	  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
	  gsl_monte_vegas_free (s);
	  gsl_rng_free (r);
	  integral = std::log (res);
	}
    }
    double seconds3=get_wall_time();


    // free stuff
    GEFUNCS ge_free (gms->msnorms);

    if(0)
      std::cout << "eval times " << seconds1-seconds0 << " " << 
      seconds2-seconds1 << " " << seconds3-seconds2 << " " << 
      galind << " " << integral << std::endl;
    // error check
    if (status || ((galstruct*)parms_)->ges->fatalerror)
    {
      std::cerr << "QAGS failed: inner integral " << 
	status << " " << galind << std::endl;
	((galstruct*)parms_)->ges->fatalerror = 
	  0==((galstruct*)parms_)->ges->fatalerror ? 1 : ((galstruct*)parms_)->ges->fatalerror;
	integral = 0;
    }

    return integral;    
}

double GECALCS lnprobability (void *vcomm, void *parms_)
{
  MPI_Comm comm = *((MPI_Comm*)vcomm);
  int size, rank;
  char pname[MPI_MAX_PROCESSOR_NAME]; int len;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  MPI_Get_processor_name(pname, &len);
  pname[len] = 0;
  
  galstruct *parms    = (galstruct*)parms_;
    galevolstruct *ges  = parms->ges;

    gsl_set_error_handler_off();
    
    double integral = 0;
    ges->fatalerror = 0;

    // skip if this rank isn't spawning the pthreads and get hostindex
    int hostind = -1;
  char hostname[HOST_NAME_MAX];
  char username[LOGIN_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  getlogin_r(username, LOGIN_NAME_MAX);

  for (int j=0; j<ges->numhosts; ++j)
    if (!strcmp(ges->hostnames+j*100,hostname) && rank!=ges->hostranks[j])
      return 0;
    else if (!strcmp(ges->hostnames+j*100,hostname))
      {
	hostind = j;
	break;
      }

  // find out if this is the node with rank 0
  int ismasternode=0;
  if (!strcmp(hostname,ges->hostnames+ges->numhosts*100))
    ismasternode=1;

    for (int cat=0; cat<catind_num_cat; ++cat)
    {
      if (tabulate_rein_bool)
	pthread_mutex_init (&rein_mutex, 0);

	// setup multithreading and launch threads
      int numthreads = size/ges->numhosts;
      int interval = (int)(std::ceil((double)ges->numlenses[cat]/(size-1)))*numthreads;
      numthreads *=2;
      
	pthread_t    threads[numthreads];
	threadstruct tstruct[numthreads];
	pthread_attr_t attrjoinable;
	pthread_attr_init (&attrjoinable);
	pthread_attr_setdetachstate (&attrjoinable, PTHREAD_CREATE_JOINABLE);
	for (int thr=0; thr<numthreads; ++thr)
	  {
	    tstruct[thr].catind     = cat;
	    tstruct[thr].startlens  = hostind*interval+thr;//hostind*(ges->numlenses[cat]/ges->numhosts)+thr;//rank-1;//thr;
	    tstruct[thr].numlenses  = (hostind+1)*interval;//(hostind+1)*(ges->numlenses[cat]/ges->numhosts);//ges->numlenses[cat];
	    tstruct[thr].numthreads = numthreads;//size-1;//numthreads;
	    tstruct[thr].gestruct   = ges;

	    if (ges->numhosts-1==hostind)
	      tstruct[thr].numlenses  = ges->numlenses[cat];
	    if (tstruct[thr].numlenses>ges->numlenses[cat])
	      tstruct[thr].numlenses = ges->numlenses[cat];

	    pthread_create (&threads[thr], &attrjoinable,
	    	    &(GECALCS lnprob_thread), &tstruct[thr]);
	    //GECALCS lnprob_thread (&tstruct[thr]);
	}

	// wait for threads to finish
	for (int thr=0; thr<numthreads; ++thr)
	    pthread_join (threads[thr], NULL);
	
	if (tabulate_rein_bool)
	  exit(0);
	
	// add results from threads
	for (int thr=0; thr<numthreads; ++thr)
	{
	    integral += tstruct[thr].integral;
	}
    }

    if (ges->fatalerror)
	return GEDEFS huge_penalty;

    return integral;
}

void* GECALCS lnprob_thread (void *args_)
{
    // init
    threadstruct *args      = (threadstruct*)args_;
    args->integral = 0;
    
    if (args->startlens>=args->numlenses)
      {
	usleep(1000);
	return 0;
      }
    
    // create galaxy model struct
    galmodelstruct gmstruct;
    gmstruct.catind = args->catind;
    GEFUNCS ge_malloc (&gmstruct.xi,  xiind_num_xi);
    GEFUNCS ge_malloc (&gmstruct.eta, etaind_num_eta);

    // create global struct
    galstruct gstruct;
    gstruct.ges = args->gestruct;
    gstruct.gms = &gmstruct;

    GEFUNCS ge_malloc (&gmstruct.data,         catind_num_cat);
    GEFUNCS ge_malloc (&gmstruct.datatablower, catind_num_cat);
    GEFUNCS ge_malloc (&gmstruct.datatabupper, catind_num_cat);
    GEFUNCS ge_malloc (&gmstruct.data[args->catind],         gstruct.ges->numlenses[args->catind]);
    GEFUNCS ge_malloc (&gmstruct.datatablower[args->catind], gstruct.ges->numlenses[args->catind]);
    GEFUNCS ge_malloc (&gmstruct.datatabupper[args->catind], gstruct.ges->numlenses[args->catind]);
    std::fill (gmstruct.data[args->catind], gmstruct.data[args->catind]+gstruct.ges->numlenses[args->catind], (double*)0);
    std::fill (gmstruct.datatablower[args->catind], gmstruct.datatablower[args->catind]+gstruct.ges->numlenses[args->catind], (double*)0);
    std::fill (gmstruct.datatabupper[args->catind], gmstruct.datatabupper[args->catind]+gstruct.ges->numlenses[args->catind], (double*)0);
    for (int lens=args->startlens; lens<args->numlenses; lens+=args->numthreads)
      {
	GEFUNCS ge_malloc (&gmstruct.data[args->catind][lens],         dataind_num_observ);
	GEFUNCS ge_malloc (&gmstruct.datatablower[args->catind][lens], dataind_num_observ);
	GEFUNCS ge_malloc (&gmstruct.datatabupper[args->catind][lens], dataind_num_observ);
	std::copy (gstruct.ges->dataorig[args->catind][lens], gstruct.ges->dataorig[args->catind][lens]+dataind_num_observ, gmstruct.data[args->catind][lens]);
	std::copy (gstruct.ges->dataorig[args->catind][lens], gstruct.ges->dataorig[args->catind][lens]+dataind_num_observ, gmstruct.datatablower[args->catind][lens]);
	std::copy (gstruct.ges->dataorig[args->catind][lens], gstruct.ges->dataorig[args->catind][lens]+dataind_num_observ, gmstruct.datatabupper[args->catind][lens]);
      }
    
    if (tabulate_rein_bool)
      {
	gmstruct.galind = args->startlens;
	tabulate_r_ein ((void*)&gstruct);
	return 0;
      }

    for (int lens=args->startlens; lens<args->numlenses; lens+=args->numthreads)
    {
	gmstruct.galind = lens;
	
	// do calculation
	if (!gstruct.ges->fatalerror)
	  {
	    double eval = GECALCS integrate_one_lens_allerr ((void*)&gstruct);
	    args->integral += eval;
	  }
    }
    
    // cleanup
    GEFUNCS ge_free (gmstruct.xi);
    GEFUNCS ge_free (gmstruct.eta);

    for (int j=0; j<gstruct.ges->numlenses[args->catind]; ++j)
      GEFUNCS ge_free (gmstruct.data[args->catind][j]);
    GEFUNCS ge_free (gmstruct.data[args->catind]);
    GEFUNCS ge_free (gmstruct.data);

    for (int j=0; j<gstruct.ges->numlenses[args->catind]; ++j)
      GEFUNCS ge_free (gmstruct.datatablower[args->catind][j]);
    GEFUNCS ge_free (gmstruct.datatablower[args->catind]);
    GEFUNCS ge_free (gmstruct.datatablower);

    for (int j=0; j<gstruct.ges->numlenses[args->catind]; ++j)
      GEFUNCS ge_free (gmstruct.datatabupper[args->catind][j]);
    GEFUNCS ge_free (gmstruct.datatabupper[args->catind]);
    GEFUNCS ge_free (gmstruct.datatabupper);

    return 0;
}

