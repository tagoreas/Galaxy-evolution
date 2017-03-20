#include <cstring>
#include <algorithm>
#include <glob.h>
#include <unistd.h>
#include <limits.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <iomanip>
#include "galevol_calcs.hpp"
#include "galevol_defs.hpp"
#include "galevol_funcs_templates.cpp"
#include <cmath>

#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX  1
#include <mpi.h>

int debugmode = 0;
int tabulate_rein_bool_main = 0; // tabulate Einsteiin radii or don't (only needed once)

inline std::vector<std::string> glob(const std::string& pat){
  using namespace std;
  glob_t glob_result;
  glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> ret;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    ret.push_back(string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return ret;
}

extern "C" void galevol_init1 (void *vcomm)
{   
  std::cout << "Initializing" << std::endl;

    MPI_Comm comm = *((MPI_Comm*)vcomm);
    int size, rank;
    char pname[MPI_MAX_PROCESSOR_NAME]; int len;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Get_processor_name(pname, &len);
    pname[len] = 0;

    char hostname[HOST_NAME_MAX];
    char username[LOGIN_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
    getlogin_r(username, LOGIN_NAME_MAX);

    int numcat = 1;

    // write out hostname and rank
    string fn = GEFUNCS ge_tostring(rank) + ".rank";
    std::ofstream myfile;
    myfile.open (fn.c_str());
    myfile << hostname << " " << rank << std::endl;
    myfile.close();
}

extern "C" void* galevol_init (void *vcomm)
{   
    MPI_Comm comm = *((MPI_Comm*)vcomm);
    int size, rank;
    char pname[MPI_MAX_PROCESSOR_NAME]; int len;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Get_processor_name(pname, &len);
    pname[len] = 0;

    char hostname[HOST_NAME_MAX];
    char username[LOGIN_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
    getlogin_r(username, LOGIN_NAME_MAX);

    int numcat = 1;

    // find one rank per host
    vector <string> hostnames;
    vector <int> hostranks;
    vector <string> fns = glob ("*.rank");
    int numlines;
    char **filein, *tok;
    string thishn;
    int thisrn;
    char masterrank[100];
    for (int j=0; j<(int)fns.size(); ++j)
      {
	GEFUNCS readfile (fns[j].c_str(), &filein, &numlines);
	tok = strtok (filein[0], " ");
	thishn = string (tok);
	tok = strtok (0, " ");
	thisrn = atoi (tok);
	if (std::find (hostnames.begin(), hostnames.end(), thishn) == hostnames.end() && 0!=thisrn)
	  {
	    hostnames.push_back (thishn);
	    hostranks.push_back (thisrn);
	  }
	GEFUNCS ge_free (filein, numlines);
	
	if (0==thisrn)
	  strcpy(masterrank,thishn.c_str());
      }
    
    // move masternode to end
    for (int j=0; j<(int)hostnames.size(); ++j)
      {
	if (!strcmp(masterrank,hostnames[j].c_str()))
	  {
	    hostnames.push_back(hostnames[j]);
	    hostranks.push_back(hostranks[j]);
	    hostnames.erase(hostnames.begin()+j);
	    hostranks.erase(hostranks.begin()+j);
	    break;
	  }
      }
    
    // quick conversion
    char *hostnames2;
    int  *hostranks2;
    GEFUNCS ge_malloc (&hostnames2, (hostnames.size()+1)*100);
    GEFUNCS ge_malloc (&hostranks2, hostnames.size());
    for (int j=0; j<(int)hostnames.size(); ++j)
      {
	strcpy (hostnames2+j*100,hostnames[j].c_str());
	hostranks2[j] = hostranks[j];
      }
    strcpy(hostnames2+hostnames.size()*100,masterrank); // copy name or host with rank=0

    if (1==rank)
      {
	string fn4 = "master-ranks.dat";
	std::ofstream myfile1;
	myfile1.open (fn4.c_str());
	for (int j=0; j<(int)hostranks.size(); ++j)
	  myfile1 << hostranks[j] << std::endl;
	myfile1.close();
      }

    // load Einstein radii and velocity dispersion lookup tables
    double *erlut=0, erlutlims[15], erlutints[5];
    double *vdlut=0, vdlutlims[9];
    double *ADD=0;
    for (int j=0; j<(int)hostranks.size(); ++j)
      {
	if (rank==hostranks[j])
	  {
	    std::cout << "loading Rein" << std::endl;
	    // R_Ein
	    if (!tabulate_rein_bool_main)
	      {
		FILE *fbinerlim = fopen ("erlut/limits_zd_zs_ra_ms_md.bin", "rb");
		fread ((void*)erlutlims, sizeof(erlutlims[0]), 15, fbinerlim);
		fclose (fbinerlim);
		size_t sizes[5] = {(size_t)erlutlims[10], (size_t)erlutlims[11], (size_t)erlutlims[12], (size_t)erlutlims[13], (size_t)erlutlims[14]};
		size_t totalsize = sizes[0]*sizes[1]*sizes[2]*sizes[3]*sizes[4];
		erlutints[0] = (erlutlims[5]-erlutlims[0])/(erlutlims[10]-1);
		erlutints[1] = (erlutlims[6]-erlutlims[1])/(erlutlims[11]-1);
		erlutints[2] = (erlutlims[7]-erlutlims[2])/(erlutlims[12]-1);
		erlutints[3] = (erlutlims[8]-erlutlims[3])/(erlutlims[13]-1);
		erlutints[4] = (erlutlims[9]-erlutlims[4])/(erlutlims[14]-1);
		
		GEFUNCS ge_malloc (&erlut, totalsize);
		std::fill (erlut, erlut+totalsize, 1);
		size_t currpos = 0;
		for (size_t jj=0; jj<sizes[0]; ++jj)
		  {
		    string jjstr = GEFUNCS ge_tostring (jj);
		    double *reintab;
		    size_t totalsize2 = 1;
		    for (int j=0; j<5; ++j) totalsize2 *= (size_t)sizes[j];
		    string fn5 = "erlut/erlut."+jjstr+".bin";
		    std::cout << "loading " << fn5 << std::endl;
		    FILE *fbiner = fopen (fn5.c_str(), "rb");
		    if (!debugmode)
		      fread ((void*)(erlut+currpos), sizeof(erlut[0]), totalsize2/sizes[0], fbiner);
		    fclose (fbiner);
		    currpos += totalsize2/sizes[0];
		  }
		std::cout << "done loading Rein. Read in " << currpos << " items out of " << totalsize << "." << std::endl;

		/*
		const gsl_interp2d_type *T = gsl_interp2d_bilinear;
		gsl_spline2d *spline = gsl_spline2d_alloc(T, 76,76);
		gsl_interp_accel *xacc = gsl_interp_accel_alloc();
		gsl_interp_accel *yacc = gsl_interp_accel_alloc();
		double xas[76],yas[76];
		for (int kkk=0; kkk<76;+kkk)
		  {
		    xas[kkk] = 10.5+2./75.*kkk;
		    yas[kkk] = 10.0+2./75.*kkk;
		    for (int kkk2=0; kkk2<76;+kkk2)
		      {
			gsl_spline2d_set(spline, erlut, kkk, kkk2, erlut[kkk*76+kkk2]);
		      }
		  }
		gsl_spline2d_init(spline, xas, yas, erlut, 76, 76);
		*/
	      }
	    
	    // vel disp
	    if (!tabulate_rein_bool_main)
	      {
		FILE *fbinvdlim = fopen ("vdlut/limits_rp_ms_md.bin", "rb");
		fread ((void*)vdlutlims, sizeof(vdlutlims[0]), 9, fbinvdlim);
		fclose (fbinvdlim);
		size_t sizes[3] = {(size_t)vdlutlims[6], (size_t)vdlutlims[7], (size_t)vdlutlims[8]};
		size_t totalsize = sizes[0]*sizes[1]*sizes[2];
		
		GEFUNCS ge_malloc (&vdlut, totalsize);
		string fn6 = "vdlut/vdlut.bin";
		FILE *fbinvd = fopen (fn6.c_str(), "rb");
		fread ((void*)vdlut, sizeof(vdlut[0]), totalsize, fbinvd);
		fclose (fbinvd);
	      }
	    
	    // angular diameter distance
	    {
	      GEFUNCS ge_malloc (&ADD, 100000);
	      FILE *fbinadd = fopen ("data/ADDs.bin", "rb");
	      fread ((void*)ADD, sizeof(ADD[0]), 100000, fbinadd);
	      fclose (fbinadd);
	    }
	  }
      }
    
    // load data
    string fn22 = "data/lenses.euclid.pp.dat";
    int numlines2;
    char **filein2;
    GEFUNCS readfile (fn22.c_str(), &filein2, &numlines2);
    //numlines2 = 100;//16*10*5;//800;
    int *numlenses;
    double /* ***data, */ ***dataerr, ***dataorig/*, ***datatablower, ***datatabupper */;
    GEFUNCS ge_malloc (&numlenses, numcat);
    //GEFUNCS ge_malloc (&data,      numcat, numlines2, 10);
    GEFUNCS ge_malloc (&dataerr,   numcat, numlines2, 10);
    GEFUNCS ge_malloc (&dataorig,  numcat, numlines2, 10);
    //GEFUNCS ge_malloc (&datatablower,  numcat, numlines2, 10);
    //GEFUNCS ge_malloc (&datatabupper,  numcat, numlines2, 10);
    numlenses[0] = numlines2;
    char **splitline2 = 0;
    int numsplitline2 = 0;
    int dataindex[10] = {0,2,4,6,8,10};
    int dataerrindex[10] = {1,3,5,7,9,11};
    for (int line=0; line<numlines2; ++line)
    {
	GEFUNCS ge_free (splitline2, numsplitline2);
	GEFUNCS split (filein2[line], " ", &splitline2, &numsplitline2);
	//for (int sline=0; sline<6; ++sline)
	//    data[0][line][sline]    = GEFUNCS convert_string (splitline2[dataindex[sline]]);
	for (int sline=0; sline<6; ++sline)
	    dataorig[0][line][sline]= GEFUNCS convert_string (splitline2[dataindex[sline]]);
	for (int sline=0; sline<6; ++sline)
	    dataerr[0][line][sline] = GEFUNCS convert_string (splitline2[dataerrindex[sline]]);
    }
    GEFUNCS ge_free (splitline2, numsplitline2);
    GEFUNCS ge_free (filein2, numlines2);

    // create data structures
    galevolstruct *gestruct;
    galstruct *gstruct;
    GEFUNCS ge_malloc (&gestruct, 1);
    GEFUNCS ge_malloc (&gstruct,  1);
    gstruct->ges = gestruct;
    gestruct->numlenses = numlenses;
    //gestruct->data      = data;
    gestruct->dataorig  = dataorig;
    gestruct->dataerr   = dataerr;
    //gestruct->datatablower   = datatablower;
    //gestruct->datatabupper   = datatabupper;
    gestruct->vdlut     = vdlut;
    gestruct->erlut     = erlut;
    gestruct->ADD       = ADD;
    gestruct->numhosts  = (int)hostnames.size();
    gestruct->hostnames = hostnames2;
    gestruct->hostranks = hostranks2;
    std::copy (vdlutlims,vdlutlims+9, gestruct->vdlutlims);
    std::copy (erlutlims,erlutlims+15,gestruct->erlutlims);
    std::copy (erlutints,erlutints+5, gestruct->erlutints);
    gestruct->vdlutsize[0] = (int)vdlutlims[6];
    gestruct->vdlutsize[1] = (int)vdlutlims[7];
    gestruct->vdlutsize[2] = (int)vdlutlims[8];
    gestruct->erlutsize[0] = (size_t)erlutlims[10];
    gestruct->erlutsize[1] = (size_t)erlutlims[11];
    gestruct->erlutsize[2] = (size_t)erlutlims[12];
    gestruct->erlutsize[3] = (size_t)erlutlims[13];
    gestruct->erlutsize[4] = (size_t)erlutlims[14];
    gestruct->erlutsizemul[0] = gestruct->erlutsize[1]*gestruct->erlutsize[2]*gestruct->erlutsize[3]*gestruct->erlutsize[4];
    gestruct->erlutsizemul[1] = gestruct->erlutsize[2]*gestruct->erlutsize[3]*gestruct->erlutsize[4];
    gestruct->erlutsizemul[2] = gestruct->erlutsize[3]*gestruct->erlutsize[4];
    gestruct->erlutsizemul[3] = gestruct->erlutsize[4];

    // create function pointers
    GEFUNCS ge_malloc (&gestruct->catfunc, numcat, catind_num_cat);
    //gestruct->catfunc[0][0] = (void*)(&GECALCS mu_star_slacs);
    //gestruct->catfunc[0][1] = (void*)(&GECALCS mu_r_slacs);

    // create a startup model
    GEFUNCS ge_malloc (&gestruct->psi,    numcat, psiind_num_psi);
    GEFUNCS ge_malloc (&gestruct->lambda, numcat, lambdaind_num_lambda);
    GEFUNCS ge_malloc (&gestruct->theta,  thetaind_num_theta);
    
    gestruct->psi[0][psiind_mu_star_0]  = std::pow (10.0, 11.66);
    gestruct->psi[0][psiind_zeta_star]  = 2.36;
    gestruct->psi[0][psiind_sigma_star] = 0.23;
    gestruct->psi[0][psiind_mu_r_0]     = std::pow (10.0, 0.70);
    gestruct->psi[0][psiind_zeta_r]     = 0.07;
    gestruct->psi[0][psiind_beta_r]     = 0.64;
    gestruct->psi[0][psiind_sigma_r]    = 0.07;

    gestruct->theta[thetaind_gamma_0]     = 0.80;
    gestruct->theta[thetaind_sigma_gamma] = 0.34;
    gestruct->theta[thetaind_zeta_dm]     = 0.86;
    gestruct->theta[thetaind_beta_dm]     = 0.05;
    gestruct->theta[thetaind_xi_dm]       = -0.49;
    gestruct->theta[thetaind_m_dm_0]      = std::pow (10.0, 10.69);
    gestruct->theta[thetaind_sigma_dm]    = 0.26;
    gestruct->theta[thetaind_zeta_imf]    = -0.11;
    gestruct->theta[thetaind_beta_imf]    = 0.19;
    gestruct->theta[thetaind_xi_imf]      = 0.09;
    gestruct->theta[thetaind_alpha_imf_0] = std::pow (10.0, 0.04);
    gestruct->theta[thetaind_sigma_imf]   = 0.02;

    gestruct->lambda[0][lambdaind_r_sel]     = 0.96;
    gestruct->lambda[0][lambdaind_sigma_sel] = 0.30;

    // save setup
    return (void*)gstruct;
}

extern "C" void galevol_killme ()
{
  exit(0);
}

extern "C" double galevol_sample (void *vcomm, void *gstruct_, double *geparms)
{
    //double *geparms = (double*)geparms_;
    galstruct *gstruct = (galstruct*)gstruct_;
    galevolstruct *gestruct = gstruct->ges;
    //std::cout << "grepme "; for (int i=0; i<19; ++i) std::cout << geparms[i] << " ";
    int ind=0;
    gestruct->psi[0][psiind_mu_star_0]       = geparms[ind++]; 
    gestruct->psi[0][psiind_zeta_star]       = geparms[ind++]; 
    gestruct->psi[0][psiind_sigma_star]      = geparms[ind++]; 
    gestruct->psi[0][psiind_mu_r_0]          = geparms[ind++]; 
    gestruct->psi[0][psiind_zeta_r]          = geparms[ind++]; 
    gestruct->psi[0][psiind_beta_r]          = geparms[ind++]; 
    gestruct->psi[0][psiind_sigma_r]         = geparms[ind++]; 
    //gestruct->theta[thetaind_gamma_0]        = geparms[ind++]; 
    //gestruct->theta[thetaind_sigma_gamma]    = geparms[ind++]; 
    gestruct->theta[thetaind_zeta_dm]        = geparms[ind++]; 
    gestruct->theta[thetaind_beta_dm]        = geparms[ind++]; 
    gestruct->theta[thetaind_xi_dm]          = geparms[ind++]; 
    gestruct->theta[thetaind_m_dm_0]         = geparms[ind++]; 
    gestruct->theta[thetaind_sigma_dm]       = geparms[ind++]; 
    gestruct->theta[thetaind_zeta_imf]       = geparms[ind++]; 
    gestruct->theta[thetaind_beta_imf]       = geparms[ind++]; 
    gestruct->theta[thetaind_xi_imf]         = geparms[ind++]; 
    gestruct->theta[thetaind_alpha_imf_0]    = geparms[ind++]; 
    gestruct->theta[thetaind_sigma_imf]      = geparms[ind++]; 
    gestruct->lambda[0][lambdaind_r_sel]     = geparms[ind++]; 
    gestruct->lambda[0][lambdaind_sigma_sel] = geparms[ind++]; 

    // set bounds
    ind=0;
    double **bounds;
    GEFUNCS ge_malloc (&bounds, 21-2, 2);
    bounds[ind][0] = GEDEFS mstarlim[0]; bounds[ind++][1] = GEDEFS mstarlim[1];
    bounds[ind][0] = -1;                 bounds[ind++][1] = 10;
    bounds[ind][0] = -1.5;               bounds[ind++][1] = 0;
    bounds[ind][0] = -4;                 bounds[ind++][1] = 1;
    bounds[ind][0] = -6.5;               bounds[ind++][1] = 5;
    bounds[ind][0] = -10;                bounds[ind++][1] = 4.5;
    bounds[ind][0] = -1.5;               bounds[ind++][1] = 0;
    //bounds[ind][0] = -1e300;             bounds[ind++][1] = 1e300;
    //bounds[ind][0] = -1e300;             bounds[ind++][1] = 1e300;
    bounds[ind][0] = -2.5;               bounds[ind++][1] = 10;
    bounds[ind][0] = -6;                bounds[ind++][1] = 3.5;
    bounds[ind][0] = -2;                bounds[ind++][1] = 2;
    bounds[ind][0] = GEDEFS mdm5lim[0];  bounds[ind++][1] = GEDEFS mdm5lim[1];
    bounds[ind][0] = -1.5;               bounds[ind++][1] = 0.5;
    bounds[ind][0] = -2.5;                bounds[ind++][1] = 2;
    bounds[ind][0] = -2;                bounds[ind++][1] = 4;
    bounds[ind][0] = -1;                bounds[ind++][1] = 5;
    bounds[ind][0] = -1;                 bounds[ind++][1] = 1;
    bounds[ind][0] = -2;                 bounds[ind++][1] = 0;
    bounds[ind][0] = -1;                 bounds[ind++][1] = 1;
    bounds[ind][0] = -1;                 bounds[ind++][1] = 1;
    
    ind=0;
    
    bounds[ind][0] = GEDEFS mstarlim[0];    bounds[ind++][1] = GEDEFS mstarlim[1];
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    //bounds[ind][0] = -1e300;             bounds[ind++][1] = 1e300;
    //bounds[ind][0] = -1e300;             bounds[ind++][1] = 1e300;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = GEDEFS mdm5lim[0];    bounds[ind++][1] = GEDEFS mdm5lim[1];
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -100;    bounds[ind++][1] = 2;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    bounds[ind][0] = -10;    bounds[ind++][1] = 10;
    

    
    
    double penaltydist = 0;
    for (int m=0; m<21-2; ++m)
    {
	double diff[2] = {bounds[m][1]-geparms[m], geparms[m]-bounds[m][0]};
	for (int n=0; n<2; ++n)
	    if (diff[n]<0)
	      {
		//std::cout << "bounds" << m << ": " << bounds[m][0] << " " << geparms[m] << " " << bounds[m][1] << std::endl;
		penaltydist += diff[n]*diff[n];
	    }
    }
    GEFUNCS ge_free (bounds, 21-2);
    if (penaltydist)
      {
	//return -std::numeric_limits<double>::infinity();
	return GEDEFS huge_penalty *(1.0+penaltydist);
	//	usleep((unsigned long long int)(10000000000));
      }

    //gestruct->psi[0][psiind_mu_r_0]       = std::pow (10.0, gestruct->psi[0][psiind_mu_r_0]);
    //gestruct->psi[0][psiind_mu_star_0]    = std::pow (10.0, gestruct->psi[0][psiind_mu_star_0]);
    //gestruct->theta[thetaind_m_dm_0]      = std::pow (10.0, gestruct->theta[thetaind_m_dm_0]);
    //gestruct->theta[thetaind_alpha_imf_0] = std::pow (10.0, gestruct->theta[thetaind_alpha_imf_0]);
    gestruct->psi[0][psiind_sigma_star]   = std::pow (10.0, gestruct->psi[0][psiind_sigma_star]);
    gestruct->psi[0][psiind_sigma_r]      = std::pow (10.0, gestruct->psi[0][psiind_sigma_r]);
    gestruct->theta[thetaind_sigma_dm]    = std::pow (10.0, gestruct->theta[thetaind_sigma_dm]);
    gestruct->theta[thetaind_sigma_imf]   = std::pow (10.0, gestruct->theta[thetaind_sigma_imf]);
    gestruct->lambda[0][lambdaind_r_sel]  = std::pow (10.0, gestruct->lambda[0][lambdaind_r_sel]);
    gestruct->lambda[0][lambdaind_sigma_sel] = std::pow (10.0, gestruct->lambda[0][lambdaind_sigma_sel]);
    //gestruct->psi[0][psiind_mu_r_0]          = std::log10(gestruct->psi[0][psiind_mu_r_0]);
    

    
    // distribute memory to other nodes in cluster
    MPI_Comm comm = *((MPI_Comm*)vcomm);
    int size, rank;
    char pname[MPI_MAX_PROCESSOR_NAME]; int len;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Get_processor_name(pname, &len);
    pname[len] = 0;
    

    double lnprob = GECALCS lnprobability (vcomm, (void*)gstruct);
    if (!lnprob==lnprob)
      lnprob = GEDEFS huge_penalty;
    //lnprob += -0.5 *std::pow(gestruct->lambda[0][lambdaind_r_sel]-1,2.0) /(0.2*0.2);
    //lnprob += -0.5 *std::pow(gestruct->lambda[0][lambdaind_sigma_sel]-1,2.0) /(0.2*0.2);
    //std::cout << "ln prob = " << lnprob << std::endl;

    //std::cout << "Hello, World! I am process " << rank << " of " 
    //<< size << " on " << pname << ": ans = " << lnprob << std::endl;

    return lnprob;
}
