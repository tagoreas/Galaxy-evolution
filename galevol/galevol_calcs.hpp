#ifndef GALEVOL_CALCS_HPP_
#define GALEVOL_CALCS_HPP_

#include "galevol_defs.hpp"

#define GECALCS galevol_calcs::
class galevol_calcs
{
public:
    
    static double ein_radius           (void *parms);
    static double ein_radius_min_func  (const gsl_vector *v, void *parms_);
    static double vel_dispersion       (void *parms);
    static double integral_mdm5_norm   (double var, void *parms);
    static double integral_mstar       (double var, void *parms);
    static double integral_mdm5        (double var, void *parms);
    static double integral_mdm5_print  (double var, void *parms);
    static double integrate_one_lens   (void *parms);
    static double integrate_one_lens_allerr (void*);
    static double integrate_one_lens_2 (double*, size_t, void*);
    static double lnprobability        (void *vcomm, void *parms);
    static void*  lnprob_thread        (void *parms);
};


#endif
