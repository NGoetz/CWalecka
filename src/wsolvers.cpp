#include <wfunctions.h>
#include <wsolvers.h>
#include <wplot.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <matplot/matplot.h>
#include <vector>
#include <cassert>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>

double round_to_n_digits(double x, int n)
{ 
    double scale = pow(10.0, fabs(ceil(log10(fabs(x)))) + n);

    return round(x * scale) / scale;
}

//single step of the root solving for the self-consistent equation of m* and scalar density at T=0; x is the vector of values to be determined,
// params the parameters and f the value of the function to be reduced to zero
int mass_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f){
        vector<double> scalar_coeff = ((struct mass_eff_scalar_density_params *) params)->scalar_coeff;
        vector<double> scalar_exp = ((struct mass_eff_scalar_density_params *) params)->scalar_exp;
        double number_density=((struct mass_eff_scalar_density_params *) params)->number_density;
        double nucleon_mass=((struct mass_eff_scalar_density_params *) params)->nucleon_mass;
        double degeneracy=(( struct mass_eff_scalar_density_params *) params )->degeneracy;
        //x0 is the prediction of the effective mass
        const double x0 = gsl_vector_get (x, 0);
        // the calculated minus the predicted value for the mass has to be zero
        const double y0 =x0-mass_eff_from_scalar_density(nucleon_mass, scalar_coeff, scalar_exp, scalar_density(degeneracy,x0, number_density));
        gsl_vector_set (f, 0, y0);
        return GSL_SUCCESS;
}

//printout of the status for debug
int print_state_mass_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s){
    printf ("iter = %3lu x = % .3f  "
            "f(x) = % .3e \n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->f, 0));
    return 0;
}
// solve the self consistent equation. This takes the scalar coupling constant
// as a parameter!
pair<double, double> get_mass_eff_scalar_density(double degeneracy, double nucleon_mass,vector<double> scalar_coeff ,vector<double> scalar_exp,double number_density,bool print) {

    struct mass_eff_scalar_density_params p ={scalar_coeff,scalar_exp,number_density,nucleon_mass, degeneracy,{0},{0}};
    
    int status;
    size_t iter = 0;

    const size_t num = 1;
    
    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (T, num);
    
    
    gsl_multiroot_function f = {&mass_eff_scalar_density_root, num, &p};
    // this is the initialization of the effective mass
    // it should maybe be handeled centreally
    double x_init[1] = {0.5*nucleon_mass};
    gsl_vector *x = gsl_vector_alloc (num);

    gsl_vector_set (x, 0, x_init[0]);
    
    gsl_multiroot_fsolver_set (s, &f, x);

    if (print) {
        print_state_mass_eff_scalar_density (iter, s);
    }
    // perform the root solving
    do
        {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (print) {
            print_state_mass_eff_scalar_density (iter, s);
        }

        if (status)   /* check if solver is stuck */
            break;

        status =
            gsl_multiroot_test_residual (s->f, 1e-7);
        }
    while (status == GSL_CONTINUE && iter < 1000);
    if (print){
        printf ("status = %s\n", gsl_strerror (status));
    }
    // we extract the values and calculate ns from this
    // both is returned in a pair
    double mass_eff=gsl_vector_get(s->x,0);
    double scalar_density_val=scalar_density(degeneracy, mass_eff, number_density);
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    pair<double, double> result={mass_eff,scalar_density_val};
    return result;
}

//single step of the root solving for the self-consistent equation at finite T; x is the vector of values to be determined,
// params the parameters and f the value of the function to be reduced to zero
int T_mass_eff_mu_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f){
        vector<double> scalar_coeff  = ((struct mass_eff_scalar_density_params *) params)->scalar_coeff;
       vector<double> scalar_exp  = ((struct mass_eff_scalar_density_params *) params)->scalar_exp;
        double nucleon_mass=((struct mass_eff_scalar_density_params *) params)->nucleon_mass;
        double number_density=((struct mass_eff_scalar_density_params *) params)->number_density;
        double temperature=((struct mass_eff_scalar_density_params *) params)->temperature;
        double degeneracy=(( struct mass_eff_scalar_density_params *) params )->degeneracy;
        const double x0 = abs(gsl_vector_get (x, 0));
        const double x1 = abs(gsl_vector_get (x, 1));//prediction of effective baryochemical potential
        const double y0 =x0-mass_eff_from_scalar_density(nucleon_mass,scalar_coeff,scalar_exp,T_scalar_density(degeneracy,temperature,x0,x1));
        const double y1 =number_density-T_number_density(degeneracy,temperature,x0,x1);
        gsl_vector_set (f, 0, y0);
        gsl_vector_set (f, 1, y1);
        return GSL_SUCCESS;
}

//printout of the status for debug
int print_state_T_mass_eff_mu_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s){
    printf ("iter = %3lu x = % .3f % .3f "
            "f(x) = % .3e % .3e \n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1));
    return 0;
}
// solve the self consistent equation for m*,mu* and ns at finite T
tuple<double, double,double,bool> get_T_mass_eff_mu_eff_scalar_density(vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double nucleon_mass,double number_density, double temperature, double degeneracy, bool print) {
    struct mass_eff_scalar_density_params p ={scalar_coeff,scalar_exp,number_density,nucleon_mass,degeneracy,  temperature, vec_coeff,vec_exp};
    int status;
    size_t iter = 0;
    const size_t num = 2;
    gsl_vector *x = gsl_vector_alloc (num);
    const gsl_multiroot_fsolver_type * K = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (K, num);
  
    gsl_multiroot_function f = {&T_mass_eff_mu_eff_scalar_density_root, num, &p};

    double x_init[2] = {nucleon_mass,nucleon_mass};
    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    
    gsl_multiroot_fsolver_set (s, &f, x);

    if (print) {
        print_state_T_mass_eff_mu_eff_scalar_density (iter, s);
    }
    // perform the root solving
    do{
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (print) {
            print_state_T_mass_eff_mu_eff_scalar_density (iter, s);
        }

        if (status)   /* check if solver is stuck */
            break;

        status =
            gsl_multiroot_test_residual (s->f, 1e-7);
        }while (status == GSL_CONTINUE && iter < 1000);
    if (print){
        printf ("status = %s\n", gsl_strerror (status));
    
    }

    // we extract the values and calculate ns from this
    // the result is returned in a tuple
    double mass_eff=abs(gsl_vector_get(s->x,0));
    double mu_eff=abs(gsl_vector_get(s->x,1));
    double scalar_density=T_scalar_density(degeneracy,temperature, mass_eff, mu_eff);
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    tuple<double, double,double,bool> result={mass_eff,mu_eff,scalar_density,status==GSL_SUCCESS};
    //cout<<"finished "<<"\n";
    return result;
}

// the single root solving step for the coupling constants at T=0
// this only works for two terms due to lack of conditions!
int coeff_root(const gsl_vector * x, void *params, gsl_vector * f){
    double saturation_density=((struct coeff_params *) params)->saturation_densitiy;
    double nucleon_mass=((struct coeff_params *) params)->nucleon_mass;
    double binding_energy=((struct coeff_params *) params)->binding_energy;
    double degeneracy=((struct coeff_params *) params)->degeneracy;
    vector<double> scalar_exp=((struct coeff_params *) params)->scalar_exp;
    vector<double> vec_exp=((struct coeff_params *) params)->vec_exp;
    unsigned int* terms =((struct coeff_params *) params)->terms;
    
    vector<double>xs,xv;
    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    switch(terms[0]){
        case 0:
            xs={};
            xv={x0,x1};
            break;
        case 1:
            xs={x0};
            xv={x1};
            break;
        case 2:
            xs={x0,x1};
            xv={};
            break;
    }
    // we have to first calculate the effective mass with the
    // guess we have for the coefficients
    pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,xs, scalar_exp, saturation_density, false);
    // here,  the minimum condition for stable nuclei is given
    const double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, xs,  scalar_exp, xv,  vec_exp, output.second);
    const double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density, xs,  scalar_exp, xv,  vec_exp, output.second);
    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);
    return GSL_SUCCESS;
}
// same as before
int print_state_coeff (size_t iter, gsl_multiroot_fsolver * s){
    printf ("iter = %3lu x = % .3f % .3f "
            "f(x) = % .3e % .3e \n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1));
    return 0;
}
// same as before
pair<double, double> get_coeff(double nucleon_mass, double saturation_density, double binding_energy, double degeneracy, vector<double> scalar_exp, vector<double> vec_exp, bool print, unsigned int terms[2], bool * success){
    assert(terms[0]+terms[1]==2);
    struct coeff_params p ={saturation_density,nucleon_mass, binding_energy,degeneracy,scalar_exp, vec_exp, terms};
    int status;
    size_t iter = 0;

    const size_t num = 2;
    
    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (T, num);
    
    
    gsl_multiroot_function f = {&coeff_root, num, &p};

    double x_init[2] = {10*4*M_PI/pow(550,2.0),10*4*M_PI/pow(750,2.0)};
    gsl_vector *x = gsl_vector_alloc (num);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    
    gsl_multiroot_fsolver_set (s, &f, x);
    if (print){
        print_state_coeff (iter, s);
    }

    do
        {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if(print){
            print_state_coeff (iter, s);
        }

        if (status)   /* check if solver is stuck */
            break;

        status =
            gsl_multiroot_test_residual (s->f, 1e-13);
        }
    while (status == GSL_CONTINUE && iter < 300);
    if (print){
        printf ("status = %s\n", gsl_strerror (status));
    }
    double a=gsl_vector_get(s->x,0);
    double b=gsl_vector_get(s->x,1);
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    pair<double, double> result={a,b};
    if(status==GSL_SUCCESS&&success!=NULL){
        *success=true;
    }
    return result;
}

int crit_root(const gsl_vector * x, void *params, gsl_vector * f){
    vector<double> scalar_coeff=((struct crit_params *) params)->scalar_coeff;
    vector<double> scalar_exp=((struct crit_params *) params)->scalar_exp;
    vector<double> vec_coeff=((struct crit_params *) params)->vec_coeff;
    vector<double> vec_exp=((struct crit_params *) params)->vec_exp;
    double nucleon_mass=((struct crit_params *) params)->nucleon_mass;
    double degeneracy=((struct crit_params *) params)->degeneracy;
    
    
    const double x0=abs(gsl_vector_get(x,0));//number density
    const double x1=abs(gsl_vector_get(x,1));//temperature - both need to be positive
    
    // find the inflection point
    const double y0 =T_press_dn_solv(degeneracy, nucleon_mass, x0, scalar_coeff, scalar_exp, vec_coeff, vec_exp, x1);
    const double y1 =T_press_dn2_solv(degeneracy, nucleon_mass, x0,scalar_coeff, scalar_exp, vec_coeff, vec_exp, x1);
    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);
    return GSL_SUCCESS;
}
//find the critical point
pair<double, double> get_crit(double nucleon_mass,double saturation_density, double binding_energy, double degeneracy, vector<double> scalar_exp, vector<double> vec_exp, bool print, double crit_T, double crit_n, vector<double> scalar_coeff, vector<double> vec_coeff){
    struct crit_params p ;
    //if no coefficients are given, the model was only determined from the saturation point
    if(abs(scalar_coeff[0])<1e-10&&abs(vec_coeff[0])<1e-10){
        unsigned int terms [2]={scalar_exp.size(),vec_exp.size()};
        pair<double, double> coeff=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, scalar_exp, vec_exp, false,terms);
        p ={{coeff.first}, scalar_exp, {coeff.second}, vec_exp, nucleon_mass, degeneracy};
    }else{
         p ={scalar_coeff, scalar_exp, vec_coeff, vec_exp, nucleon_mass, degeneracy};
    }
    int status;
    size_t iter = 0;

    const size_t num = 2;
    
    const gsl_multiroot_fsolver_type * Q = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (Q, num);
    
    
    gsl_multiroot_function f = {&crit_root, num, &p};

    double x_init[num] = {crit_n,crit_T};
    gsl_vector *x = gsl_vector_alloc (num);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    
    gsl_multiroot_fsolver_set (s, &f, x);
    if (print){
        print_state_coeff (iter, s);
    }

    do
        {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if(print){
            print_state_coeff (iter, s);
        }

        if (status)   /* check if solver is stuck */
            break;

        status =
            gsl_multiroot_test_residual (s->f, 1e-13);
        }
    while (status == GSL_CONTINUE && iter < 1000);
    if (print){
        printf ("status = %s\n", gsl_strerror (status));
    }
    double temperature_crit=abs(gsl_vector_get(s->x,1));
    double number_density_crit=abs(gsl_vector_get(s->x,0));
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    pair<double, double> result={temperature_crit,number_density_crit};
    return result;
}
//currently only for 2 terms
int interaction_root(const gsl_vector * x, void *params, gsl_vector * f){
    double nucleon_mass=((struct interaction_params *) params)->nucleon_mass;
    double degeneracy=((struct interaction_params *) params)->degeneracy;
    double saturation_density=((struct interaction_params *) params)->saturation_density;
    double binding_energy=((struct interaction_params *) params)->binding_energy;
    double critical_temperature=((struct interaction_params *) params)->critical_temperature;
    double critical_density=((struct interaction_params *) params)->critical_density;
    unsigned int * terms =((struct interaction_params *) params)->terms;
    
    vector<double> scalar_coeff, scalar_exp, vector_coeff, vector_exp;
  
    const double x0=gsl_vector_get(x,0);//coefficient 1
    const double x1=(gsl_vector_get(x,1));//coefficient 2
    const double x2=gsl_vector_get(x,2);//exponent 1
    const double x3=gsl_vector_get(x,3);//exponent 2
    switch(terms[0]){
        case 0:
            scalar_coeff={};
            scalar_exp={};
            vector_coeff={x0,x1};
            vector_exp={x2,x3};
            break;
        case 1:
            scalar_coeff={x0};
            scalar_exp={x2};
            vector_coeff={x1};
            vector_exp={x3};
            break;
        case 2:
            scalar_coeff={x0,x1};
            scalar_exp={};
            vector_coeff={};
            vector_exp={x2,x3};
            break;
    }
    pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,scalar_coeff, scalar_exp, saturation_density, false);
    //predict binding at demanded saturation density
    const double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, scalar_coeff, scalar_exp,vector_coeff, vector_exp, output.second);
    const double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,scalar_coeff, scalar_exp,vector_coeff, vector_exp, output.second);
    //there is an inflection point at the critical temperature/density
    pair<double, double> pressures=T_press_dn12_solv( degeneracy,  nucleon_mass, critical_density,  scalar_coeff, scalar_exp,vector_coeff, vector_exp, critical_temperature);
    const double y2 =pressures.first;
    const double y3 =pressures.second;
    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);
    gsl_vector_set (f, 2, y2);
    gsl_vector_set (f, 3, y3);
    return GSL_SUCCESS;
}

int print_state_interaction (size_t iter, gsl_multiroot_fsolver * s){
    printf ("iter = %3lu x = % .3f % .3f % .3f % .3f "
            "f(x) = % .3e % .3e % .3e % .3e\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1),
            gsl_vector_get (s->f, 2),
            gsl_vector_get (s->f, 3));
    return 0;
}
//find the coefficients and exponents of 2 terms with the nuclear saturation point and the nuclear liquid/gas transition
tuple<double, double,double,double, double, bool> get_interaction_4D(void * p, bool print, vector<double> init_scalar_exp, vector<double> init_vec_exp, vector<double> init_scalar_coeff, vector<double> init_vec_coeff) {
    gsl_set_error_handler(&gsl_handler);
    int status;
    size_t iter = 0;
    const size_t num = 4;
    gsl_vector *x = gsl_vector_alloc (num);
    const gsl_multiroot_fsolver_type * K = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (K, num);
  
    gsl_multiroot_function f = {&interaction_root, num, p};
    unsigned int * terms=((struct interaction_params *) p)->terms;
    assert(terms[0]+terms[1]==2);
    
    double x_init[4];
    for(unsigned int i=0; i<terms[0];i++){
        x_init[i]=init_scalar_coeff[i];
        x_init[2+i]=init_scalar_exp[i];
    }
    for(unsigned int i=terms[0]; i<terms[0]+terms[1];i++){
        x_init[i]=init_vec_coeff[i-terms[0]];
        x_init[2+i]=init_vec_exp[i-terms[0]];
    }
    
    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    gsl_vector_set (x, 2, x_init[2]);
    gsl_vector_set (x, 3, x_init[3]);
    
    gsl_multiroot_fsolver_set (s, &f, x);

    if (print) {
        print_state_interaction (iter, s);
    }
    // perform the root solving
    do{
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (print) {
            print_state_interaction (iter, s);
        }
        
        if (status)   /* check if solver is stuck */
            break;

        status =
            gsl_multiroot_test_residual (s->f, 1e-12);
        
        }while (status == GSL_CONTINUE && iter < 1000);
    if (print){
        printf ("status = %s\n", gsl_strerror (status));
    
    }
    double error=abs(gsl_vector_get(s->f,2))+abs(gsl_vector_get(s->f,3));
    tuple<double, double,double,double, double, bool> result={gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),gsl_vector_get(s->x,2),gsl_vector_get(s->x,3),error,status==GSL_SUCCESS};
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    return result;
}
// solve the energy density at T=0 with given coupling constants
double energy_pp_minus_mass_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, scalar_coeff, scalar_exp, number_density, false);
    return eps_over_n(degeneracy, x.first,number_density,scalar_coeff, scalar_exp, vec_coeff, vec_exp, x.second)-nucleon_mass;
}

// derivative of this
double energy_pp_minus_mass_dn_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, scalar_coeff, scalar_exp, number_density, false);
    return eps_over_n_dn(degeneracy, x.first,number_density,scalar_coeff, scalar_exp, vec_coeff, vec_exp, x.second);
}
//same for pressure as T=0
double press_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, scalar_coeff, scalar_exp, number_density, false);
    return press(degeneracy, x.first,number_density,scalar_coeff, scalar_exp, vec_coeff, vec_exp,x.second);
}
//incompressibility at T=0
double incsolv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double delta){
    pair<double, double> x1=get_mass_eff_scalar_density(degeneracy, nucleon_mass, scalar_coeff, scalar_exp, number_density+delta, false);
    pair<double, double> x2=get_mass_eff_scalar_density(degeneracy, nucleon_mass, scalar_coeff, scalar_exp, number_density-delta, false);
    return 9*number_density*number_density*(eps_over_n_dn(degeneracy, x1.first,number_density+delta,scalar_coeff, scalar_exp, vec_coeff, vec_exp, x1.second)-eps_over_n_dn(degeneracy, x2.first,number_density-delta,scalar_coeff, scalar_exp, vec_coeff, vec_exp, x2.second))/(2*delta);
}
//pressure at finite T
double T_press_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density, temperature, degeneracy,   false);
    return T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x), number_density);
}
//derivative of pressure at finite T after mu
double  T_press_dmu_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp, nucleon_mass,number_density, temperature, degeneracy,   false);
    return (T_pressure(degeneracy, temperature,get<0>(x),get<1>(x)+delta,scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x), number_density)-T_pressure(degeneracy, temperature,get<0>(x),get<1>(x)-delta,scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x), number_density))/2*delta;
}
//derivative of pressure at finite T after T
double  T_press_dT_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density, temperature+delta, degeneracy,   false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density, temperature-delta, degeneracy,   false);
    return (T_pressure(degeneracy, temperature+delta,get<0>(x1),get<1>(x1),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x1), number_density)-T_pressure(degeneracy, temperature-delta,get<0>(x2),get<1>(x2),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x2), number_density))/2*delta;
}
//derivative of pressure at finite T after the number density
double  T_press_dn_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density+delta, temperature, degeneracy,   false);
    
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density-delta, temperature, degeneracy,   false);
    
    return (T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x1), number_density+delta)-T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x2), number_density-delta))/(2*delta);
}
//second derivative of pressure at finite T after the number density
double  T_press_dn2_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density+delta, temperature, degeneracy,   false);
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density, temperature, degeneracy,   false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density-delta, temperature, degeneracy,   false);
    return (T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x1), number_density+delta)+T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x2), number_density-delta)-2*T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x), number_density))/(delta*delta);
}
// energy density at finite T
double T_eps_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density, temperature, degeneracy,   false);
    double vec=0;
    for(unsigned int i=0;i<vec_coeff.size(); i++){
        vec+=vec_coeff[i]*pow(number_density,vec_exp[i]-1);
    }
    double mu=get<1>(x)+vec;//get mu from n
    return -T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x), number_density)+mu*number_density+temperature*T_press_dT_solv(degeneracy, nucleon_mass, number_density,  scalar_coeff, scalar_exp, vec_coeff, vec_exp, temperature);
}
//get both the first and second derivate of the pressure at the same time
pair<double, double> T_press_dn12_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density+delta, temperature, degeneracy,   false);
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density, temperature, degeneracy,   false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass,number_density-delta, temperature, degeneracy,   false);
    double pressure1=T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x1), number_density+delta);
    double pressure2=T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x2), number_density-delta);
    double pressure=T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),scalar_coeff, scalar_exp, vec_coeff, vec_exp,get<2>(x), number_density);
    pair<double,double> res={(-pressure2+pressure1)/(2*delta),(pressure1+pressure2-2*pressure)/(delta*delta)};
    return res;
}
//gsl handler for dealing with negative mass events
void gsl_handler (const char * reason,
              const char * file,
              int line,
              int gsl_errno){
    cout<<reason<<" in "<<file<<" at line "<<line<<" with code "<<gsl_errno<<"\n";
    throw 42;
            }

//validate interaction terms found
bool validate_interaction(double error, vector<double> coeff_guess, vector<double> exp_guess,vector<double> coeff, vector<double> exp, string filename, int timestamp, void* params){
    double eps=numeric_limits<double>::epsilon();
    double nucleon_mass=((struct interaction_params *) params)->nucleon_mass;
    double degeneracy=((struct interaction_params *) params)->degeneracy;
    double saturation_density=((struct interaction_params *) params)->saturation_density;
    double binding_energy=((struct interaction_params *) params)->binding_energy;
    double critical_temperature=((struct interaction_params *) params)->critical_temperature;
    double critical_density=((struct interaction_params *) params)->critical_density;
    unsigned int * terms =((struct interaction_params *) params)->terms;
    
    vector<double> scalar_coeff, scalar_exp, vec_coeff, vec_exp,scalar_coeff_guess, scalar_exp_guess, vec_coeff_guess, vec_exp_guess;
    
    for(unsigned int i=0; i<terms[0]; i++){
       scalar_coeff.push_back(coeff[i]); 
       scalar_exp.push_back(exp[i]); 
       scalar_coeff_guess.push_back(coeff_guess[i]); 
       scalar_exp_guess.push_back(exp_guess[i]); 
    }
    for(unsigned int i=terms[0]; i<terms[0]+terms[1]; i++){
       vec_coeff.push_back(coeff[i]); 
       vec_exp.push_back(exp[i]); 
       vec_coeff_guess.push_back(coeff_guess[i]); 
       vec_exp_guess.push_back(exp_guess[i]); 
    }
   
    pair<double, double> mstarns=get_mass_eff_scalar_density( degeneracy,  nucleon_mass, scalar_coeff,  scalar_exp,  saturation_density, false);
    
    double new_binding_energy=nucleon_mass-eps_over_n(degeneracy, mstarns.first, saturation_density, scalar_coeff,  scalar_exp, vec_coeff, vec_exp, mstarns.second);
    double binding_energy_minimum =eps_over_n_dn(degeneracy, mstarns.first, saturation_density, scalar_coeff,  scalar_exp, vec_coeff, vec_exp, mstarns.second);
    double incompress=incsolv( degeneracy,  nucleon_mass,saturation_density, scalar_coeff,  scalar_exp, vec_coeff, vec_exp);
    pair<double, double> CP= get_crit( nucleon_mass, saturation_density,  binding_energy,  degeneracy,  scalar_exp,  vec_exp, false, 19,0.44*saturation_density,scalar_coeff,vec_coeff);
    
    double CP_P=T_press_solv(degeneracy, nucleon_mass,CP.second,  scalar_coeff, scalar_exp, vec_coeff, vec_exp, CP.first)/pow(197.3,3);
    
    if(abs(round_to_n_digits(CP.first,5)-round_to_n_digits(critical_temperature,5))>eps||abs(round_to_n_digits( CP.second/saturation_density,3)-round_to_n_digits( critical_density/saturation_density,3))>eps){
        return false;
    }
    FILE * pFile;
    pFile = fopen ((filename+".csv").c_str(),"a");
    fprintf(pFile,"%d, ",timestamp);
    for(unsigned int i=0; i<scalar_coeff_guess.size();i++){
        fprintf(pFile, "%.16e, ",scalar_coeff_guess[i]);
    }
    for(unsigned int i=0; i<scalar_exp_guess.size();i++){
        fprintf(pFile, "%.16f, ",scalar_exp_guess[i]);
    }
    for(unsigned int i=0; i<vec_coeff_guess.size();i++){
        fprintf(pFile, "%.16e, ",vec_coeff_guess[i]);
    }
    for(unsigned int i=0; i<vec_exp_guess.size();i++){
        fprintf(pFile, "%.16f, ",vec_exp_guess[i]);
    }
    for(unsigned int i=0; i<scalar_coeff.size();i++){
        fprintf(pFile, "%.16e, ",scalar_coeff[i]);
    }
    for(unsigned int i=0; i<scalar_exp.size();i++){
        fprintf(pFile, "%.16f, ",scalar_exp[i]);
    }
    for(unsigned int i=0; i<vec_coeff.size();i++){
        fprintf(pFile, "%.16e, ",vec_coeff[i]);
    }
    for(unsigned int i=0; i<vec_exp.size();i++){
        fprintf(pFile, "%.16f, ",vec_exp[i]);
    }
    fprintf(pFile, "%.5e, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16F; \n ", error, mstarns.first, new_binding_energy, binding_energy_minimum, incompress, CP.first, CP.second/saturation_density, CP_P);
   
    fclose(pFile);
    return true;
    
}

//find interactions on a grid of initial values
void interaction_4D_grid(double nucleon_mass,double critical_temperature, double critical_density, double binding_energy, double saturation_density, double degeneracy, double boundaries  [4][3], unsigned int terms [2], string filename, bool print, int num_sol){
    gsl_set_error_handler(&gsl_handler);
    double eps=numeric_limits<double>::epsilon();
    int timestamp=time(0);
    if(print){
        FILE * pFile;
        pFile = fopen ((filename+".csv").c_str(),"w");
        fprintf(pFile,"timestamp, ");
        for(unsigned int i=0; i<terms[0];i++){
            fprintf(pFile, "scalar_coeff_guess[%d], ",i);
        }
        for(unsigned int i=0; i<terms[0];i++){
            fprintf(pFile, "scalar_exp_guess[%d], ",i);
        }
        for(unsigned int i=0; i<terms[1];i++){
            fprintf(pFile, "vector_coeff_guess[%d], ", i);
        }
        for(unsigned int i=0; i<terms[1];i++){
            fprintf(pFile, "vector_exp_guess[%d], ",i);
        }
        for(unsigned int i=0; i<terms[0];i++){
            fprintf(pFile, "scalar_coeff[%d], ",i);
        }
        for(unsigned int i=0; i<terms[0];i++){
            fprintf(pFile, "scalar_exp[%d], ",i);
        }
        for(unsigned int i=0; i<terms[1];i++){
            fprintf(pFile, "vector_coeff[%d], ", i);
        }
        for(unsigned int i=0; i<terms[1];i++){
            fprintf(pFile, "vector_exp[%d], ",i);
        }
            
        fprintf(pFile,"GSL_error, mstar, binding_energy, binding_energy_minimum, incompress, CP_T, CP_n, CP_P; \n");
        fclose(pFile); 
    }
    
    
    struct interaction_params p={nucleon_mass, degeneracy, saturation_density,binding_energy, critical_temperature, critical_density,  terms };
    
    unsigned int count=0;
    vector<vector<double>> solution_storage;
    int num_coeff1_guess=ceil((boundaries[0][1]-boundaries[0][0])/boundaries[0][2]);
    int num_threads=ceil(omp_get_max_threads()-1);
    
#pragma omp parallel num_threads(num_threads) shared(p, solution_storage, num_sol, count)
{
#pragma omp for
    for(int i=0; i<=num_coeff1_guess; i+=1 ){
        double coeff1_guess=boundaries[0][0]+i*boundaries[0][2];
        
        for(double coeff2_guess=boundaries[1][0]; coeff2_guess<=boundaries[1][1]; coeff2_guess+=boundaries[1][2] ){

            for(double exp1_guess=boundaries[2][0]; exp1_guess<=boundaries[2][1]; exp1_guess+=boundaries[2][2] ){
                for(double exp2_guess=boundaries[3][0]; exp2_guess<=boundaries[3][1]; exp2_guess+=boundaries[3][2] ){
                    try{
                         
                    tuple<double, double,double,double, double, bool> result= get_interaction_4D(&p, false,{exp1_guess, exp2_guess}, {coeff1_guess, coeff2_guess});
                    
                    if(get<4>(result)<1e-5&&count<num_sol){
                        
                        double coeff1_cut=round_to_n_digits((get<0>(result)),3);
                        double coeff2_cut=round_to_n_digits((get<1>(result)),3);
                        double exp1_cut=round_to_n_digits((get<2>(result)),3);
                        double exp2_cut=round_to_n_digits((get<3>(result)),3);
                        bool found=false;
                        
                        for(unsigned int j=0; j<solution_storage.size(); j++){
                            if(abs(solution_storage[j][0]-coeff1_cut)<eps &&abs(solution_storage[j][1]-coeff2_cut)<eps&&abs(solution_storage[j][2]-exp1_cut)<eps&&abs(solution_storage[j][3]-exp2_cut)<eps){
                                found=true;
                               
                                break;
                            }
                        }
                        
                        bool validated=validate_interaction(get<4>(result), {coeff1_guess, coeff2_guess}, {exp1_guess, exp2_guess}, {(get<0>(result)), (get<1>(result))}, {abs(get<2>(result)), abs(get<3>(result))}, filename, timestamp, &p);
                        if(!validated){
                          
                            continue;
                        }
                        
                        if(!found){
                            #pragma omp critical
                            {
                            solution_storage.push_back({coeff1_cut, coeff2_cut, exp1_cut, exp2_cut});
                            count+=1;
                           
                            }
                            if(count>=num_sol){
                            
                                for(unsigned int j=0; j<solution_storage.size(); j++){
                                }
                                #pragma omp cancel for
                            }
                        }
                    }else{
                        continue;
                    }
                    }catch(int exception){
                        continue;
                    }
                    
                }
                #pragma omp cancellation point for
            }
        }
    }
}
  
    return;
}
//find interactions for a grid of critical points
void interaction_4D_crit_grid(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy, double boundaries_model  [4][3], double boundaries_crit [2][3], unsigned int terms [2],string filename, bool print, int num_sol){
    FILE * pFile;
    pFile = fopen ((filename+".csv").c_str(),"w");
    fprintf(pFile,"timestamp, ");
    for(unsigned int i=0; i<terms[0];i++){
        fprintf(pFile, "scalar_coeff_guess[%d], ",i);
    }
    for(unsigned int i=0; i<terms[0];i++){
        fprintf(pFile, "scalar_exp_guess[%d], ",i);
    }
    for(unsigned int i=0; i<terms[1];i++){
        fprintf(pFile, "vector_coeff_guess[%d], ", i);
    }
    for(unsigned int i=0; i<terms[1];i++){
        fprintf(pFile, "vector_exp_guess[%d], ",i);
    }
    for(unsigned int i=0; i<terms[0];i++){
        fprintf(pFile, "scalar_coeff[%d], ",i);
    }
    for(unsigned int i=0; i<terms[0];i++){
        fprintf(pFile, "scalar_exp[%d], ",i);
    }
    for(unsigned int i=0; i<terms[1];i++){
        fprintf(pFile, "vector_coeff[%d], ", i);
    }
    for(unsigned int i=0; i<terms[1];i++){
        fprintf(pFile, "vector_exp[%d], ",i);
    }
        
    fprintf(pFile,"GSL_error, mstar, binding_energy, binding_energy_minimum, incompress, CP_T, CP_n, CP_P; \n");
    fclose(pFile); 
    for(double crit_T=boundaries_crit[0][0]; crit_T<=boundaries_crit[0][1];crit_T+=boundaries_crit[0][2] ){
        for(double crit_n=boundaries_crit[1][0]; crit_n<=boundaries_crit[1][1];crit_n+=boundaries_crit[1][2] ){
            interaction_4D_grid( nucleon_mass, crit_T,  crit_n,  binding_energy,  saturation_density,  degeneracy,boundaries_model, terms, filename, false,num_sol);
        }
    }
    return;
}

