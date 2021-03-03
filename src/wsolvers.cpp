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


struct mass_eff_scalar_density_params
{
    double a;
    double exp_a;
    double number_density;
    double nucleon_mass;
    double degeneracy;
    double mass_scalar;
    double temperature;
    double b;
    double exp_b;
};
//struct for the parameters for root solving for the coupling constants
struct coeff_params
{
    double saturation_densitiy;
    double nucleon_mass;
    double binding_energy;
    double degeneracy;
    double mass_scalar;
    double mass_vector;
    double exp_a;
    double exp_b;
};
//struct for the parameters for finding the critical point
struct crit_params
{
    double a;
    double exp_a;
    double b;
    double exp_b;
    double nucleon_mass;
    double degeneracy;
    double mass_scalar;
    double mass_vector;
};
//struct for the parameters for finding the interaction terms
struct interaction_params
{
    double nucleon_mass;
    double degeneracy;
    double saturation_density;
    double binding_energy;
    double critical_temperature;
    double critical_density;
    double mass_scalar;
    double mass_vector;
};

//single step of the root solving for the self-consistent equation of m* and scalar density at T=0; x is the vector of values to be determined,
// params the parameters and f the value of the function to be reduced to zero
int mass_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f){
        double a = ((struct mass_eff_scalar_density_params *) params)->a;
        double exp_a = ((struct mass_eff_scalar_density_params *) params)->exp_a;
        double number_density=((struct mass_eff_scalar_density_params *) params)->number_density;
        double nucleon_mass=((struct mass_eff_scalar_density_params *) params)->nucleon_mass;
        double degeneracy=(( struct mass_eff_scalar_density_params *) params )->degeneracy;
        double mass_scalar=(( struct mass_eff_scalar_density_params *) params )->mass_scalar;
        //x0 is the prediction of the effective mass
        const double x0 = gsl_vector_get (x, 0);
        // the calculated minus the predicted value for the mass has to be zero
        const double y0 =x0-mass_eff_from_scalar_density(nucleon_mass, a, exp_a, scalar_density(degeneracy,x0, number_density),mass_scalar);
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
pair<double, double> get_mass_eff_scalar_density(double degeneracy, double nucleon_mass,double a,double exp_a,double number_density,double mass_scalar, bool print) {

    struct mass_eff_scalar_density_params p ={a,exp_a,number_density,nucleon_mass, degeneracy,mass_scalar,0,0};
    
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
        double a = ((struct mass_eff_scalar_density_params *) params)->a;
        double exp_a = ((struct mass_eff_scalar_density_params *) params)->exp_a;
        double nucleon_mass=((struct mass_eff_scalar_density_params *) params)->nucleon_mass;
        double number_density=((struct mass_eff_scalar_density_params *) params)->number_density;
        double temperature=((struct mass_eff_scalar_density_params *) params)->temperature;
        double degeneracy=(( struct mass_eff_scalar_density_params *) params )->degeneracy;
        double mass_scalar=(( struct mass_eff_scalar_density_params *) params )->mass_scalar;
        const double x0 = abs(gsl_vector_get (x, 0));
        const double x1 = abs(gsl_vector_get (x, 1));//prediction of effective baryochemical potential
        const double y0 =x0-mass_eff_from_scalar_density(nucleon_mass,a,exp_a,T_scalar_density(degeneracy,temperature,x0,x1),mass_scalar);
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
tuple<double, double,double,bool> get_T_mass_eff_mu_eff_scalar_density(double a, double exp_a, double b, double exp_b, double nucleon_mass,double number_density, double temperature, double degeneracy, double mass_scalar, bool print) {
    struct mass_eff_scalar_density_params p ={a,exp_a,number_density,nucleon_mass,degeneracy, mass_scalar, temperature, b,exp_b};
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
// the coefficients are g**2/4*M_PI
int coeff_root(const gsl_vector * x, void *params, gsl_vector * f){
    double saturation_density=((struct coeff_params *) params)->saturation_densitiy;
    double nucleon_mass=((struct coeff_params *) params)->nucleon_mass;
    double binding_energy=((struct coeff_params *) params)->binding_energy;
    double degeneracy=((struct coeff_params *) params)->degeneracy;
    double mass_scalar=((struct coeff_params *) params)->mass_scalar;
    double mass_vector=((struct coeff_params *) params)->mass_vector;
    double exp_a=((struct coeff_params *) params)->exp_a;
    double exp_b=((struct coeff_params *) params)->exp_b;
    
    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    // we have to first calculate the effective mass with the
    // guess we have for the coefficients
    pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,x0, exp_a, saturation_density,mass_scalar, false);
    // here,  the minimum condition for stable nuclei is given
    const double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, x0, exp_a, x1, exp_b, output.second, mass_vector, mass_scalar);
    const double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,x0, exp_a, x1, exp_b, output.second, mass_vector, mass_scalar);
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
pair<double, double> get_coeff(double nucleon_mass, double saturation_density, double binding_energy, double degeneracy, double mass_scalar, double mass_vector, double exp_a, double exp_b, bool print, bool * success){

    struct coeff_params p ={saturation_density,nucleon_mass, binding_energy,degeneracy,mass_scalar,mass_vector, exp_a, exp_b};
    int status;
    size_t iter = 0;

    const size_t num = 2;
    
    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (T, num);
    
    
    gsl_multiroot_function f = {&coeff_root, num, &p};

    double x_init[2] = {10,10};
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
    double a=((struct crit_params *) params)->a;
    double exp_a=((struct crit_params *) params)->exp_a;
    double b=((struct crit_params *) params)->b;
    double exp_b=((struct crit_params *) params)->exp_b;
    double nucleon_mass=((struct crit_params *) params)->nucleon_mass;
    double degeneracy=((struct crit_params *) params)->degeneracy;
    double mass_vector=((struct crit_params *) params)->mass_vector;
    double mass_scalar=((struct crit_params *) params)->mass_scalar;
    
    const double x0=abs(gsl_vector_get(x,0));//number density
    const double x1=abs(gsl_vector_get(x,1));//temperature - both need to be positive
    
    // find the inflection point
    const double y0 =T_press_dn_solv(degeneracy, nucleon_mass, x0, mass_vector, mass_scalar, a, exp_a, b, exp_b, x1);
    const double y1 =T_press_dn2_solv(degeneracy, nucleon_mass, x0, mass_vector, mass_scalar, a, exp_a,b, exp_b, x1);
    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);
    return GSL_SUCCESS;
}
//find the critical point
pair<double, double> get_crit(double nucleon_mass,double saturation_density, double binding_energy, double mass_scalar, double mass_vector, double degeneracy, double exp_a, double exp_b, bool print, double crit_T, double crit_n, double a, double b){
    struct crit_params p ;
    //if no coefficients are given, the model was only determined from the saturation point
    if(a==0.0&&b==0.0){
        pair<double, double> coeff=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, exp_a, exp_b, false);
        p ={coeff.first, exp_a, coeff.second, exp_b, nucleon_mass, degeneracy, mass_scalar, mass_vector};
    }else{
         p ={a, exp_a, b, exp_b, nucleon_mass, degeneracy, mass_scalar, mass_vector};
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

int interaction_root(const gsl_vector * x, void *params, gsl_vector * f){
    double nucleon_mass=((struct interaction_params *) params)->nucleon_mass;
    double degeneracy=((struct interaction_params *) params)->degeneracy;
    double mass_vector=((struct interaction_params *) params)->mass_vector;
    double mass_scalar=((struct interaction_params *) params)->mass_scalar;
    double saturation_density=((struct interaction_params *) params)->saturation_density;
    double binding_energy=((struct interaction_params *) params)->binding_energy;
    double critical_temperature=((struct interaction_params *) params)->critical_temperature;
    double critical_density=((struct interaction_params *) params)->critical_density;
    const double x0=gsl_vector_get(x,0);//coefficient a (scalar term)
    const double x1=(gsl_vector_get(x,1));//coefficient b (vector term)
    const double x2=gsl_vector_get(x,2);//exponent for scalar term
    const double x3=gsl_vector_get(x,3);//exponent for vector term
    pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,x0, x2, saturation_density,mass_scalar, false);
    //predict binding at demanded saturation density
    const double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, x0, x2, x1, x3, output.second, mass_vector, mass_scalar);
    const double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,x0, x2, x1, x3, output.second, mass_vector, mass_scalar);
    //there is an inflection point at the critical temperature/density
    pair<double, double> pressures=T_press_dn12_solv( degeneracy,  nucleon_mass, critical_density,  mass_vector,  mass_scalar, x0, x2, x1, x3,  critical_temperature);
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

// solve the energy density at T=0 with given coupling constants
double energy_pp_minus_mass_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double exp_a, double b, double exp_b){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a, exp_a, number_density, mass_scalar, false);
    return eps_over_n(degeneracy, x.first,number_density,a, exp_a, b, exp_b, x.second,mass_vector,mass_scalar)-nucleon_mass;
}

// derivative of this
double energy_pp_minus_mass_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a, exp_a, number_density, mass_scalar, false);
    return eps_over_n_dn(degeneracy, x.first,number_density,a,exp_a, b, exp_b, x.second,mass_vector,mass_scalar);
}
//same for pressure as T=0
double press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a,exp_a, number_density, mass_scalar, false);
    return press(degeneracy, x.first,number_density,a, exp_a, b, exp_b,x.second,mass_vector,mass_scalar);
}
//incompressibility at T=0
double incsolv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double delta){
    pair<double, double> x1=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a,exp_a, number_density+delta, mass_scalar, false);
    pair<double, double> x2=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a,exp_a, number_density-delta, mass_scalar, false);
    return 9*number_density*number_density*(eps_over_n_dn(degeneracy, x1.first,number_density+delta, a, exp_a, b, exp_b, x1.second, mass_vector,mass_scalar)-eps_over_n_dn(degeneracy, x2.first,number_density-delta, a, exp_a, b, exp_b, x2.second, mass_vector,mass_scalar))/(2*delta);
}
//pressure at finite T
double T_press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double exp_a, double b, double exp_b, double temperature){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    return T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),a, exp_a, b, exp_b,get<2>(x), number_density,mass_vector, mass_scalar);
}
//derivative of pressure at finite T after mu
double  T_press_dmu_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double exp_a, double b, double exp_b, double temperature, double delta){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a,exp_a, b, exp_b, nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    return (T_pressure(degeneracy, temperature,get<0>(x),get<1>(x)+delta,a, exp_a, b, exp_b,get<2>(x), number_density,mass_vector, mass_scalar)-T_pressure(degeneracy, temperature,get<0>(x),get<1>(x)-delta,a, exp_a, b, exp_b,get<2>(x), number_density,mass_vector, mass_scalar))/2*delta;
}
//derivative of pressure at finite T after T
double  T_press_dT_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density, temperature+delta, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density, temperature-delta, degeneracy,  mass_scalar, false);
    return (T_pressure(degeneracy, temperature+delta,get<0>(x1),get<1>(x1),a, exp_a, b, exp_b,get<2>(x1), number_density,mass_vector, mass_scalar)-T_pressure(degeneracy, temperature-delta,get<0>(x2),get<1>(x2),a, exp_a, b, exp_b,get<2>(x2), number_density,mass_vector, mass_scalar))/2*delta;
}
//derivative of pressure at finite T after the number density
double  T_press_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density+delta, temperature, degeneracy,  mass_scalar, false);
    
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density-delta, temperature, degeneracy,  mass_scalar, false);
    
    return (T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),a, exp_a, b, exp_b,get<2>(x1), number_density+delta,mass_vector, mass_scalar)-T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),a, exp_a, b, exp_b,get<2>(x2), number_density-delta,mass_vector, mass_scalar))/(2*delta);
}
//second derivative of pressure at finite T after the number density
double  T_press_dn2_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(a,exp_a, b, exp_b,nucleon_mass,number_density+delta, temperature, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density-delta, temperature, degeneracy,  mass_scalar, false);
    return (T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),a, exp_a, b, exp_b,get<2>(x1), number_density+delta,mass_vector, mass_scalar)+T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),a, exp_a, b, exp_b,get<2>(x2), number_density-delta,mass_vector, mass_scalar)-2*T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),a, exp_a, b, exp_b,get<2>(x), number_density,mass_vector, mass_scalar))/(delta*delta);
}
// energy density at finite T
double T_eps_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double exp_a, double b, double exp_b, double temperature){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    double mu=get<1>(x)+4*M_PI*b/pow(mass_vector,2.0)*number_density;//get mu from n
    return -T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),a, exp_a, b, exp_b,get<2>(x), number_density,mass_vector, mass_scalar)+mu*number_density+temperature*T_press_dT_solv(degeneracy, nucleon_mass, number_density, mass_vector, mass_scalar, a, exp_a, b, exp_b, temperature);
}
//get both the first and second derivate of the pressure at the same time
pair<double, double> T_press_dn12_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(a,exp_a, b, exp_b,nucleon_mass,number_density+delta, temperature, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(a, exp_a, b, exp_b,nucleon_mass,number_density-delta, temperature, degeneracy,  mass_scalar, false);
    double pressure1=T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),a, exp_a, b, exp_b,get<2>(x1), number_density+delta,mass_vector, mass_scalar);
    double pressure2=T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),a, exp_a, b, exp_b,get<2>(x2), number_density-delta,mass_vector, mass_scalar);
    double pressure=T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),a, exp_a, b, exp_b,get<2>(x), number_density,mass_vector, mass_scalar);
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
bool validate_interaction(double error, double a_guess, double b_guess, double exp_a_guess, double exp_b_guess, double a_in, double b_in, double exp_a, double exp_b, string filename, int timestamp, void* params){
    double nucleon_mass=((struct interaction_params *) params)->nucleon_mass;
    double degeneracy=((struct interaction_params *) params)->degeneracy;
    double mass_vector=((struct interaction_params *) params)->mass_vector;
    double mass_scalar=((struct interaction_params *) params)->mass_scalar;
    double saturation_density=((struct interaction_params *) params)->saturation_density;
    double binding_energy=((struct interaction_params *) params)->binding_energy;
    double critical_temperature=((struct interaction_params *) params)->critical_temperature;
    double critical_density=((struct interaction_params *) params)->critical_density;
   
    pair<double, double> mstarns=get_mass_eff_scalar_density( degeneracy,  nucleon_mass, a_in,  exp_a,  saturation_density, mass_scalar,false);
    
    double new_binding_energy=nucleon_mass-eps_over_n(degeneracy, mstarns.first, saturation_density, a_in, exp_a, b_in, exp_b, mstarns.second, mass_vector, mass_scalar);
    double binding_energy_minimum =eps_over_n_dn(degeneracy, mstarns.first, saturation_density,a_in, exp_a, b_in, exp_b, mstarns.second, mass_vector, mass_scalar);
    double incompress=incsolv( degeneracy,  nucleon_mass,saturation_density,  mass_vector,  mass_scalar,  a_in,  exp_a,  b_in,  exp_b);
    
    pair<double, double> CP= get_crit( nucleon_mass, saturation_density,  binding_energy,  mass_scalar,  mass_vector,  degeneracy,  exp_a,  exp_b, false, 19,0.44*saturation_density,a_in,b_in);
    
    
    double CP_P=T_press_solv(degeneracy, nucleon_mass,CP.second, mass_vector, mass_scalar, a_in, exp_a, b_in, exp_b, CP.first)/pow(197.3,3);
    
    if(abs(CP.first-critical_temperature)>1e-2||abs(CP.second/saturation_density-critical_density/saturation_density)>1e-3){
        return false;
    }
    FILE * pFile;
    pFile = fopen ((filename+".csv").c_str(),"a");
    fprintf(pFile, "%d, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.5e, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16F; \n ",timestamp, a_guess, b_guess, exp_a_guess, exp_b_guess, a_in, b_in, exp_a, exp_b, error, mstarns.first, new_binding_energy, binding_energy_minimum, incompress, CP.first, CP.second/saturation_density, CP_P);
   
    fclose(pFile);
    return true;
    
}
//find the coefficients and exponents of the saclar and vector term
tuple<double, double,double,double, double, bool> get_interaction_4D(void * p, bool print, double int_exp_a, double int_exp_b, double int_a, double int_b) {
     gsl_set_error_handler(&gsl_handler);
    int status;
    size_t iter = 0;
    const size_t num = 4;
    gsl_vector *x = gsl_vector_alloc (num);
    const gsl_multiroot_fsolver_type * K = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (K, num);
  
    gsl_multiroot_function f = {&interaction_root, num, p};
    double x_init[4] = {int_a,int_b,int_exp_a,int_exp_b};
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
//find interactions on a grid of initial values
void interaction_4D_grid(double nucleon_mass,double critical_temperature, double critical_density, double binding_energy, double saturation_density, double degeneracy, double mass_scalar, double mass_vector, double boundaries  [4][3],string filename, bool print, int num_sol){
    gsl_set_error_handler(&gsl_handler);
    int timestamp=time(0);
    if(num_sol==0){
        FILE * pFile;
        pFile = fopen ((filename+".csv").c_str(),"w");
        fprintf(pFile,"timestamp, guess_a, guess_b, guess_exp_a, guess_exp_b, a, b, exp_a, exp_b, GSL_error, mstar, binding_energy, binding_energy_minimum, incompress, CP_T, CP_n, CP_P; \n");
        fclose(pFile);
    }
    struct interaction_params p={nucleon_mass, degeneracy, saturation_density,binding_energy, critical_temperature, critical_density, mass_scalar, mass_vector };
    double min, min_a, min_b, min_exp_a, min_exp_b, guess_min_a, guess_min_b, guess_min_exp_a, guess_min_exp_b;
    min=1e20;
    int count=0;
    vector<vector<double>> solution_storage;
    int num_a_guess=ceil((boundaries[0][1]-boundaries[0][0])/boundaries[0][2]);
    int num_threads=ceil(omp_get_max_threads()-1);
    
#pragma omp parallel num_threads(num_threads) shared(p, solution_storage, num_sol, count)
{
#pragma omp for
    for(int i=0; i<=num_a_guess; i+=1 ){
        double a_guess=boundaries[0][0]+i*boundaries[0][2];
        
        for(double b_guess=boundaries[1][0]; b_guess<=boundaries[1][1]; b_guess+=boundaries[1][2] ){

            for(double exp_a_guess=boundaries[2][0]; exp_a_guess<=boundaries[2][1]; exp_a_guess+=boundaries[2][2] ){
                for(double exp_b_guess=boundaries[3][0]; exp_b_guess<=boundaries[3][1]; exp_b_guess+=boundaries[3][2] ){
                    try{
                         
                    tuple<double, double,double,double, double, bool> result= get_interaction_4D(&p, false,exp_a_guess,exp_b_guess,a_guess, b_guess);
                    
                    if(get<4>(result)<1e-5&&count<num_sol){
                        
                        bool validated=validate_interaction(get<4>(result), a_guess, b_guess,  exp_a_guess,  exp_b_guess, (get<0>(result)), (get<1>(result)), abs(get<2>(result)), abs(get<3>(result)), filename, timestamp, &p);
                        if(!validated){
                          
                            continue;
                        }
                        double a_cut=round((get<0>(result))*100)/100.0;
                        double b_cut=round((get<1>(result))*100)/100.0;
                        double exp_a_cut=round(abs(get<2>(result))*100)/100.0;
                        double exp_b_cut=round(abs(get<3>(result))*100)/100.0;
                        bool found=false;
                        
                        for(int j=0; j<solution_storage.size(); j++){
                            if(abs(solution_storage[j][0]-a_cut)<1e-1 &&abs(solution_storage[j][1]-b_cut)<1e-1&&abs(solution_storage[j][2]-exp_a_cut)<1e-2&&abs(solution_storage[j][3]-exp_b_cut)<1e-2){
                                found=true;
                               
                                break;
                            }
                        }
                        if(!found){
                            #pragma omp critical
                            {
                            solution_storage.push_back({a_cut,b_cut,exp_a_cut, exp_b_cut});
                            count+=1;
                           
                            }
                            if(count>=num_sol){
                            
                                cout<<solution_storage.size()<<"\n";
                                for(int j=0; j<solution_storage.size(); j++){
                                }
                                #pragma omp cancel for
                            }
                        }
                    }else{
                    }
                    if(get<4>(result)<min){
                        min=get<4>(result);
                        min_a=get<0>(result);
                        min_b=get<1>(result);
                        min_exp_a=get<2>(result);
                        min_exp_b=get<3>(result);
                        guess_min_a=a_guess;
                        guess_min_b=b_guess;
                        guess_min_exp_a=exp_a_guess;
                        guess_min_exp_b=exp_b_guess;
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
    cout<<"error "<<min<<" for guess (a,b, exp_a, exp_b) "<<guess_min_a<<" "<<guess_min_b<<" "<<guess_min_exp_a<<" "<<guess_min_exp_b<<" and with gsl output "<<min_a<<" "<<min_b<<" "<<min_exp_a<<" "<<min_exp_b<<"\n";
    return;
}
//find interactions for a grid of critical points
void interaction_4D_crit_grid(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy, double mass_scalar, double mass_vector, double boundaries_model  [4][3], double boundaries_crit [2][3], string filename, bool print, int num_sol){
    FILE * pFile;
    pFile = fopen ((filename+".csv").c_str(),"w");
    fprintf(pFile,"timestamp, guess_a, guess_b, guess_exp_a, guess_exp_b, a, b, exp_a, exp_b, GSL_error, mstar, binding_energy, binding_energy_minimum, incompress, CP_T, CP_n, CP_P; \n");
    fclose(pFile);
    for(double crit_T=boundaries_crit[0][0]; crit_T<=boundaries_crit[0][1];crit_T+=boundaries_crit[0][2] ){
        for(double crit_n=boundaries_crit[1][0]; crit_n<=boundaries_crit[1][1];crit_n+=boundaries_crit[1][2] ){
            interaction_4D_grid( nucleon_mass, crit_T,  crit_n,  binding_energy,  saturation_density,  degeneracy,  mass_scalar,  mass_vector,boundaries_model,filename, print,num_sol);
        }
    }
    return;
}

