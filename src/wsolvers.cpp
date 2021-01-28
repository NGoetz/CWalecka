#include <../../src/wfunctions.h>
#include <../../src/wsolvers.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <matplot/matplot.h>
#include <vector>

//struct of the parameters for solving the self consistent effective mass-scalar density-effective chemical potential equations 
struct mass_eff_scalar_density_params
{
    double a;
    double number_density;
    double nucleon_mass;
    double degeneracy;
    double mass_scalar;
    double temperature;
    double b;
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
};
//struct for the parameters for finding the critical point
struct crit_params
{
    double a;
    double b;
    double nucleon_mass;
    double degeneracy;
    double mass_scalar;
    double mass_vector;
};

//single step of the root solving for the self-consistent equation of m* and scalar density at T=0; x is the vector of values to be determined,
// params the parameters and f the value of the function to be reduced to zero
int mass_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f){
        double a = ((struct mass_eff_scalar_density_params *) params)->a;
        double number_density=((struct mass_eff_scalar_density_params *) params)->number_density;
        double nucleon_mass=((struct mass_eff_scalar_density_params *) params)->nucleon_mass;
        double degeneracy=(( struct mass_eff_scalar_density_params *) params )->degeneracy;
        double mass_scalar=(( struct mass_eff_scalar_density_params *) params )->mass_scalar;
        //x0 is the prediction of the effective mass
        const double x0 = gsl_vector_get (x, 0);
        // the calculated minus the predicted value for the mass has to be zero
        const double y0 =x0-mass_eff_from_scalar_density(nucleon_mass, a, scalar_density(degeneracy,x0, number_density),mass_scalar);
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
pair<double, double> get_mass_eff_scalar_density(double degeneracy, double nucleon_mass,double a,double number_density,double mass_scalar, bool print) {

    struct mass_eff_scalar_density_params p ={a,number_density,nucleon_mass, degeneracy,mass_scalar,0,0};
    
    int status;
    size_t iter = 0;

    const size_t num = 1;
    
    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (T, num);
    
    
    gsl_multiroot_function f = {&mass_eff_scalar_density_root, num, &p};
    // this is the initialization of the effective mass
    // it should maybe be handeled centreally
    double x_init[1] = {0.6*nucleon_mass};
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
        double b = ((struct mass_eff_scalar_density_params *) params)->b;
        double nucleon_mass=((struct mass_eff_scalar_density_params *) params)->nucleon_mass;
        double number_density=((struct mass_eff_scalar_density_params *) params)->number_density;
        double temperature=((struct mass_eff_scalar_density_params *) params)->temperature;
        double degeneracy=(( struct mass_eff_scalar_density_params *) params )->degeneracy;
        double mass_scalar=(( struct mass_eff_scalar_density_params *) params )->mass_scalar;
        const double x0 = gsl_vector_get (x, 0);//prediction of effective mass
        const double x1 = gsl_vector_get (x, 1);//prediction of effective baryochemical potential
        const double y0 =x0-mass_eff_from_scalar_density(nucleon_mass,a,T_scalar_density(degeneracy,temperature,x0,x1),mass_scalar);
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
tuple<double, double,double,bool> get_T_mass_eff_mu_eff_scalar_density(double a,double b,double nucleon_mass,double number_density, double temperature, double degeneracy, double mass_scalar, bool print) {
    struct mass_eff_scalar_density_params p ={a,number_density,nucleon_mass,degeneracy, mass_scalar, temperature, b};
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
        }while (status == GSL_CONTINUE && iter < 2000);
    if (print){
        printf ("status = %s\n", gsl_strerror (status));
    
    }

    // we extract the values and calculate ns from this
    // the result is returned in a tuple
    double mass_eff=gsl_vector_get(s->x,0);
    double mu_eff=gsl_vector_get(s->x,1);
    double scalar_density=T_scalar_density(degeneracy,temperature, mass_eff, mu_eff);
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    tuple<double, double,double,bool> result={mass_eff,mu_eff,scalar_density,status==GSL_SUCCESS};
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
    
    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    // we have to first calculate the effective mass with the
    // guess we have for the coefficients
    pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,x0,saturation_density,mass_scalar, false);
    // here,  the minimum condition for stable nuclei is given
    const double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, x0, x1, output.second, mass_vector, mass_scalar);
    const double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,x0,x1, output.second, mass_vector, mass_scalar);
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
pair<double, double> get_coeff(double nucleon_mass, double saturation_density, double binding_energy, double degeneracy, double mass_scalar, double mass_vector, bool print){

    struct coeff_params p ={saturation_density,nucleon_mass, binding_energy,degeneracy,mass_scalar,mass_vector};
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
    while (status == GSL_CONTINUE && iter < 1000);
    if (print){
        printf ("status = %s\n", gsl_strerror (status));
    }
    double a=gsl_vector_get(s->x,0);
    double b=gsl_vector_get(s->x,1);
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    pair<double, double> result={a,b};
    return result;
}

int crit_root(const gsl_vector * x, void *params, gsl_vector * f){
    double a=((struct crit_params *) params)->a;
    double b=((struct crit_params *) params)->b;
    double nucleon_mass=((struct crit_params *) params)->nucleon_mass;
    double degeneracy=((struct crit_params *) params)->degeneracy;
    double mass_vector=((struct crit_params *) params)->mass_vector;
    double mass_scalar=((struct crit_params *) params)->mass_scalar;
    
    const double x0=abs(gsl_vector_get(x,0));//number density
    const double x1=abs(gsl_vector_get(x,1));//temperature - both need to be positive
    
    // find the inflection point
    const double y0 =T_press_dn_solv(degeneracy, nucleon_mass, x0, mass_vector, mass_scalar, a, b, x1);
    const double y1 =T_press_dn2_solv(degeneracy, nucleon_mass, x0, mass_vector, mass_scalar, a, b, x1);
    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);
    return GSL_SUCCESS;
}
//find the critical point
pair<double, double> get_crit(double nucleon_mass,double saturation_density, double binding_energy, double mass_scalar, double mass_vector, double degeneracy, bool print){
    pair<double, double> coeff=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    struct crit_params p ={coeff.first, coeff.second, nucleon_mass, degeneracy, mass_scalar, mass_vector};
    int status;
    size_t iter = 0;

    const size_t num = 2;
    
    const gsl_multiroot_fsolver_type * Q = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (Q, num);
    
    
    gsl_multiroot_function f = {&crit_root, num, &p};

    double x_init[num] = {saturation_density,10};
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

// solve the energy density at T=0 with given coupling constants
double energy_pp_minus_mass_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a, number_density, mass_scalar, false);
    return eps_over_n(degeneracy, x.first,number_density,a,b,x.second,mass_vector,mass_scalar)-nucleon_mass;
}

// derivative of this
double energy_pp_minus_mass_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a, number_density, mass_scalar, false);
    return eps_over_n_dn(degeneracy, x.first,number_density,a,b,x.second,mass_vector,mass_scalar);
}
//same for pressure as T=0
double press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a, number_density, mass_scalar, false);
    return press(degeneracy, x.first,number_density,a,b,x.second,mass_vector,mass_scalar);
}
//incompressibility at T=0
double incsolv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b){
    pair<double, double> x=get_mass_eff_scalar_density(degeneracy, nucleon_mass, a, number_density, mass_scalar, false);
    return incompress(degeneracy, x.first,number_density,a,b,x.second,mass_vector,mass_scalar);
}
//pressure at finite T
double T_press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double b, double temperature){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    return T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),a, b,get<2>(x), number_density,mass_vector, mass_scalar);
}
//derivative of pressure at finite T after mu
double  T_press_dmu_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double b, double temperature, double delta){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    return (T_pressure(degeneracy, temperature,get<0>(x),get<1>(x)+delta,a, b,get<2>(x), number_density,mass_vector, mass_scalar)-T_pressure(degeneracy, temperature,get<0>(x),get<1>(x)-delta,a, b,get<2>(x), number_density,mass_vector, mass_scalar))/2*delta;
}
//derivative of pressure at finite T after T
double  T_press_dT_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density, temperature+delta, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density, temperature-delta, degeneracy,  mass_scalar, false);
    return (T_pressure(degeneracy, temperature+delta,get<0>(x1),get<1>(x1),a, b,get<2>(x1), number_density,mass_vector, mass_scalar)-T_pressure(degeneracy, temperature-delta,get<0>(x2),get<1>(x2),a, b,get<2>(x2), number_density,mass_vector, mass_scalar))/2*delta;
}
//derivative of pressure at finite T after the number density
double  T_press_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double b, double temperature, double delta){
    
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density+delta, temperature, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density-delta, temperature, degeneracy,  mass_scalar, false);
    return (T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),a, b,get<2>(x1), number_density+delta,mass_vector, mass_scalar)-T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),a, b,get<2>(x2), number_density-delta,mass_vector, mass_scalar))/2*delta;
}
//second derivative of pressure at finite T after the number density
double  T_press_dn2_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double b, double temperature, double delta){
    tuple<double, double,double,bool>x1= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density+delta, temperature, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    tuple<double, double,double,bool>x2= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density-delta, temperature, degeneracy,  mass_scalar, false);
    return (T_pressure(degeneracy, temperature,get<0>(x1),get<1>(x1),a, b,get<2>(x1), number_density+delta,mass_vector, mass_scalar)+T_pressure(degeneracy, temperature,get<0>(x2),get<1>(x2),a, b,get<2>(x2), number_density-delta,mass_vector, mass_scalar)-2*T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),a, b,get<2>(x), number_density,mass_vector, mass_scalar))/delta*delta;
}
// energy density at finite T
double T_eps_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double b, double temperature){
    tuple<double, double,double,bool>x= get_T_mass_eff_mu_eff_scalar_density(a,b,nucleon_mass,number_density, temperature, degeneracy,  mass_scalar, false);
    double mu=get<1>(x)+4*M_PI*b/pow(mass_vector,2.0)*number_density;//get mu from n
    return -T_pressure(degeneracy, temperature,get<0>(x),get<1>(x),a, b,get<2>(x), number_density,mass_vector, mass_scalar)+mu*number_density+temperature*T_press_dT_solv(degeneracy, nucleon_mass, number_density, mass_vector, mass_scalar, a,b, temperature);
}


