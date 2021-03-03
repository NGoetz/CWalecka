#include <wfunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>

using namespace std;

//holds the parameters of integrations for finite temperature 
struct integration_params
{
    double temperature;
    double mass_eff;
    double mu_eff;
    double degeneracy;
};


//calculates the momentum from density
double fmom(double degeneracy, double number_density){
    return pow((((6.0*pow(M_PI,2.0))/degeneracy)*number_density),(1.0/3.0));
}
//calculates density from fermi momentum at T=0
double inv_p(double degeneracy, double p){
    return pow(p,(3.0))*(degeneracy/(6.0*pow(M_PI,2.0)));
}
//calculates single particle energy from momentum
double E_eff(double p, double mass_eff){
    return pow((pow(p,2.0)+pow(mass_eff,2.0)),(1.0/2.0));
}
//the scalar density at T=0
double scalar_density(double degeneracy, double mass_eff, double number_density){
    return degeneracy*(mass_eff/(4*pow(M_PI,2.0)))*(E_eff(fmom(degeneracy,number_density),mass_eff)*fmom(degeneracy,number_density)-0.5*pow(mass_eff,2.0)*log(((1+fmom(degeneracy,number_density)/E_eff(fmom(degeneracy,number_density),mass_eff))/(1-fmom(degeneracy, number_density)/E_eff(fmom(degeneracy,number_density),mass_eff)))));
}
// its derivative after the number density
double scalar_density_dn(double degeneracy, double mass_eff, double number_density, double a, double exp_a, double mass_scalar){
    double scalar_density_val=scalar_density(degeneracy, mass_eff, number_density);
    return (2.0*mass_eff*pow(mass_scalar,2.0)*M_PI)/(12.0*a*number_density*(exp_a-1.0)*pow(scalar_density_val, exp_a-2.0)*pow(M_PI,2.0) + 6.0*a*(exp_a-1.0)*pow(scalar_density_val, exp_a-2.0)*degeneracy*pow(mass_eff,2.0)*fmom(degeneracy,number_density) + 
    (2.0*pow(mass_scalar,2.0)*M_PI - 3*a*(exp_a-1.0)*pow(scalar_density_val, exp_a-2.0)*degeneracy*pow(mass_eff,2.0)*log(1.0 + 12.0/((-6.0*pow(number_density/degeneracy,1.0/3.0) + (pow(6.0/M_PI, 2.0/3.0)*E_eff(fmom(degeneracy,number_density), mass_eff)))/pow(number_density/degeneracy,1.0/3.0))))*E_eff(fmom(degeneracy,number_density), mass_eff));
}

//the effective mass from scalar density
double mass_eff_from_scalar_density(double nucleon_mass, double a, double exp_a, double scalar_density, double mass_scalar){
    double renorm=((4*M_PI*a)/pow(mass_scalar,2))*pow(scalar_density,exp_a-1.0);
    if (renorm>nucleon_mass){
        return 0.0;
    }
    return nucleon_mass-renorm;
}
// its derivative after n
double mass_eff_from_scalar_density_dn(double mass_eff,  double number_density,  double a, double exp_a, double mass_scalar, double degeneracy){
    double scalar_density_val=scalar_density(degeneracy, mass_eff, number_density);
    return -((4*M_PI*a)/pow(mass_scalar,2))*scalar_density_dn(degeneracy, mass_eff, number_density, a, exp_a, mass_scalar)*(exp_a-1.0)*pow(scalar_density_val,exp_a-2.0);
}

//the energy density at T=0 divided by the density.
double eps_over_n(double degeneracy, double mass_eff, double number_density, double a, double exp_a, double b, double exp_b, double scalar_density, double mass_vector, double mass_scalar){
    return(((degeneracy/(16.0*pow(M_PI,2)))*((2.0)*E_eff(fmom(degeneracy,number_density), mass_eff)*(pow(fmom(degeneracy,number_density),3)) + (pow(mass_eff,2))*(E_eff(fmom(degeneracy,number_density),mass_eff)*fmom(degeneracy,number_density) + pow(mass_eff,2)*log(mass_eff/(E_eff(fmom(degeneracy,number_density), mass_eff) + fmom(degeneracy,number_density))))) + ((exp_b-1)/(exp_b))*(4*M_PI*b)/(pow(mass_vector,2))*pow(number_density,exp_b) + ((exp_a-1)/(exp_a))*(4*M_PI*a)/(pow(mass_scalar,2))*(pow(scalar_density,exp_a)))/(number_density));
}

//the derivate of the energy density divided by density at T=0. 
double eps_over_n_dn(double degeneracy, double mass_eff, double number_density, double a, double exp_a, double b, double exp_b, double scalar_density, double mass_vector, double mass_scalar){
    
    double scalar_density_dn_val = scalar_density_dn(degeneracy, mass_eff, number_density, a, exp_a, mass_scalar);
    double mass_eff_dn_val = mass_eff_from_scalar_density_dn( mass_eff,  number_density,  a,exp_a, mass_scalar, degeneracy);
   
    return ((1.0/(16.0*pow(number_density*M_PI,2.0)*E_eff(fmom(degeneracy,number_density), mass_eff) ))*(-2.0*number_density*pow(M_PI*mass_eff,2.0) - degeneracy*fmom(degeneracy,number_density)*pow(mass_eff,4.0)+ 4.0*number_density*pow(fmom(degeneracy,number_density)*M_PI,2.0) + 
   16.0*(pow(exp_b-1,2.0)/exp_b)*pow(number_density,exp_b)*pow(M_PI,2.0)*4*M_PI*(b/pow(mass_vector,2.0))* E_eff(fmom(degeneracy,number_density), mass_eff) + 
   24.0*pow(number_density*M_PI,2.0)*mass_eff*mass_eff_dn_val + 4*degeneracy*fmom(degeneracy,number_density)*number_density*pow(mass_eff,3.0)*mass_eff_dn_val  - degeneracy*log(mass_eff/(fmom(degeneracy,number_density) + E_eff(fmom(degeneracy,number_density), mass_eff)))*pow(mass_eff,3.0)* E_eff(fmom(degeneracy,number_density), mass_eff)* (mass_eff - 4*number_density*mass_eff_dn_val) + 16.0*((exp_a-1.0)/exp_a)*E_eff(fmom(degeneracy,number_density), mass_eff)*4*pow(M_PI,3.0)*(a/pow(mass_scalar,2.0))*pow(scalar_density, exp_a-1.0)* (-scalar_density + exp_a*number_density* scalar_density_dn_val)));
 
}
//the pressure at T=0. 
double press(double degeneracy, double mass_eff, double number_density, double a, double exp_a, double b, double exp_b, double scalar_density, double mass_vector, double mass_scalar){
    return(((degeneracy/(16.0*pow(M_PI,2)))*((2.0/3.0)*E_eff(fmom(degeneracy,number_density),mass_eff)*(pow(fmom(degeneracy,number_density),3)) - (pow(mass_eff,2))*(E_eff(fmom(degeneracy,number_density),mass_eff)*fmom(degeneracy,number_density) + pow(mass_eff,2)*log(mass_eff/(E_eff(fmom(degeneracy,number_density),mass_eff) + fmom(degeneracy,number_density))))) + ((exp_b-1)/(exp_b))*(4*M_PI*b)/(pow(mass_vector,2))*pow(number_density,exp_b) - ((exp_a-1)/(exp_a))*(4*M_PI*a)/(pow(mass_scalar,2))*(pow(scalar_density,exp_a))));
}
//Finite temperature properties
//fermi distribution with -mu
double fermiminus(double temperature, double mass_eff, double p, double mu_eff){
    return 1.0/(exp((1/temperature)*(E_eff(p,mass_eff)-mu_eff))+1.0);
}
//fermi distribution with -mu
double fermiplus(double temperature, double mass_eff, double p, double mu_eff){
    return 1.0/(exp((1/temperature)*(E_eff(p,mass_eff)+mu_eff))+1.0);
}
// chemical potential from density
double muB_from_number_density(double degeneracy, double mass_eff, double number_density, double mass_vector, double b, double exp_b){
    return sqrt(pow(mass_eff,2.0)+pow(fmom(degeneracy,number_density),2.0))+((4*M_PI*b)/pow(mass_vector,2))*pow(number_density,exp_b-1);
}
//the effective chemical potential at finite temperature (at zero, it is just the energy)
double T_mu_eff(double number_density, double muB, double b, double exp_b, double mass_vector){
    return muB-((4*M_PI*b)/pow(mass_vector,2))*pow(number_density,exp_b-1);
}
//integrand of pressure at finite T, to be integrated over momentum p
double T_pressure_integrand (double p, void * params){
    double temperature = ((struct integration_params *) params)->temperature;
    double mass_eff=((struct integration_params *) params)->mass_eff;
    double mu_eff=((struct integration_params *) params)->mu_eff;
    double degeneracy=((struct integration_params *) params)->degeneracy;
    return (degeneracy*temperature/(2*pow(M_PI,2.0))*pow(p,2.0)*(log(1/fermiminus(-temperature,mass_eff,p,mu_eff))+log(1/fermiplus(-temperature,mass_eff,p,mu_eff))));
}
//integrand of scalar density at finite T, to be integrated over momentum p
double T_scalar_density_integrand (double p, void * params){
    double temperature = ((struct integration_params *) params)->temperature;
    double mass_eff=((struct integration_params *) params)->mass_eff;
    double mu_eff=((struct integration_params *) params)->mu_eff;
    double degeneracy=((struct integration_params *) params)->degeneracy;
    return (mass_eff*degeneracy/(2.0*pow(M_PI,2.0))*pow(p,2.0)/E_eff(p,mass_eff)*fermiminus(temperature,mass_eff,p,mu_eff)+mass_eff*degeneracy/(2.0*pow(M_PI,2.0))*pow(p,2.0)/E_eff(p,mass_eff)*fermiplus(temperature,mass_eff,p,mu_eff));
}
//integrand of number density at finite T, to be integrated over momentum p
double T_number_density_integrand (double p, void * params){
    double temperature = ((struct integration_params *) params)->temperature;
    double mass_eff=((struct integration_params *) params)->mass_eff;
    double mu_eff=((struct integration_params *) params)->mu_eff;
    double degeneracy=((struct integration_params *) params)->degeneracy;
    return (degeneracy/(2.0*pow(M_PI,2.0))*pow(p,2.0)*fermiminus(temperature,mass_eff,p,mu_eff)-degeneracy/(2.0*pow(M_PI,2.0))*pow(p,2.0)*fermiplus(temperature,mass_eff,p,mu_eff));
}
//pressure at finite temperature
double T_pressure(double degeneracy, double temperature, double mass_eff, double mu_eff, double a, double exp_a, double b, double exp_b, double scalar_density, double number_density, double mass_vector, double mass_scalar){
    struct integration_params param ={temperature,mass_eff,mu_eff,degeneracy};
    gsl_integration_workspace * w= gsl_integration_workspace_alloc (700);
    double result, error;

    gsl_function F;
    F.function = &T_pressure_integrand;
    F.params = &param;
    gsl_integration_qagiu (&F, 0, 0, 1e-10, 700,
                    w, &result, &error);
    gsl_integration_workspace_free (w);
    return result- ((exp_a-1)/(exp_a))*4.0*M_PI*a/pow(mass_scalar,exp_a)*pow(scalar_density,2.0)+((exp_b-1)/(exp_b))*4.0*M_PI*b/pow(mass_vector,2.0)*pow(number_density,exp_b);
}
//the scalar density at finite T
double T_scalar_density(double degeneracy, double temperature, double mass_eff, double mu_eff){
    struct integration_params param ={temperature,mass_eff,mu_eff,degeneracy};
    gsl_integration_workspace * w= gsl_integration_workspace_alloc (700);
    double result, error;
    //cout<<"Tscalardenstiy mass eff "<<mass_eff<<" mueff"<<mu_eff<<" temperature"<<temperature<<"\n";
    gsl_function F;
    F.function = &T_scalar_density_integrand;
    F.params = &param;
    gsl_integration_qagiu (&F, 0, 0, 1e-10, 700,
                    w, &result, &error);
    gsl_integration_workspace_free (w);
    return result;
}
//the number density at finite T
double T_number_density(double degeneracy, double temperature, double mass_eff, double mu_eff){
     struct integration_params param ={temperature,mass_eff,mu_eff,degeneracy};
    gsl_integration_workspace * w= gsl_integration_workspace_alloc (700);
    double result, error;

    gsl_function F;
    F.function = &T_number_density_integrand;
    F.params = &param;
    gsl_integration_qagiu (&F, 0, 0, 1e-10, 700,
                    w, &result, &error);
    gsl_integration_workspace_free (w);
    return result;
}
