#include <wfunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <vector>

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
double scalar_density_dn(double degeneracy, double mass_eff, double number_density, vector<double> scalar_coeff, vector<double> scalar_exp){
    double scalar_density_val=scalar_density(degeneracy, mass_eff, number_density);
    double numerator=(8.0*mass_eff*pow(M_PI,2.0));
    double denom1=0;
    double denom2=0;
    double denom3=0;
     for(int i=0;i<scalar_coeff.size(); i++){
        denom1+=12.0*scalar_coeff[i]*number_density*(scalar_exp[i]-1.0)*pow(scalar_density_val, scalar_exp[i]-2.0)*pow(M_PI,2.0);
        denom2+=6.0*scalar_coeff[i]*(scalar_exp[i]-1.0)*pow(scalar_density_val, scalar_exp[i]-2.0)*degeneracy*pow(mass_eff,2.0)*fmom(degeneracy,number_density);
        denom3+=-3*scalar_coeff[i]*(scalar_exp[i]-1.0)*pow(scalar_density_val, scalar_exp[i]-2.0)*degeneracy*pow(mass_eff,2.0)*log(1.0 + 2*fmom(degeneracy,number_density)/(E_eff(fmom(degeneracy,number_density),mass_eff)-fmom(degeneracy,number_density)))*E_eff(fmom(degeneracy,number_density), mass_eff);
    }
    double denom4= (8.0*pow(M_PI,2.0))*E_eff(fmom(degeneracy,number_density), mass_eff);
    
    return numerator/(denom1+denom2+denom3+denom4);
}

//the effective mass from scalar density
double mass_eff_from_scalar_density(double nucleon_mass, vector<double> scalar_coeff, vector<double> scalar_exp, double scalar_density){
    double renorm=0;
    for(int i=0;i<scalar_coeff.size(); i++){
        renorm+=scalar_coeff[i]*pow(scalar_density,scalar_exp[i]-1.0);
       // cout<<scalar_coeff[i]<<" "<<scalar_exp[i]<<endl;
    }
    //cout<<"renorm "<<renorm<<endl;
    
    if (renorm>nucleon_mass){
        return 0.0;
    }
    return nucleon_mass-renorm;
}
// its derivative after n
double mass_eff_from_scalar_density_dn(double mass_eff,  double number_density,  vector<double> scalar_coeff, vector<double> scalar_exp, double degeneracy){
    double scalar_density_val=scalar_density(degeneracy, mass_eff, number_density);
    double returner=0;
    for(int i=0;i<scalar_coeff.size(); i++){
        returner-=(scalar_exp[i]-1.0)*scalar_coeff[i]*pow(scalar_density_val,scalar_exp[i]-2.0);
    }
    return returner*scalar_density_dn(degeneracy, mass_eff, number_density, scalar_coeff, scalar_exp);
}

//the energy density at T=0 divided by the density.
double eps_over_n(double degeneracy, double mass_eff, double number_density,vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density){
    double eps=(degeneracy/(16.0*pow(M_PI,2)))*((2.0)*E_eff(fmom(degeneracy,number_density), mass_eff)*(pow(fmom(degeneracy,number_density),3)) + (pow(mass_eff,2))*(E_eff(fmom(degeneracy,number_density),mass_eff)*fmom(degeneracy,number_density) + pow(mass_eff,2)*log(mass_eff/(E_eff(fmom(degeneracy,number_density), mass_eff) + fmom(degeneracy,number_density)))));
    double scalar=0;
    for(int i=0;i<scalar_coeff.size(); i++){
        scalar+=((scalar_exp[i]-1.0)/scalar_exp[i])*scalar_coeff[i]*pow(scalar_density,scalar_exp[i]);
    }
    double vec=0;
    for(int i=0;i<vec_coeff.size(); i++){
        vec+=((vec_exp[i]-1.0)/vec_exp[i])*vec_coeff[i]*pow(number_density,vec_exp[i]);
    }
    eps+=scalar+vec;
    return eps/number_density;
}

//the derivate of the energy density divided by density at T=0. 
double eps_over_n_dn(double degeneracy, double mass_eff, double number_density, vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density){
    
    double scalar_density_dn_val = scalar_density_dn(degeneracy, mass_eff, number_density,  scalar_coeff, scalar_exp);
    double mass_eff_dn_val = mass_eff_from_scalar_density_dn( mass_eff,  number_density, scalar_coeff,  scalar_exp, degeneracy);
   
    double scalar=0;
    for(int i=0;i<scalar_coeff.size(); i++){
        scalar+= 16.0*((scalar_exp[i]-1.0)/scalar_exp[i])*E_eff(fmom(degeneracy,number_density), mass_eff)*pow(M_PI,2.0)*(scalar_coeff[i])*pow(scalar_density, scalar_exp[i]-1.0)* (-scalar_density + scalar_exp[i]*number_density* scalar_density_dn_val);
    }
    double vec=0;
    for(int i=0;i<vec_coeff.size(); i++){
        vec+= 16.0*(pow(vec_exp[i]-1,2.0)/vec_exp[i])*pow(number_density,vec_exp[i])*pow(M_PI,2.0)*(vec_coeff[i])* E_eff(fmom(degeneracy,number_density), mass_eff);
    }
    return ((1.0/(16.0*pow(number_density*M_PI,2.0)*E_eff(fmom(degeneracy,number_density), mass_eff) ))*(-2.0*number_density*pow(M_PI*mass_eff,2.0) - degeneracy*fmom(degeneracy,number_density)*pow(mass_eff,4.0)+ 4.0*number_density*pow(fmom(degeneracy,number_density)*M_PI,2.0) + 
   24.0*pow(number_density*M_PI,2.0)*mass_eff*mass_eff_dn_val + 4*degeneracy*fmom(degeneracy,number_density)*number_density*pow(mass_eff,3.0)*mass_eff_dn_val  - degeneracy*log(mass_eff/(fmom(degeneracy,number_density) + E_eff(fmom(degeneracy,number_density), mass_eff)))*pow(mass_eff,3.0)* E_eff(fmom(degeneracy,number_density), mass_eff)* (mass_eff - 4*number_density*mass_eff_dn_val) +scalar+vec));
}
//the pressure at T=0. 
double press(double degeneracy, double mass_eff, double number_density,vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density){
    double press=(((degeneracy/(16.0*pow(M_PI,2)))*((2.0/3.0)*E_eff(fmom(degeneracy,number_density),mass_eff)*(pow(fmom(degeneracy,number_density),3)) - (pow(mass_eff,2))*(E_eff(fmom(degeneracy,number_density),mass_eff)*fmom(degeneracy,number_density) + pow(mass_eff,2)*log(mass_eff/(E_eff(fmom(degeneracy,number_density),mass_eff) + fmom(degeneracy,number_density)))))));
    double scalar=0;
    for(int i=0;i<scalar_coeff.size(); i++){
        scalar-=((scalar_exp[i]-1.0)/scalar_exp[i])*scalar_coeff[i]*pow(scalar_density,scalar_exp[i]);
    }
    double vec=0;
    for(int i=0;i<vec_coeff.size(); i++){
        vec+=((vec_exp[i]-1.0)/vec_exp[i])*vec_coeff[i]*pow(number_density,vec_exp[i]);
    }
    return press+scalar+vec;
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
double muB_from_number_density(double degeneracy, double mass_eff, double number_density, vector<double> vec_coeff, vector<double> vec_exp){
    double vec=0;
    for(int i=0;i<vec_coeff.size(); i++){
        vec+=vec_coeff[i]*pow(number_density,vec_exp[i]-1);
    }
    return sqrt(pow(mass_eff,2.0)+pow(fmom(degeneracy,number_density),2.0))+vec;
}
//the effective chemical potential at finite temperature (at zero, it is just the energy)
double T_mu_eff(double number_density, double muB, vector<double> vec_coeff, vector<double> vec_exp){
    double vec=0;
    for(int i=0;i<vec_coeff.size(); i++){
        vec+=vec_coeff[i]*pow(number_density,vec_exp[i]-1);
    }
    return muB-vec;
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
double T_pressure(double degeneracy, double temperature, double mass_eff, double mu_eff, vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density, double number_density){
    struct integration_params param ={temperature,mass_eff,mu_eff,degeneracy};
    gsl_integration_workspace * w= gsl_integration_workspace_alloc (700);
    double result, error;

    gsl_function F;
    F.function = &T_pressure_integrand;
    F.params = &param;
    gsl_integration_qagiu (&F, 0, 0, 1e-10, 700,
                    w, &result, &error);
    gsl_integration_workspace_free (w);
    double scalar=0;
    for(int i=0;i<scalar_coeff.size(); i++){
        scalar-=((scalar_exp[i]-1.0)/scalar_exp[i])*scalar_coeff[i]*pow(scalar_density,scalar_exp[i]);
    }
    double vec=0;
    for(int i=0;i<vec_coeff.size(); i++){
        vec+=((vec_exp[i]-1.0)/vec_exp[i])*vec_coeff[i]*pow(number_density,vec_exp[i]);
    }
    return result+scalar+vec;
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
