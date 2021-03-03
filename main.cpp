#include <wplot.h>
#include <wsolvers.h>
#include <wfunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <tuple>
using namespace std;
//the value of the degeneracy
double const degeneracy=4.0;
//the nucleon mass
double const nucleon_mass=938;
//the mass of the scalar particle
double const mass_scalar=550;
//the mass of the vector particle
double const mass_vector=783;
//the conversion factor from fm to MeV
double const conv=197.3;
//position of the energy density minimum in fm -3
double const minpos=0.16;
// in MeV
double const saturation_density = minpos*pow(conv, 3);
//value of the minimum energy density minus the nucleus mass
double const binding_energy=16.3;
//exponent of first interaction term
double const exp_a=2.0;
//exponent of second interaction term
double const exp_b=2.0;

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

int main(int argc, char** argv) {
    pair<double, double>x=get_coeff(nucleon_mass,  saturation_density,  binding_energy,  degeneracy,  mass_scalar,  mass_vector,  exp_a,  exp_b, false);
    cout << x.first << " "<<x.second<<"\n";
    pair<double, double> q=get_crit( nucleon_mass, saturation_density,  binding_energy,  mass_scalar,  mass_vector,  degeneracy, 2.0, 2.0, false);
    cout<<q.first<<" "<<q.second/saturation_density<<"\n";
    double critical_temperature=q.first;
    double critical_density=q.second;
    double x0=x.first;
    double x1=x.second;
    double x2=2.;
    double x3=2.;
    pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,x0, x2, saturation_density,mass_scalar, false);
    double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, x0, x2, x1, x3, output.second, mass_vector, mass_scalar);
    double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,x0, x2,x1, x3, output.second, mass_vector, mass_scalar);
    double y2 =T_press_dn_solv(degeneracy, nucleon_mass, critical_density, mass_vector, mass_scalar, x0, x2, x1, x3, critical_temperature);
    double y3 =T_press_dn2_solv(degeneracy, nucleon_mass, critical_density, mass_vector, mass_scalar, x0, x2,x1, x3, critical_temperature);
    cout<<"y0 "<<y0<<" y1 "<<y1<<" y2 "<<y2<<" y3 "<<y3<<"\n";
    x0=3.98235;
    x1=0.0891071;
    x2=2.00645;
    x3=2.28271;
    output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,x0, x2, saturation_density,mass_scalar, false);
    y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, x0, x2, x1, x3, output.second, mass_vector, mass_scalar);
    y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,x0, x2,x1, x3, output.second, mass_vector, mass_scalar);
    y2 =T_press_dn_solv(degeneracy, nucleon_mass, critical_density, mass_vector, mass_scalar, x0, x2, x1, x3, critical_temperature);
    y3 =T_press_dn2_solv(degeneracy, nucleon_mass, critical_density, mass_vector, mass_scalar, x0, x2,x1, x3, critical_temperature);
    cout<<"y0 "<<y0<<" y1 "<<y1<<" y2 "<<y2<<" y3 "<<y3<<"\n";
    cout<<"Test "<<"\n";
    double boundaries_model [4][3]={{5,15,1},{5,15,1},{1.7,2.3,0.1},{1.7,2.3,0.1}};
    double boundaries_crit [2][3]={{17.5,18.3,0.1},{(0.05/0.16)*saturation_density,(0.07/0.16)*saturation_density,(0.0025/0.16)*saturation_density}};
    struct interaction_params p={nucleon_mass, degeneracy, saturation_density,binding_energy, 17.5, (0.05/0.16)*saturation_density, mass_scalar, mass_vector };
    tuple<double, double,double,double, double, bool> test= get_interaction_4D(& p,true, 2.1, 2.1, 6, 6);
    interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  mass_scalar,  mass_vector,  boundaries_model  ,  boundaries_crit , "solutions_scan_test_v2", false,2);
    
    return 0;
}

