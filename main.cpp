#include <../../src/wplot.h>
#include <../../src/wsolvers.h>
#include <../../src/wfunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
//the value of the degeneracy
double const degeneracy=4.0;
//the nucleon mass
double const nucleon_mass=939;
//the mass of the scalar particle
double const mass_scalar=550;
//the mass of the vector particle
double const mass_vector=783;
//the convergence factor from fm to MeV
double const conv=197.3;
//positon of the energy density minimum in fm -3
double const minpos=0.153;
// in MeV
double const saturation_density = minpos*pow(conv, 3);
//value of the minimum energy density minus the nucleus mass
double const binding_energy=16.3;

int main() {
    plot_params p={binding_energy,saturation_density,nucleon_mass, degeneracy,mass_scalar,mass_vector,1};
    //pair<double, double> x1= get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, m_scalar, m_vec, false);
    //pair<double, double> x2=get_mass_eff_scalar_density(degeneracy, nucleon_mass,x1.first,saturation_density,m_scalar, false);
    //cout << x1.first << " "<<x1.second<<"\n";
   // pair<double, double> x3=get_mass_eff_scalar_density(degeneracy, nucleon_mass,9.537,saturation_density,m_scalar, true);
    //cout << x3.first << " "<<x3.second<<"\n";
   // cout <<incompress(degeneracy, x2.first, saturation_density,x1.first,x1.second,x2.second, m_vec,m_scalar)<<"\n";
    //energy_pp_plot(0.1,0.01,1.3,&p);
    //T_energy_pp_plot(0.1,0.01,1.3,&p);
    //T_press_plot(0.1,0.05,1.3,&p);
    pair<double, double> result=get_crit(nucleon_mass,saturation_density, binding_energy, mass_scalar, mass_vector, degeneracy,true);
    cout<<result.first<<" "<<result.second<<"\n";
    return 0;
}
