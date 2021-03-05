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
unsigned int terms [2]={1,1};


int main(int argc, char** argv) {
    //originial minimum
    double a=9.07314815896341*4*M_PI/pow(550,2.0);
    double b=13.897700851874*4*M_PI/pow(783,2.0);
    
    pair<double, double> q=get_crit( nucleon_mass, saturation_density,  binding_energy, degeneracy, {2.0}, {2.0}, false);
    cout<<q.first<<" "<<q.second/saturation_density<<"\n";
    double critical_temperature=q.first;
    double critical_density=q.second;
    
    double boundaries_model [4][3]={{5,15,1},{5,15,1},{1.7,2.3,0.1},{1.7,2.3,0.1}};
    double boundaries_crit [2][3]={{17.5,18.3,0.1},{(0.05/0.16)*saturation_density,(0.07/0.16)*saturation_density,(0.0025/0.16)*saturation_density}};
    struct interaction_params p={nucleon_mass, degeneracy, saturation_density,binding_energy, 17.5, (0.05/0.16)*saturation_density,terms };
    struct interaction_params p_org={nucleon_mass, degeneracy, saturation_density,binding_energy, critical_temperature, critical_density, terms};
    
    tuple<double, double,double,double, double, bool> test= get_interaction_4D(& p,false, {2.1}, {2.1}, {6*4*M_PI/pow(550,2.0)}, {6*4*M_PI/pow(783,2.0)});
    //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model  ,  boundaries_crit ,terms, "solutions_scan_test_v2", false,2);
   
    return 0;
}
