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
unsigned int terms11 [2]={1,1};
unsigned int terms20 [2]={2,0};
unsigned int terms02 [2]={0,2};
unsigned int termsag[2]={0,4};

int main(int argc, char** argv) {
  //originial minimum
  double a=9.07314815896341*4*M_PI/pow(550,2.0);
  double b=13.897700851874*4*M_PI/pow(783,2.0);

  pair<double, double> q=get_crit( nucleon_mass, saturation_density,  binding_energy, degeneracy, {2.0}, {2.0}, false);
  // cout<<q.first<<" "<<q.second/saturation_density<<"\n";
  //cout<<"NEW CRIT"<<endl;
  //q=get_crit( nucleon_mass, saturation_density,  binding_energy, degeneracy, {2.0}, {2.0}, true,19.09,540000,{a},{b});
  ////cout<<q.first<<" "<<q.second/saturation_density<<"\n";
  double critical_temperature=q.first;
  double critical_density=q.second;
  //double boundaries_model_test[4][3]={{-19*4*M_PI/pow(550,2.0),20*4*M_PI/pow(550,2.0),4*4*M_PI/pow(550,2.0)},{-19*4*M_PI/pow(783,2.0),20*4*M_PI/pow(783,2.0),4*4*M_PI/pow(783,2.0)},{1,3,0.2},{1,3,0.2}};
  //double boundaries_model_test2[4][3]={{-19*4*M_PI/pow(550,2.0),20*4*M_PI/pow(550,2.0),2*4*M_PI/pow(550,2.0)},{-19*4*M_PI/pow(783,2.0),20*4*M_PI/pow(783,2.0),2*4*M_PI/pow(783,2.0)},{1,4,0.1},{1,4,0.1}};


  //double boundaries_model[4][3]={{2*4*M_PI/pow(550,2.0),20*4*M_PI/pow(550,2.0),1*4*M_PI/pow(550,2.0)},{2*4*M_PI/pow(783,2.0),20*4*M_PI/pow(783,2.0),1*4*M_PI/pow(783,2.0)},{1,3,0.1},{1,3,0.1}};
  //double boundaries_model [4][3]={{9*4*M_PI/pow(550,2.0),16*4*M_PI/pow(550,2.0),1*4*M_PI/pow(550,2.0)},{3*4*M_PI/pow(783,2.0),19*4*M_PI/pow(783,2.0),2*4*M_PI/pow(783,2.0)},{1.5,2.5,0.1},{1.5,2.5,0.1}};
  //double boundaries_crit_ex [2][3]={{17.5,18.3,0.1},{(0.05/0.16)*saturation_density,(0.07/0.16)*saturation_density,(0.0025/0.16)*saturation_density}};
  //double boundaries_crit_test [2][3]={{17.5,17.5,0.1},{(0.05/0.16)*saturation_density,(0.05/0.16)*saturation_density,(0.0025/0.16)*saturation_density}};
  //double boundaries_crit_wal [2][3]={{critical_temperature,critical_temperature,0.1},{critical_density,critical_density,(0.0025/0.16)*saturation_density}};
  //struct interaction_params p={nucleon_mass, degeneracy, saturation_density,binding_energy, 17.5, (0.05/0.16)*saturation_density,terms };
  //struct interaction_params p_org={nucleon_mass, degeneracy, saturation_density,binding_energy, critical_temperature, critical_density, terms};
  //vector<double> ag_coeff_1={pow(0.16,2.0)*pow(197.3,3.0)*-83.15987,pow(0.16,2.0)*pow(197.3,3.0)*61.44706,pow(0.16,2.0)*pow(197.3,3.0)*-31.08395,pow(0.16,2.0)*pow(197.3,3.0)*0.3127069};
  vector<double> ag_coeff_1={-83.15987,61.44706,-31.08395,0.3127069};
  vector<double> ag_exp_1={1.7614679,3.8453863,4.4772660,6.7707861};
  for (int i=0;i<ag_coeff_1.size();i++){
    ag_coeff_1[i]=ag_coeff_1[i]/(pow(saturation_density,ag_exp_1[i]-1));
    cout<<ag_coeff_1[i]<<" ";
  }
  cout<<endl;
  double boundaries_model[8][2]={{1,7},{1,7},{1,7},{1,7},{-0.002,0.002},{-0.002,0.002},{-0.002,0.002},{-0.002,0.002}};
  double boundaries_crit[4][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{50,50,0.1},{(3)*saturation_density,(3)*saturation_density,1}};
  interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.7*saturation_density, 3.22*saturation_density, boundaries_model, boundaries_crit ,  termsag ,"ag_test", false, 0,1600);
  //struct plot_params p_sc={binding_energy,saturation_density,nucleon_mass,degeneracy,{2,2.5},{},{-1e-5,1e-7},{},0};
  //struct plot_params p_vec={binding_energy,saturation_density,nucleon_mass,degeneracy,{},{2,2.9},{},{1e-4,-1e-9},0};
  //struct plot_params p_w={binding_energy,saturation_density,nucleon_mass,degeneracy,{2},{2},{a},{b},0};
  /*
  struct plot_params p_w_ag2={binding_energy,saturation_density,nucleon_mass,degeneracy,{2},{2},{74.1048/(0.16*pow(saturation_density,1.0))},{ 56.0021/(0.16*pow(saturation_density,1.0))},0};//convention from 17.02.21
  struct plot_params p_ag1_1={binding_energy,saturation_density, nucleon_mass, degeneracy,{},{1.7614679,3.8453863,4.4772660,6.7707861},{},ag_coeff_1,0.21};
  struct plot_params p_ag1_2={binding_energy,saturation_density, nucleon_mass, degeneracy,{},{1.7614679,3.8453863,4.4772660,6.7707861},{},ag_coeff_1,18};
  struct plot_params p_ag1_3={binding_energy,saturation_density, nucleon_mass, degeneracy,{},{1.7614679,3.8453863,4.4772660,6.7707861},{},ag_coeff_1,50};
  struct plot_params ps [3]={p_ag1_1,p_ag1_2,p_ag1_3};



  struct plot_params s1v1_sol1={binding_energy,saturation_density, nucleon_mass, degeneracy,{2},{2},{0.000376914183565},{0.00028485858464},q.first};
  struct plot_params s1v1_sol2={binding_energy,saturation_density, nucleon_mass, degeneracy,{1.93514001415468},{1.91667693248702},{0.000962171732948},{0.000945866431751},q.first};
  struct plot_params s0v2_sol={binding_energy,saturation_density, nucleon_mass, degeneracy,{},{1.71268642792956,4.27860033978592},{},{-0.003376104298288,2.3161985281533E-19},q.first};
  struct plot_params pw [2]={s1v1_sol1,s0v2_sol};

  pair<double, double> s=get_coeff( nucleon_mass,  saturation_density,  binding_energy,  degeneracy, {1.5,1.6}, {}, false,terms20);
  cout<<s.first<<" "<<s.second<<endl;
  struct plot_params s2v0_test1={binding_energy,saturation_density, nucleon_mass, degeneracy,{1.5,1.6},{},{s.first, s.second},{},0.21};
  struct plot_params s2v0_test2={binding_energy,saturation_density, nucleon_mass, degeneracy,{1.5,1.6},{},{s.first, s.second},{},5};
  struct plot_params s2v0_test3={binding_energy,saturation_density, nucleon_mass, degeneracy,{1.5,1.6},{},{s.first, s.second},{},10};
  struct plot_params p0 [3]={s2v0_test1,s2v0_test2,s2v0_test3};*/
  //q=get_crit( nucleon_mass, saturation_density,  binding_energy, degeneracy, {2.0}, {2.0}, true,19.09,540000,{s.first},{s.second});
  //cout<<q.first<<" "<<q.second<<endl;
  //
  /*
  struct plot_params s1v1_comp={binding_energy,saturation_density, nucleon_mass, degeneracy,{1.89977208377173},{2.7714313971606},{0.000387261300706},{0.000000000921380362759002},17.5};
  struct plot_params s2v0_comp={binding_energy,saturation_density, nucleon_mass, degeneracy,{1.89715562976373  , 2.70220323171616},{},{0.000422148625599,-0.00000000293734967138934},{},17.5};
  struct plot_params s0v2_comp={binding_energy,saturation_density, nucleon_mass, degeneracy,{},{1.83643855767161,2.79768468004429},{},{-0.000834587370244,0.00000000057654923114384},17.5};
  struct plot_params pex [3]={s1v1_comp, s0v2_comp, s2v0_comp};*/
  //energy_pp_plot(0.01, 0.01, 2, pex, 3);
  //eff_mass_plot_p(0.01, 0.01, 4, p0, 1);
  //press_plot(0.01,0.01,1,&s2v0_test1);
  //scalar_density_plot(0.01, 0.01, 4, p0);
  //T_press_plot(0.1, 0.01,1,pex, 3);
  //eff_mass_plot_p(0.01, 0.01, 2, pw, 2);
  //T_press_plot(0.25, 0.02, 0.6,pw, 2);
  //energy_pp_dn_plot(0.05, 0.01, 0.4, &p_sc);
  //tuple<double, double,double,double, double, bool> test= get_interaction_4D(& p,false, {2.1}, {2.1}, {6*4*M_PI/pow(550,2.0)}, {6*4*M_PI/pow(783,2.0)});
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_w  ,  boundaries_crit_wal ,terms, "solution_scan_walecka_crit_all_v2s0", false,0);
  //cout<<"FIRST DONE"<<endl;
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_test2  ,  boundaries_crit_test ,terms02, "new_w_s0v2_test_low_n", false,0);
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_test  ,  boundaries_crit_wal ,terms11, "w_s1v1_all", false,0);
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_test  ,  boundaries_crit_wal ,terms20, "w_s2v0_all", false,0);
  // interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_test  ,  boundaries_crit_wal ,terms02, "w_s0v2_all", false,0);
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_test  ,  boundaries_crit_ex ,terms11, "ex_s1v1_1", false,1);
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_test  ,  boundaries_crit_ex ,terms20, "ex_s2v0_1", false,1);
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model_test  ,  boundaries_crit_ex ,terms02, "ex_s0v2_1", false,1);
  //interaction_4D_crit_grid( nucleon_mass,  binding_energy,  saturation_density,  degeneracy,  boundaries_model ,  boundaries_crit_ex ,terms, "solution_scan_crit_ex_range", false,1);
  return 0;
  //find for CP 17.5/0.34375 (probably needs wider range)
}
