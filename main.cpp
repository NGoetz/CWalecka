#include <wplot.h>
#include <wsolvers.h>
#include <wfunctions.h>
#include <w2crit.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <tuple>
#include <time.h>
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

int main(int argc, char* argv[]) {

  vector<double> ag_coeff_1={-83.15987,61.44706,-31.08395,0.3127069};
  vector<double> ag_exp_1={1.7614679,3.8453863,4.4772660,6.7707861};
  vector<double> ag_coeff_1b;
  for (int i=0;i<ag_coeff_1.size();i++){
    ag_coeff_1b.push_back(ag_coeff_1[i]/(pow(saturation_density,ag_exp_1[i]-1)));
  }
  vector<double> ag_coeff_2={-9.204350e+01,3.968766e+01,-1.306487e-01,2.434034e-03};
  vector<double> ag_exp_2={1.8033077,3.0693813,7.9232548,10.7986978};
  vector<double> ag_coeff_3={-9.224000e+01,3.986263e+01,-1.066766e-01,2.160279e-11};
  vector<double> ag_exp_3={1.8042024,3.0631798,6.6860893,20.7276154};

  vector<double> ag_coeff_3b;
  vector<double> ag_coeff_2b;
  for (int i=0;i<ag_coeff_3.size();i++){
    ag_coeff_3b.push_back(ag_coeff_3[i]/(pow(saturation_density,ag_exp_3[i]-1)));
  }
  for (int i=0;i<ag_coeff_2.size();i++){
    ag_coeff_2b.push_back(ag_coeff_2[i]/(pow(saturation_density,ag_exp_2[i]-1)));
  }
  vector<double> ag_coeff_4={-8.450948e+01,3.843139e+01,-7.958557e+00,1.552593e+00};
  vector<double> ag_exp_4={1.7681391,3.5293515,5.4352787,6.3809823};

  vector<double> ag_coeff_4b;
  for (int i=0;i<ag_coeff_4.size();i++){
    ag_coeff_4b.push_back(ag_coeff_4[i]/(pow(saturation_density,ag_exp_4[i]-1)));
  }

  vector<double> ag_coeff_5={-8.627959e+01,4.786488e+01,-1.406946e+01,1.182795e-04};
  vector<double> ag_exp_5={1.7782362,3.4936863,4.2528897,10.3240297};

  vector<double> ag_coeff_5b;
  for (int i=0;i<ag_coeff_5.size();i++){
    ag_coeff_5b.push_back(ag_coeff_5[i]/(pow(saturation_density,ag_exp_5[i]-1)));
  }
  vector<double> ag_coeff_6={-9.101665e+01,3.899891e+01,-4.856681e-01,1.935808e-02};
  vector<double> ag_exp_6={1.7989835,3.1098389,6.3017683,8.0937872};

  vector<double> ag_coeff_6b;
  for (int i=0;i<ag_coeff_6.size();i++){
    ag_coeff_6b.push_back(ag_coeff_6[i]/(pow(saturation_density,ag_exp_6[i]-1)));
  }

  double boundaries_model[8][2]={{1,2},{2.5,4},{4,8},{6,22},{50,150},{30,90},{0.1,50},{1e-11,10}};
  double boundaries_model_III[8][2]={{1.8042023,1.8042025},{3.0631797,3.0631799},{6.6860892,6.6860894},{20.7276153,20.7276155},{9.224000e+01,9.224001e+01},{3.986262e+01,3.986264e+01},{1.066765e-01,1.066767e-01},{2.160278e-11,2.16028e-11}};
  //double boundaries_model_I[8][2]={{1.7614679,1.7614679},{3.8453863,3.8453863},{4.4772660,4.4772660},{6.7707861,6.7707861},{83.15987,83.15987},{61.44706,61.44706},{31.08395,31.08395},{0.3127069,0.3127069}};
  //double boundaries_model_II[8][2]={{1.8033076,1.8033078},{3.0693812,3.0693814},{7.9232547,7.9232549},{10.7986977,10.7986979},{9.204349e+01,9.204351e+01},{3.968765e+01,3.968767e+01},{1.306486e-01,1.306488e-01},{2.434033e-03,2.434035e-03}};
  //double boundaries_model_IV[8][2]={{1.768139,1.7681392},{3.5293514,3.5293516},{5.4352786,5.4352788},{6.3809822,6.3809824},{8.450947e+01,8.450949e+01},{3.843138e+01,3.84314e+01},{7.958556e+00,7.958558e+00},{1.552592e+00,1.552595e+00}};
  //double boundaries_model_IV_2[8][2]={{1.7681391,1.7681391},{3.5293515,3.5293515},{5.4352787,5.4352787},{6.3809823,6.3809823},{8.450948e+01,8.450948e+01},{3.843139e+01,3.843139e+01},{7.958557e+00,7.958557e+00},{1.552593e+00,1.552593e+00}};
  //double
  //crit II: spinodial_r=3.116
  double boundaries_crit[4][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{50,50,0.1},{(3)*saturation_density,(3)*saturation_density,1}};
  double boundaries_crit_III[4][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{50,50,0.1},{(4)*saturation_density,(4)*saturation_density,1}};
  double boundaries_crit_IV[4][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{100,100,0.1},{(3)*saturation_density,(3)*saturation_density,1}};
  double boundaries_crit_V[4][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{100,100,0.1},{(4)*saturation_density,(4)*saturation_density,1}};
  double boundaries_crit_VI[4][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{125,125,0.1},{(4)*saturation_density,(4)*saturation_density,1}};
  //struct interaction_params_2crit p={nucleon_mass, degeneracy, saturation_density,binding_energy, 18, (0.06/0.16)*saturation_density,50, 3*saturation_density,2.7*saturation_density,3.22*saturation_density,  termsag };
  //interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.7*saturation_density, 3.22*saturation_density, boundaries_model, boundaries_crit ,  termsag ,"ag_test_sol_firstAg", false, 0,8000);
  //interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.85*saturation_density, 3.116*saturation_density, boundaries_model, boundaries_crit ,  termsag ,"ag_test_sol_II", false, 0,200000);
  time_t timer;
  time(&timer);

  string str;
  if(argc>1){
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 3.90*saturation_density, 4.082*saturation_density, boundaries_model_III, boundaries_crit_III ,  termsag ,"run_test_III"+str.assign(argv[1]), false, 0,800,str.assign(argv[1]));
  }else{
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 3.90*saturation_density, 4.082*saturation_density, boundaries_model_III, boundaries_crit_III ,  termsag ,"run_test_III", false, 0,800,"");
  }

  time_t timer2;
  time(&timer2);
  cout<<difftime(timer, timer2)<<endl;

  /*struct plot_params p_ag_III_0={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ng_exp_3,{},ng_coeff_3b,0.21};//4.082
  struct plot_params p_ag_III_1={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ng_exp_3,{},ng_coeff_3b,50};//4.082
  struct plot_params p_ag_III_2={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ng_exp_3,{},ng_coeff_3b,279};//4.082*/

  //struct plot_params p_ag_IV={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_4,{},ag_coeff_4b,0.21};//3.315
  //struct plot_params p_ag_V={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_5,{},ag_coeff_5b,0.21};//4.278
  //struct plot_params p_ag_VI={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_6,{},ag_coeff_6b,0.21};
  return 0;
  
}
