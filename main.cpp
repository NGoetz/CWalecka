#include <wplot.h>
#include <wsolvers.h>
#include <wfunctions.h>
#include <w2crit.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <tuple>
#include <time.h>
#include <filesystem>
#include <cstring>
#include <cassert>
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
  double boundaries_model_full[8][2]={{1,22},{1,22},{1,22},{1,22},{1e-11,150},{1e-11,150},{1e-11,150},{1e-11,150}};
  //double boundaries_model_III[8][2]={{1.8042023,1.8042025},{3.0631797,3.0631799},{6.6860892,6.6860894},{20.7276153,20.7276155},{9.224000e+01,9.224001e+01},{3.986262e+01,3.986264e+01},{1.066765e-01,1.066767e-01},{2.160278e-11,2.16028e-11}};
  //double boundaries_model_I[8][2]={{1.7614679,1.7614679},{3.8453863,3.8453863},{4.4772660,4.4772660},{6.7707861,6.7707861},{83.15987,83.15987},{61.44706,61.44706},{31.08395,31.08395},{0.3127069,0.3127069}};
  //double boundaries_model_II[8][2]={{1.8033076,1.8033078},{3.0693812,3.0693814},{7.9232547,7.9232549},{10.7986977,10.7986979},{9.204349e+01,9.204351e+01},{3.968765e+01,3.968767e+01},{1.306486e-01,1.306488e-01},{2.434033e-03,2.434035e-03}};
  //double boundaries_model_IV[8][2]={{1.768139,1.7681392},{3.5293514,3.5293516},{5.4352786,5.4352788},{6.3809822,6.3809824},{8.450947e+01,8.450949e+01},{3.843138e+01,3.84314e+01},{7.958556e+00,7.958558e+00},{1.552592e+00,1.552595e+00}};
  //double boundaries_model_IV_2[8][2]={{1.7681391,1.7681391},{3.5293515,3.5293515},{5.4352787,5.4352787},{6.3809823,6.3809823},{8.450948e+01,8.450948e+01},{3.843139e+01,3.843139e+01},{7.958557e+00,7.958557e+00},{1.552593e+00,1.552593e+00}};
  //double
  //crit II: spinodial_r=3.116
  double boundaries_crit_I[6][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{50,50,0.1},{(3)*saturation_density,(3)*saturation_density,1},{2.7*saturation_density,2.7*saturation_density,1},{3.12*saturation_density,3.32*saturation_density,0.1*saturation_density}};
  double boundaries_crit_II[6][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{50,50,0.1},{(3)*saturation_density,(3)*saturation_density,1},{2.85*saturation_density,2.85*saturation_density,1},{3.016027*saturation_density,3.216027*saturation_density,0.1*saturation_density}};
  double boundaries_crit_III[6][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{50,50,0.1},{(4)*saturation_density,(4)*saturation_density,1},{3.9*saturation_density,3.9*saturation_density,1},{3.982027*saturation_density,4.182028*saturation_density,0.1*saturation_density}};
  double boundaries_crit_IV[6][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{100,100,0.1},{(3)*saturation_density,(3)*saturation_density,1},{2.5*saturation_density,2.5*saturation_density,1},{3.21503*saturation_density,3.41503*saturation_density,0.1*saturation_density}};
  double boundaries_crit_V[6][3]={{18,18,0.1},{(0.06/0.16)*saturation_density,(0.06/0.16)*saturation_density,1},{100,100,0.1},{(4)*saturation_density,(4)*saturation_density,1},{3.6*saturation_density,3.6*saturation_density,1},{4.17703*saturation_density,4.37704*saturation_density,0.1*saturation_density}};

  struct plot_params p_ag_I={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_1,{},ag_coeff_1b,0.21};
  struct plot_params p_ag_II={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_2,{},ag_coeff_2b,0.21};
  //press_plot(3.11602,1e-9,3.11603,&p_ag_II);
  //T_press_plot(3,0.01,3.2,&p_ag_II,1);

  /*cout<<press_dn_solv(degeneracy, nucleon_mass,3.116028*saturation_density, {}, {}, ag_coeff_2b, ag_exp_2)<<endl;
  cout<<press_dn_solv(degeneracy, nucleon_mass,3.116027*saturation_density, {}, {}, ag_coeff_2b, ag_exp_2)<<endl;
  cout<<press_dn_solv(degeneracy, nucleon_mass,3.116026*saturation_density, {}, {}, ag_coeff_2b, ag_exp_2)<<endl;
  cout<<press_dn_solv(degeneracy, nucleon_mass,3.116025*saturation_density, {}, {}, ag_coeff_2b, ag_exp_2)<<endl;*/
  assert(argc>=4);
  array<double [6] [3],5>crit;
  memcpy(crit[0], boundaries_crit_I, sizeof(boundaries_crit_I));
  memcpy(crit[1], boundaries_crit_II, sizeof(boundaries_crit_II));
  memcpy(crit[2], boundaries_crit_III, sizeof(boundaries_crit_III));
  memcpy(crit[3], boundaries_crit_IV, sizeof(boundaries_crit_IV));
  memcpy(crit[4], boundaries_crit_V, sizeof(boundaries_crit_V));

  unsigned int terms [5][2] ={{0,4},{1,3},{2,2},{3,1},{4,0}};
  vector<string> mode={"","latin"};
  vector<string> name={"/experiment_1/model_I","/experiment_1/model_II","/experiment_1/model_III","/experiment_1/model_IV","/experiment_1/model_V"};
  string basepath=filesystem::current_path().string() +"/experiment_1";
  filesystem::path bpath=basepath;
  filesystem::create_directory(bpath);
  //1. Argument: crit 2. terms 3. latin or not
  interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy,  boundaries_model_full, crit[atoi(argv[1])] ,  terms[atoi(argv[2])] ,name[atoi(argv[1])]+"_scalar"+to_string(atoi(argv[2]))+"_"+mode[atoi(argv[3])], false, 0,400000,mode[atoi(argv[3])]);

  /*
  string str;
  if(argc>1){
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.7*saturation_density, 3.22*saturation_density, boundaries_model, boundaries_crit ,  termsag ,"model_I"+str.assign(argv[1]), false, 0,200000,str.assign(argv[1]));
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.85*saturation_density, 3.116025*saturation_density, boundaries_model, boundaries_crit ,  termsag ,"model_II"+str.assign(argv[1]), false, 0,200000,str.assign(argv[1]));
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 3.90*saturation_density, 4.082027*saturation_density, boundaries_model, boundaries_crit_III ,  termsag ,"model_III"+str.assign(argv[1]), false, 0,200000,str.assign(argv[1]));
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.5*saturation_density, 3.31503*saturation_density, boundaries_model, boundaries_crit_IV ,  termsag ,"model_IV"+str.assign(argv[1]), false, 0,200000,str.assign(argv[1]));
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 3.60*saturation_density, 4.27703*saturation_density, boundaries_model, boundaries_crit_V ,  termsag ,"model_V"+str.assign(argv[1]), false, 0,20000,str.assign(argv[1]));
  }else{
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.7*saturation_density, 3.22*saturation_density, boundaries_model, boundaries_crit ,  termsag ,"model_I", false, 0,200000,"");
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.85*saturation_density, 3.116025*saturation_density, boundaries_model, boundaries_crit ,  termsag ,"model_II", false, 0,200000,"");
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 3.90*saturation_density, 4.082027*saturation_density, boundaries_model, boundaries_crit_III ,  termsag ,"model_III", false, 0,200000,"");
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 2.5*saturation_density, 3.31503*saturation_density, boundaries_model, boundaries_crit_IV ,  termsag ,"model_IV", false, 0,200000,"");
    interaction_2crit_grid(nucleon_mass, binding_energy, saturation_density,  degeneracy, 3.60*saturation_density, 4.27703*saturation_density, boundaries_model, boundaries_crit_V ,  termsag ,"model_V", false, 0,200000,"");
  }
  */


  /*struct plot_params p_ag_III_0={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ng_exp_3,{},ng_coeff_3b,0.21};//4.082
  struct plot_params p_ag_III_1={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ng_exp_3,{},ng_coeff_3b,50};//4.082
  struct plot_params p_ag_III_2={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ng_exp_3,{},ng_coeff_3b,279};//4.082*/

  //struct plot_params p_ag_IV={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_4,{},ag_coeff_4b,0.21};//3.315
  //struct plot_params p_ag_V={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_5,{},ag_coeff_5b,0.21};//4.278
  //struct plot_params p_ag_VI={binding_energy,saturation_density, nucleon_mass, degeneracy,{},ag_exp_6,{},ag_coeff_6b,0.21};
  return 0;

}
