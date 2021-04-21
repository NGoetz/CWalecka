#include <wfunctions.h>
#include <wsolvers.h>
#include <w2crit.h>
#include <wplot.h>
#include <stdio.h>
#include <stdlib.h>
#include <filesystem>
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
#include <string.h>
#include <omp.h>
#include <random>

//algorithm to find a EOS with a liquid gas and hadron-QGP critical point

//class to handle CSV input
class CSVRow
{
    public:
        std::string_view operator[](std::size_t index) const
        {
            return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);
            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(',', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

//define equations to be solved for 2CP condition
int interaction_root_2crit(const gsl_vector * x, void *params, gsl_vector * f){
  double nucleon_mass=((struct interaction_params_2crit *) params)->nucleon_mass;
  double degeneracy=((struct interaction_params_2crit *) params)->degeneracy;
  double saturation_density=((struct interaction_params_2crit *) params)->saturation_density;
  double binding_energy=((struct interaction_params_2crit *) params)->binding_energy;
  double critical_temperature_lg=((struct interaction_params_2crit *) params)->critical_temperature_lg;
  double critical_density_lg=((struct interaction_params_2crit *) params)->critical_density_lg;
  double critical_temperature_qgp=((struct interaction_params_2crit *) params)->critical_temperature_qgp;
  double critical_density_qgp=((struct interaction_params_2crit *) params)->critical_density_qgp;
  double spinodial_r_density=((struct interaction_params_2crit *) params)->spinodial_r_density;
  double spinodial_l_density=((struct interaction_params_2crit *) params)->spinodial_l_density;
  unsigned int * terms =((struct interaction_params_2crit *) params)->terms;

  vector<double> scalar_coeff, scalar_exp, vector_coeff, vector_exp;

  for(unsigned int i=0; i<terms[0];i++){
    scalar_exp.push_back(gsl_vector_get(x,i));
    scalar_coeff.push_back(gsl_vector_get(x,terms[0]+terms[1]+i));
  }
  for(unsigned int i=terms[0]; i<terms[0]+terms[1];i++){
    vector_exp.push_back(gsl_vector_get(x,i));
    vector_coeff.push_back(gsl_vector_get(x,terms[0]+terms[1]+i));
  }

  pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,scalar_coeff, scalar_exp, saturation_density, false);
  //predict binding at demanded saturation density
  const double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, scalar_coeff, scalar_exp,vector_coeff, vector_exp, output.second);
  const double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,scalar_coeff, scalar_exp,vector_coeff, vector_exp, output.second);
  //there is an inflection point at the critical temperature/density
  pair<double, double> pressures_lg=T_press_dn12_solv( degeneracy,  nucleon_mass, critical_density_lg,  scalar_coeff, scalar_exp,vector_coeff, vector_exp, critical_temperature_lg);
  const double y2 =pressures_lg.first;
  const double y3 =pressures_lg.second;
  pair<double, double> pressures_qgp=T_press_dn12_solv( degeneracy,  nucleon_mass, critical_density_qgp,  scalar_coeff, scalar_exp,vector_coeff, vector_exp, critical_temperature_qgp);
  const double y4 =pressures_qgp.first;
  const double y5 =pressures_qgp.second;
  //position of spinodial boundaries
  const double y6=press_dn_solv( degeneracy, nucleon_mass,spinodial_l_density, scalar_coeff, scalar_exp, vector_coeff,vector_exp);
  const double y7=press_dn_solv( degeneracy, nucleon_mass,spinodial_r_density, scalar_coeff, scalar_exp, vector_coeff,vector_exp);
  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);
  gsl_vector_set (f, 3, y3);
  gsl_vector_set (f, 4, y4);
  gsl_vector_set (f, 5, y5);
  gsl_vector_set (f, 6, y6);
  gsl_vector_set (f, 7, y7);
  return GSL_SUCCESS;
}

int print_state_interaction_2crit (size_t iter, gsl_multiroot_fsolver * s){
  printf ("iter = %3lu x = % .3f % .3f % .3f % .3f % .3e % .3e % .3e % .3e "
  "f(x) = % .3e % .3e % .3e % .3e  % .3e % .3e % .3e % .3e\n",
  iter,
  gsl_vector_get (s->x, 0),
  gsl_vector_get (s->x, 1),
  gsl_vector_get (s->x, 2),
  gsl_vector_get (s->x, 3),
  gsl_vector_get (s->x, 4),
  gsl_vector_get (s->x, 5),
  gsl_vector_get (s->x, 6),
  gsl_vector_get (s->x, 7),
  gsl_vector_get (s->f, 0),
  gsl_vector_get (s->f, 1),
  gsl_vector_get (s->f, 2),
  gsl_vector_get (s->f, 3),
  gsl_vector_get (s->f, 4),
  gsl_vector_get (s->f, 5),
  gsl_vector_get (s->f, 6),
  gsl_vector_get (s->f, 7));
  return 0;
}
//find the coefficients and exponents of 4 terms with the nuclear saturation point and the nuclear liquid/gas transition
tuple<double, vector<double>, vector<double>, bool> get_interaction_2crit(void * p, bool print, vector<double> init_exp,  vector<double> init_coeff) {
  gsl_set_error_handler(&gsl_handler);
  int status;
  size_t iter = 0;
  const size_t num = 8;
  gsl_vector *x = gsl_vector_alloc (num);
  const gsl_multiroot_fsolver_type * K = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (K, num);

  gsl_multiroot_function f = {&interaction_root_2crit, num, p};
  unsigned int * terms=((struct interaction_params_2crit *) p)->terms;
  assert(terms[0]+terms[1]==4);

  gsl_vector_set (x, 0, init_exp[0]);
  gsl_vector_set (x, 1, init_exp[1]);
  gsl_vector_set (x, 2, init_exp[2]);
  gsl_vector_set (x, 3, init_exp[3]);
  gsl_vector_set (x, 4, init_coeff[0]);
  gsl_vector_set (x, 5, init_coeff[1]);
  gsl_vector_set (x, 6, init_coeff[2]);
  gsl_vector_set (x, 7, init_coeff[3]);

  gsl_multiroot_fsolver_set (s, &f, x);

  if (print) {
    print_state_interaction_2crit (iter, s);
  }
  // perform the root solving
  do{
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    if (print) {
      print_state_interaction_2crit (iter, s);
    }

    if (status)
    break;

    status =
    gsl_multiroot_test_residual (s->f, 1e-11);

  }while (status == GSL_CONTINUE && iter < 1000);
  if (print){
    printf ("status = %s\n", gsl_strerror (status));

  }
  double error=abs(gsl_vector_get(s->f,0))+abs(gsl_vector_get(s->f,1))+abs(gsl_vector_get(s->f,2))+abs(gsl_vector_get(s->f,3))+abs(gsl_vector_get(s->f,4))+abs(gsl_vector_get(s->f,5))+abs(gsl_vector_get(s->f,6))+abs(gsl_vector_get(s->f,7));
  tuple<double, vector<double>, vector<double>, bool> result={error,{gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),gsl_vector_get(s->x,2),gsl_vector_get(s->x,3)},{gsl_vector_get(s->x,4),gsl_vector_get(s->x,5),gsl_vector_get(s->x,6),gsl_vector_get(s->x,7)},status==GSL_SUCCESS};
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return result;
}


//validate interaction terms found for 2 CPs
bool validate_interaction_2crit(double error, vector<double> exp_guess, vector<double> coeff_guess,vector<double> exp, vector<double> coeff, string filename, int timestamp, void* params, double criterium){
  double eps=numeric_limits<double>::epsilon();
  double nucleon_mass=((struct interaction_params_2crit *) params)->nucleon_mass;
  double degeneracy=((struct interaction_params_2crit *) params)->degeneracy;
  double saturation_density=((struct interaction_params_2crit *) params)->saturation_density;
  double binding_energy=((struct interaction_params_2crit *) params)->binding_energy;
  double critical_temperature_lg=((struct interaction_params_2crit *) params)->critical_temperature_lg;
  double critical_density_lg=((struct interaction_params_2crit *) params)->critical_density_lg;
  double critical_temperature_qgp=((struct interaction_params_2crit *) params)->critical_temperature_qgp;
  double critical_density_qgp=((struct interaction_params_2crit *) params)->critical_density_qgp;
  double spinodial_l_density=((struct interaction_params_2crit *) params)->spinodial_l_density;
  double spinodial_r_density=((struct interaction_params_2crit *) params)->spinodial_r_density;
  unsigned int * terms =((struct interaction_params_2crit *) params)->terms;

  vector<double> scalar_coeff, scalar_exp, vec_coeff, vec_exp,scalar_coeff_guess, scalar_exp_guess, vec_coeff_guess, vec_exp_guess;
  //validated flags wether the conditions are reproduced well enough.
  //if this is not the case, the set is not saved in the solution but in the NV file
  bool validated=true;
  for(unsigned int i=0; i<terms[0]; i++){
    scalar_coeff.push_back(coeff[i]);
    scalar_exp.push_back(exp[i]);
    scalar_coeff_guess.push_back(coeff_guess[i]);
    scalar_exp_guess.push_back(exp_guess[i]);
  }
  for(unsigned int i=terms[0]; i<terms[0]+terms[1]; i++){
    vec_coeff.push_back(coeff[i]);
    vec_exp.push_back(exp[i]);
    vec_coeff_guess.push_back(coeff_guess[i]);
    vec_exp_guess.push_back(exp_guess[i]);
  }

  pair<double, double> mstarns=get_mass_eff_scalar_density( degeneracy,  nucleon_mass, scalar_coeff,  scalar_exp,  saturation_density, false);
  double new_binding_energy=nucleon_mass-eps_over_n(degeneracy, mstarns.first, saturation_density, scalar_coeff,  scalar_exp, vec_coeff, vec_exp, mstarns.second);
  double binding_energy_minimum =eps_over_n_dn(degeneracy, mstarns.first, saturation_density, scalar_coeff,  scalar_exp, vec_coeff, vec_exp, mstarns.second);
  double incompress=incsolv( degeneracy,  nucleon_mass,saturation_density, scalar_coeff,  scalar_exp, vec_coeff, vec_exp);
  pair<double, double> CP_qgp, CP_lg;
  if(abs(binding_energy-new_binding_energy)/binding_energy>criterium){
    cout<<"not validated"<<endl;
    validated=false;
  }
  try{
    CP_lg= get_crit( nucleon_mass, saturation_density,  binding_energy,  degeneracy,  scalar_exp,  vec_exp, false, 19,0.44*saturation_density,scalar_coeff,vec_coeff);
  }catch(int exception){
    cout<<"killed in first CP"<<endl;
    return false;
  }
  double CP_lg_P=T_press_solv(degeneracy, nucleon_mass,CP_lg.second,  scalar_coeff, scalar_exp, vec_coeff, vec_exp, CP_lg.first)/pow(197.3,3);

  if(abs(CP_lg.first-critical_temperature_lg)/critical_temperature_lg>criterium||saturation_density*abs(CP_lg.second/saturation_density- critical_density_lg/saturation_density)/critical_density_lg>criterium){
    cout<<"not validated"<<endl;
    validated=false;
  }
  try{
    CP_qgp= get_crit( nucleon_mass, saturation_density,  binding_energy,  degeneracy,  scalar_exp,  vec_exp, false, critical_temperature_qgp*1.01,critical_density_qgp*1.01,scalar_coeff,vec_coeff);
  }catch(int exception){
    cout<<"killed in second CP"<<endl;
    return false;
  }
  double CP_qgp_P=T_press_solv(degeneracy, nucleon_mass,CP_qgp.second,  scalar_coeff, scalar_exp, vec_coeff, vec_exp, CP_qgp.first)/pow(197.3,3);
  if(abs(CP_qgp.first-critical_temperature_qgp)/critical_temperature_qgp>criterium||saturation_density*abs(CP_qgp.second/saturation_density- critical_density_qgp/saturation_density)/critical_density_qgp>criterium){
    cout<<"not validated"<<endl;
    validated=false;
  }

  FILE * pFile;
  if(validated){
    pFile = fopen ((filename+".csv").c_str(),"a");
  }else{
    pFile = fopen ((filename+"NV.csv").c_str(),"a");
  }

  fprintf(pFile,"%d, ",timestamp);
  for(unsigned int i=0; i<scalar_coeff_guess.size();i++){
    fprintf(pFile, "%.16e, ",scalar_coeff_guess[i]*(pow(saturation_density,scalar_exp_guess[i]-1)));
  }
  for(unsigned int i=0; i<scalar_exp_guess.size();i++){
    fprintf(pFile, "%.16f, ",scalar_exp_guess[i]);
  }
  for(unsigned int i=0; i<vec_coeff_guess.size();i++){
    fprintf(pFile, "%.16e, ",vec_coeff_guess[i]*(pow(saturation_density,vec_exp_guess[i]-1)));
  }
  for(unsigned int i=0; i<vec_exp_guess.size();i++){
    fprintf(pFile, "%.16f, ",vec_exp_guess[i]);
  }
  for(unsigned int i=0; i<scalar_coeff.size();i++){
    fprintf(pFile, "%.16e, ",scalar_coeff[i]*(pow(saturation_density,scalar_exp[i]-1)));
  }
  for(unsigned int i=0; i<scalar_exp.size();i++){
    fprintf(pFile, "%.16f, ",scalar_exp[i]);
  }
  for(unsigned int i=0; i<vec_coeff.size();i++){
    fprintf(pFile, "%.16e, ",vec_coeff[i]*(pow(saturation_density,vec_exp[i]-1)));
  }
  for(unsigned int i=0; i<vec_exp.size();i++){
    fprintf(pFile, "%.16f, ",vec_exp[i]);
  }
  fprintf(pFile, "%.5e, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.16f,%.16f, %.16f ; \n ", error, mstarns.first, new_binding_energy, binding_energy_minimum, incompress, CP_lg.first, CP_lg.second/saturation_density, CP_lg_P, CP_qgp.first, CP_qgp.second/saturation_density, CP_qgp_P,spinodial_l_density/saturation_density, spinodial_r_density/saturation_density);

  fclose(pFile);

  return validated;
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

//generate a set of initial guesses from a CSV file with latin hypercube points
vector<array<double,8> >  latin_gen(uniform_real_distribution<double> sgn_sampler [12], string latinfile, void * parameters, double boundaries_model  [8][2], double num_test, int sgns[4]){
  random_device rd;
  ranlux48_base re(rd());
  double saturation_density=((struct interaction_params_2crit *) parameters)->saturation_density;
  std::ifstream file(latinfile+".csv");
  CSVRow row;
  vector<array<double,8> >  samples;
  unsigned int count=0;
  while(count<num_test&&file >> row)
  {
    count++;

    array<double,8>sample={{stod(row[0].data()),stod(row[1].data()),stod(row[2].data()),stod(row[3].data()),stod(row[4].data()),stod(row[5].data()),stod(row[6].data()),stod(row[7].data())}};

    for(int i=0;i<4;i++){
      sample[i]=sample[i]*(boundaries_model[i][1]-boundaries_model[i][0])+boundaries_model[i][0];
    }
    double buffer;
    if (sample[3]<sample[2]){
      buffer=sample[2];
      sample[2]=sample[3];
      sample[3]=buffer;
    }
    for(int i=4;i<8;i++){
      if(!sgns||sgns[0]==0){
        sample[i]=sgn(sgn_sampler[i+4](re))*(sample[i]*(boundaries_model[i][1]-boundaries_model[i][0])+boundaries_model[i][0])/(pow(saturation_density,sample[i-4]-1));
      }else{
        sample[i]=sgn(sgns[i-4])*(sample[i]*(boundaries_model[i][1]-boundaries_model[i][0])+boundaries_model[i][0])/(pow(saturation_density,sample[i-4]-1));
      }
    }
    samples.push_back(sample);
  }
  return samples;
}

//read in one of the generated latin hypercube points and test it
//the returner conatains the GSL error, the GSL coefficients and exponents as
//well as a copy of the initial guesses
tuple<tuple<double, vector<double>, vector<double>,bool>, vector<double>, vector<double>> latin_test_8D(vector<array<double,8> > samples,void * parameters, unsigned int test_id){
  tuple<tuple<double, vector<double>,vector<double>,bool>, vector<double>, vector<double>> returner={get_interaction_2crit(parameters, false,{samples[test_id][0], samples[test_id][1],samples[test_id][2],samples[test_id][3]}, {samples[test_id][4], samples[test_id][5],samples[test_id][6],samples[test_id][7]}),{samples[test_id][0], samples[test_id][1],samples[test_id][2],samples[test_id][3]},{samples[test_id][4], samples[test_id][5],samples[test_id][6],samples[test_id][7]}};
  return returner;
}

//generate random points and test them
//the returner conatains the GSL error, the GSL coefficients and exponents as
//well as a copy of the initial guesses
tuple<tuple<double, vector<double>, vector<double>,bool>, vector<double>, vector<double>> random_test_8D(uniform_real_distribution<double> sampler [12], void * parameters, int sgns [4]){
  random_device rd;
  ranlux48_base re(rd());
  double saturation_density=((struct interaction_params_2crit *) parameters)->saturation_density;
  double exp1_guess=sampler[0](re);
  double exp2_guess=sampler[1](re);
  double exp3_guess=sampler[2](re);
  if(exp2_guess>exp3_guess){
    double buffer=exp3_guess;
    exp3_guess=exp2_guess;
    exp2_guess=buffer;
  }
  double exp4_guess=sampler[3](re);
  double coeff1_guess, coeff2_guess,coeff3_guess,coeff4_guess;

  if(!sgns||sgns[0]==0){
    coeff1_guess=sgn(sampler[8](re))*sampler[4](re)/(pow(saturation_density,exp1_guess-1));
    coeff2_guess=sgn(sampler[9](re))*sampler[5](re)/(pow(saturation_density,exp2_guess-1));
    coeff3_guess=sgn(sampler[10](re))*sampler[6](re)/(pow(saturation_density,exp3_guess-1));
    coeff4_guess=sgn(sampler[11](re))*sampler[7](re)/(pow(saturation_density,exp4_guess-1));
  }else{
    coeff1_guess=sgn(sgns[0])*sampler[4](re)/(pow(saturation_density,exp1_guess-1));
    coeff2_guess=sgn(sgns[1])*sampler[5](re)/(pow(saturation_density,exp2_guess-1));
    coeff3_guess=sgn(sgns[2])*sampler[6](re)/(pow(saturation_density,exp3_guess-1));
    coeff4_guess=sgn(sgns[3])*sampler[7](re)/(pow(saturation_density,exp4_guess-1));
  }

  tuple<tuple<double, vector<double>,vector<double>,bool>, vector<double>, vector<double>> returner={get_interaction_2crit(parameters, false,{exp1_guess, exp2_guess,exp3_guess,exp4_guess}, {coeff1_guess, coeff2_guess,coeff3_guess, coeff4_guess}),{exp1_guess,exp2_guess,exp3_guess,exp4_guess},{coeff1_guess,coeff2_guess,coeff3_guess,coeff4_guess}};

  return returner;
}



//take given input values and try to find the best parameters
void interaction_8D(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy ,double critical_temperature_lg, double critical_density_lg,double critical_temperature_qgp, double critical_density_qgp, double spinodial_l_density, double spinodial_r_density, double boundaries_model  [8][2], unsigned int terms [2],string filename, bool print, int num_sol, unsigned int num_test, double acceptance, bool del, double criterium , int sgns [4], string latinfile){

  gsl_set_error_handler(&gsl_handler);
  double eps=numeric_limits<double>::epsilon();
  int timestamp=time(0);
  if(print||del){
    FILE * pFile;
    pFile = fopen ((filename+".csv").c_str(),"w");
    fprintf(pFile,"timestamp, ");
    for(unsigned int i=0; i<terms[0];i++){
      fprintf(pFile, "scalar_coeff_guess[%d], ",i);
    }
    for(unsigned int i=0; i<terms[0];i++){
      fprintf(pFile, "scalar_exp_guess[%d], ",i);
    }
    for(unsigned int i=0; i<terms[1];i++){
      fprintf(pFile, "vector_coeff_guess[%d], ", i);
    }
    for(unsigned int i=0; i<terms[1];i++){
      fprintf(pFile, "vector_exp_guess[%d], ",i);
    }
    for(unsigned int i=0; i<terms[0];i++){
      fprintf(pFile, "scalar_coeff[%d], ",i);
    }
    for(unsigned int i=0; i<terms[0];i++){
      fprintf(pFile, "scalar_exp[%d], ",i);
    }
    for(unsigned int i=0; i<terms[1];i++){
      fprintf(pFile, "vector_coeff[%d], ", i);
    }
    for(unsigned int i=0; i<terms[1];i++){
      fprintf(pFile, "vector_exp[%d], ",i);
    }

    fprintf(pFile,"GSL_error, mstar, binding_energy, binding_energy_minimum, incompress, CP_LG_T, CP_LG_n, CP_LG_P, CP_QGP_T, CP_QGP_n, CP_QGP_P, SPIN_L_n, SPIN_R_n; \n");
    fclose(pFile);
  }
  if(num_sol!=0&&!getenv("OMP_CANCELLATION")){
    cout<<"Please do OMP_CANCELLATION=true export OMP_CANCELLATION"<<endl;
    return;
  }

  struct interaction_params_2crit p={nucleon_mass, degeneracy, saturation_density,binding_energy, critical_temperature_lg, critical_density_lg,critical_temperature_qgp, critical_density_qgp,spinodial_l_density,spinodial_r_density,  terms };

  unsigned int count_sol=0;
  unsigned int count_tests=0;
  vector<vector<double>> solution_storage;
  int num_threads;
  if(omp_get_max_threads()<20){
    num_threads=ceil(omp_get_max_threads()-3);
  }else{
     num_threads=omp_get_max_threads();
  }
  //generate the boundaries for random sampling
  uniform_real_distribution<double>  sampler [12];
  for(int i=0; i<4;i++){
    sampler[i]=uniform_real_distribution<double>(boundaries_model[i][0],boundaries_model[i][1]);
  }
  for(int i=4; i<8;i++){
    sampler[i]=uniform_real_distribution<double>(abs(boundaries_model[i][0]),abs(boundaries_model[i][1]));
  }
  for(int i=8; i<12;i++){
    sampler[i]=uniform_real_distribution<double>(-1,1);
  }
  vector<array<double,8> > samples;
  if(!latinfile.empty()){
    samples=latin_gen(sampler, latinfile, &p, boundaries_model, num_test,sgns);
  }
  #pragma omp parallel num_threads(num_threads) shared(p, solution_storage, num_sol, count_sol, count_tests,samples)
  {
    #pragma omp for
    for(int i=0; i<num_threads; i+=1 ){
      for(int k=0; k<=num_test/num_threads; k+=1){
        count_tests=count_tests+1;
        #pragma omp cancellation point for
        try{
          tuple<tuple<double,vector<double>,vector<double>,bool>, vector<double>, vector<double>>  test;
          if(latinfile.empty()){
            test=random_test_8D(sampler,&p, sgns);
          }else{
            test=latin_test_8D(samples ,&p, count_tests);

          }
          tuple<double,vector<double>,vector<double>,bool> result=get<0>(test);
          cout<<"calculated result of "<<count_tests<<" it is "<<get<0>(result)<<endl;
          //check if we have already found a similiar solution
          if(get<0>(result)<acceptance&&(count_sol<num_sol||num_sol<=0)){
            double exp1_cut=round_to_n_digits((get<1>(result)[0]),2);
            double exp2_cut=round_to_n_digits((get<1>(result))[1],3);
            double exp3_cut=round_to_n_digits((get<1>(result))[2],4);
            double exp4_cut=round_to_n_digits((get<1>(result))[3],5);
            double coeff1_cut=round_to_n_digits((get<2>(result))[0],2);
            double coeff2_cut=round_to_n_digits((get<2>(result))[1],3);
            double coeff3_cut=round_to_n_digits((get<2>(result))[2],4);
            double coeff4_cut=round_to_n_digits((get<2>(result))[3],5);

            bool found=false;
            # pragma omp critical
            {
              for(unsigned int j=0; j<solution_storage.size(); j++){
                if(abs(solution_storage[j][0]-exp1_cut)<eps &&abs(solution_storage[j][1]-exp2_cut)<eps&&abs(solution_storage[j][2]-exp3_cut)<eps&&abs(solution_storage[j][3]-exp4_cut)<eps&&abs(solution_storage[j][4]-coeff1_cut)<eps &&abs(solution_storage[j][5]-coeff2_cut)<eps&&abs(solution_storage[j][6]-coeff3_cut)<eps&&abs(solution_storage[j][7]-coeff4_cut)<eps){
                  found=true;
                  break;
                }else{

                }
              }
            }
            if(!found||num_sol<=0){
              bool validated=validate_interaction_2crit(get<0>(result), get<1>(test),get<2>(test), get<1>(result), get<2>(result), filename, timestamp, &p, criterium);
              if(!validated){
                continue;
              }

              solution_storage.push_back({exp1_cut, exp2_cut,exp3_cut, exp4_cut,coeff1_cut, coeff2_cut,coeff3_cut, coeff4_cut});
              count_sol+=1;

              if((count_sol>=num_sol&&num_sol>0)||count_tests>num_test){
                #pragma omp cancel for
              }

            }
          }else{
            continue;
          }
        }catch(int exception){
          continue;
        }
      }
    }
  }

  return;
}

  //find solutions in a grid of input values
void interaction_2crit_grid(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy, double spinodial_l_density, double spinodial_r_density, double boundaries_model  [8][2], double boundaries_crit [4][3], unsigned int terms [2],string filename, bool print, unsigned int num_sol,unsigned int num_test, string latinfile){
  int count=0;
  assert(terms[0]+terms[1]==4);
  string basepath=filesystem::current_path().string() +"/"+filename;
  filesystem::path bpath=basepath;
  filesystem::create_directory(bpath);
  for(double crit_lg_T=boundaries_crit[0][0]; crit_lg_T<=boundaries_crit[0][1];crit_lg_T+=boundaries_crit[0][2] ){
    for(double crit_lg_n=boundaries_crit[1][0]; crit_lg_n<=boundaries_crit[1][1];crit_lg_n+=boundaries_crit[1][2] ){
      for(double crit_qgp_T=boundaries_crit[2][0]; crit_qgp_T<=boundaries_crit[2][1];crit_qgp_T+=boundaries_crit[2][2] ){
        for(double crit_qgp_n=boundaries_crit[3][0]; crit_qgp_n<=boundaries_crit[3][1];crit_qgp_n+=boundaries_crit[3][2] ){
          string dirpath=filesystem::current_path().string() +"/"+filename+"/"+to_string(count);
          filesystem::path path=dirpath;
          filesystem::create_directory(path);
          interaction_8D_2step(nucleon_mass, binding_energy, saturation_density, degeneracy ,crit_lg_T,crit_lg_n,crit_qgp_T,crit_qgp_n, spinodial_l_density, spinodial_r_density,  boundaries_model,terms ,dirpath+"/solutions"+to_string(count), true, num_sol, num_test, 1e-9,latinfile);
          count++;
        }
      }
    }
  }
  return;
}

//use two additonal step to first find approximate parameters and refine them afterwards
void interaction_8D_2step(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy ,double critical_temperature_lg, double critical_density_lg,double critical_temperature_qgp, double critical_density_qgp, double spinodial_l_density, double spinodial_r_density, double boundaries_model  [8][2], unsigned int terms [2],string filename, bool print, int num_sol, unsigned int num_test, double acceptance,string latinfile ){
  bool del=true;
  interaction_8D(nucleon_mass, binding_energy, saturation_density, degeneracy, critical_temperature_lg, critical_density_lg, critical_temperature_qgp, critical_density_qgp, spinodial_l_density, spinodial_r_density, boundaries_model, terms, filename, print, max(5*num_sol,5), num_test, 1, true, 0.05, nullptr, latinfile);
  vector<array<array<double,2>,8> > init=extract_approx(filename, {0.05,0.1,3,4,0.1,0.5,30,1e5});
  vector<array<int,4>>signs=extract_signs(filename);
  int sgns[4];
  double  model [8][2];
  if(init.size()==0){
    cout<<"NO SOLUTION FOUND"<<endl;
    return;
  }
  for(int i=0; i<init.size(); i++){
    del=(i==0);
    memcpy(model, init[i].data(), sizeof(init[i]));
    memcpy(sgns, signs[i].data(), sizeof(signs[i]));
    interaction_8D(nucleon_mass, binding_energy, saturation_density, degeneracy, critical_temperature_lg, critical_density_lg, critical_temperature_qgp,  critical_density_qgp,spinodial_l_density, spinodial_r_density, model, terms, filename+"step2", false, 2, num_test/init.size(), 1e-5, del,0.005,sgns, latinfile);
  }
  if(init.size()==0){
    cout<<"NO SOLUTION FOUND"<<endl;
    return;
  }
  init=extract_approx(filename);
  signs=extract_signs(filename);
  for(int i=0; i<init.size(); i++){
    memcpy(model, init[i].data(), sizeof(init[i]));
    memcpy(sgns, signs[i].data(), sizeof(signs[i]));
    del=(i==0);
    interaction_8D(nucleon_mass, binding_energy, saturation_density, degeneracy, critical_temperature_lg, critical_density_lg, critical_temperature_qgp, critical_density_qgp, spinodial_l_density, spinodial_r_density, model, terms, filename+"step3", false, 1, num_test/init.size(), 1e-9, del,0.0001,sgns, latinfile);
  }
}

//extract approximations from acsv file
vector<array<array<double,2>,8> > extract_approx(string filename, vector<double> range){
  std::ifstream file(filename+".csv");
  CSVRow row;
  vector<array<array<double,2>,8> >  boundaries;
  file >> row;
  while(file >> row)
  {

    double coeff_1=abs(stod(row[9].data()) );
    double coeff_2=abs(stod(row[10].data()));
    double coeff_3=abs(stod(row[11].data()));
    double coeff_4=abs(stod(row[12].data()));
    array<array<double,2>,8> boundaries_model={{{stod(row[13].data())-range[0],stod(row[13].data())+range[0]},{stod(row[14].data())-range[1],stod(row[14].data())+range[1]},{stod(row[15].data())-range[2],stod(row[15].data())+range[2]},{stod(row[16].data())-range[3],stod(row[16].data())+range[3]},{coeff_1/(1+range[4]),coeff_1*(1+range[4])},{coeff_2/(1+range[5]),coeff_2*(1+range[5])},{coeff_3/(1+range[6]),coeff_3*(1+range[6])},{coeff_4/(1+range[7]),coeff_4*(1+range[7])}}};
    boundaries.push_back(boundaries_model);
  }
  return boundaries;
}

//extract signs from the solutions
vector<array<int,4>> extract_signs(string filename){
  std::ifstream file(filename+".csv");
  CSVRow row;
  vector<array<int,4> >  signs;
  file >> row;
  while(file >> row)
  {

    int sgn_1=sgn(stod(row[9].data()) );
    int sgn_2=sgn(stod(row[10].data()));
    int sgn_3=sgn(stod(row[11].data()));
    int sgn_4=sgn(stod(row[12].data()));

    array<int,4> sgn_model={{sgn_1,sgn_2,sgn_3,sgn_4}};
    signs.push_back(sgn_model);
  }
  return signs;
}
