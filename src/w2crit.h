#include <stdlib.h>
#include <utility>
#include <string>
#include <random>
#include <array>

using namespace std;

bool validate_interaction_2crit(double error, vector<double> exp_guess, vector<double> coeff_guess,vector<double> exp, vector<double> coeff, string filename, int timestamp, void* params, double criterium=0.001);
tuple<tuple<double, vector<double>, vector<double>,bool>, vector<double>, vector<double>> random_test_8D(uniform_real_distribution<double> sampler [12], void * parameters,int  sgns [4]);
void interaction_8D(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy ,double critical_temperature_lg, double critical_density_lg,double critical_temperature_qgp, double critical_density_qgp, double spinodial_l_density, double spinodial_r_density, double boundaries_model  [8][2], unsigned int terms [2],string filename, bool print, int num_sol, unsigned int tests, double acceptance=1e-3,bool del=true , double criterium=0.001, int  sgns [4] =nullptr, string latinfile="");
void interaction_2crit_grid(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy,double boundaries_model  [8][2], double boundaries_crit [6][3], unsigned int terms [2],string filename, bool print, unsigned int num_sol,unsigned int num_test, string latinfile="");
void interaction_8D_2step(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy ,double critical_temperature_lg, double critical_density_lg,double critical_temperature_qgp, double critical_density_qgp, double spinodial_l_density, double spinodial_r_density, double boundaries_model  [8][2], unsigned int terms [2],string filename, bool print, int num_sol, unsigned int num_test, double acceptance,string latinfile="" );
pair<vector<array<array<double,2>,8> >,vector<array<int,4>>> extract_approx(string filename, vector<double> range={0.001,0.005});
tuple<tuple<double, vector<double>, vector<double>,bool>, vector<double>, vector<double>> latin_test_8D(vector<array<double,8> > samples,void * parameters, unsigned int test_id);
vector<array<double,8> >  latin_gen(uniform_real_distribution<double> sgn_sampler [12], string latinfile, void * parameters, double boundaries_model  [8][2], double num_test, int sgns[4]);
