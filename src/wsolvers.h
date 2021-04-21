#include <stdlib.h>
#include <utility>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf.h>
#include <string>
#include <random>
#include <array>

using namespace std;
struct mass_eff_scalar_density_params
{
    vector<double> scalar_coeff;
    vector<double> scalar_exp;
    double number_density;
    double nucleon_mass;
    double degeneracy;
    double temperature;
    vector<double> vec_coeff;
    vector<double> vec_exp;
};
//struct for the parameters for root solving for the coupling constants
struct coeff_params
{
    double saturation_densitiy;
    double nucleon_mass;
    double binding_energy;
    double degeneracy;
    vector<double> scalar_exp;
    vector<double> vec_exp;
    unsigned int* terms;
};
//struct for the parameters for finding the critical point
struct crit_params
{
    vector<double> scalar_coeff;
    vector<double> scalar_exp;
    vector<double> vec_coeff;
    vector<double> vec_exp;
    double nucleon_mass;
    double degeneracy;
};
//struct for the parameters for finding the interaction terms for 1 CP
struct interaction_params
{
    double nucleon_mass;
    double degeneracy;
    double saturation_density;
    double binding_energy;
    double critical_temperature;
    double critical_density;
    unsigned int* terms;
};
//struct for the parameters for finding the interaction terms for 2 CPs
struct interaction_params_2crit
{
    double nucleon_mass;
    double degeneracy;
    double saturation_density;
    double binding_energy;
    double critical_temperature_lg;
    double critical_density_lg;
    double critical_temperature_qgp;
    double critical_density_qgp;
    double spinodial_l_density;
    double spinodial_r_density;
    unsigned int* terms;
};
int mass_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_mass_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s);
pair<double, double> get_mass_eff_scalar_density(double degeneracy, double nucleon_mass,vector<double> scalar_coeff ,vector<double>  scalar_exp,double number_density,bool print);
int T_mass_eff_mu_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_T_mass_eff_mu_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s);
tuple<double, double,double,bool> get_T_mass_eff_mu_eff_scalar_density(vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double nucleon_mass,double number_density, double temperature, double degeneracy, bool print);
int coeff_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_coeff (size_t iter, gsl_multiroot_fsolver * s);
pair<double, double> get_coeff(double nucleon_mass, double saturation_density, double binding_energy, double degeneracy, vector<double> scalar_exp, vector<double> vec_exp, bool print, unsigned int terms[2], bool * success=NULL);
int crit_root(const gsl_vector * x, void *params, gsl_vector * f);
pair<double, double> get_crit(double nucleon_mass,double saturation_density, double binding_energy, double degeneracy, vector<double> scalar_exp, vector<double> vec_exp, bool print, double crit_T=19.09, double crit_density=540000, vector<double> scalar_coeff={0.0}, vector<double> vec_coeff={0.0});
int interaction_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_interaction (size_t iter, gsl_multiroot_fsolver * s);
tuple<double, double,double,double, double, bool> get_interaction_4D(void * p, bool print,vector<double> init_exp={2.05,2.05},vector<double> init_coeff={10*4*M_PI/pow(550,2.0),10*4*M_PI/pow(783,2.0)});
int interaction_root_2crit(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_interaction_2crit (size_t iter, gsl_multiroot_fsolver * s);
tuple<double, vector<double>, vector<double>, bool> get_interaction_2crit(void * p, bool print,vector<double> init_exp={2.05,2.05,3,3},vector<double> init_coeff={10*4*M_PI/pow(550,2.0),10*4*M_PI/pow(783,2.0),0,0});
double energy_pp_minus_mass_solv(double degeneracy, double nucleon_mass,double number_density, vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp);
double energy_pp_minus_mass_dn_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp);
double press_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp);
double press_dn_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp);
double incsolv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp,double delta=1);
double T_press_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature);
double  T_press_dmu_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta=1.);
double  T_press_dT_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta=1e-2);
double  T_press_dn_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta=1e3);
double  T_press_dn2_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp,double temperature, double delta=1e3);
pair<double, double> T_press_dn12_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature, double delta=1e3);
double T_eps_solv(double degeneracy, double nucleon_mass,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double temperature);
void gsl_handler (const char * reason, const char * file, int line, int gsl_errno);
bool validate_interaction(double error, vector<double> coeff_guess, vector<double> exp_guess,vector<double> coeff, vector<double> exp, string filename, int timestamp, void* params);
void interaction_4D_grid(double nucleon_mass,double critical_temperature, double critical_density, double binding_energy, double saturation_density, double degeneracy, double boundaries  [4][3], unsigned int terms [2], string filename, bool print, int num_sol=0);
void interaction_4D_crit_grid(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy, double boundaries_model  [4][3], double boundaries_crit [2][3], unsigned int terms [2], string filename, bool print, int num_sol);
double round_to_n_digits(double x, int n);
