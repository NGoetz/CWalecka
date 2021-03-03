#include <stdlib.h>
#include <utility>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf.h>
#include <string>
using namespace std;
int mass_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_mass_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s);
pair<double, double> get_mass_eff_scalar_density(double degeneracy, double nucleon_mass,double a, double exp_a, double number_density,double mass_scalar, bool print);
int T_mass_eff_mu_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_T_mass_eff_mu_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s);
tuple<double, double,double,bool> get_T_mass_eff_mu_eff_scalar_density(double a,double exp_a, double b,double exp_b,double nucleon_mass,double number_density, double temperature, double degeneracy, double mass_scalar, bool print);
int coeff_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_coeff (size_t iter, gsl_multiroot_fsolver * s);
pair<double, double> get_coeff(double nucleon_mass, double saturation_density, double binding_energy, double degeneracy, double mass_scalar, double mass_vector, double exp_a, double exp_b, bool print, bool * success=NULL);
int crit_root(const gsl_vector * x, void *params, gsl_vector * f);
pair<double, double> get_crit(double nucleon_mass,double saturation_density, double binding_energy, double mass_scalar, double mass_vector, double degeneracy, double exp_a, double exp_b, bool print, double crit_T=19.09, double crit_density=540000, double a=0.0, double b=0.0);
int interaction_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_interaction (size_t iter, gsl_multiroot_fsolver * s);
tuple<double, double,double,double, double, bool> get_interaction_4D(void * p, bool print,double int_exp_a=2.05, double int_exp_b=2.05,double int_a=10, double int_b=10);
double energy_pp_minus_mass_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a,double b, double exp_b);
double energy_pp_minus_mass_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b);
double press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b);
double incsolv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b,double delta=1);
double T_press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double exp_a, double b, double exp_b, double temperature);
double  T_press_dmu_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta=1.);
double  T_press_dT_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta=1e-2);
double  T_press_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta=1e3);
double  T_press_dn2_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b,double temperature, double delta=1e3);
pair<double, double> T_press_dn12_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double exp_a, double b, double exp_b, double temperature, double delta=1e3);
double T_eps_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double exp_a, double b, double exp_b, double temperature);
void gsl_handler (const char * reason, const char * file, int line, int gsl_errno);
bool validate_interaction(double error, double a_guess, double b_guess, double exp_a_guess, double exp_b_guess, double a_in, double b_in, double exp_a, double exp_b, string filename, int timestamp, void* params);
void interaction_4D_grid(double nucleon_mass,double critical_temperature, double critical_density, double binding_energy, double saturation_density, double degeneracy, double mass_scalar, double mass_vector, double  boundaries [4][3],string filename, bool print, int num_sol=0);
void interaction_4D_crit_grid(double nucleon_mass, double binding_energy, double saturation_density, double degeneracy, double mass_scalar, double mass_vector, double boundaries_model  [4][3], double boundaries_crit [2][3], string filename, bool print, int num_sol);
