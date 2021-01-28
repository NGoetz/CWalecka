#include <stdlib.h>
#include <utility>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf.h>
using namespace std;
int mass_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_mass_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s);
pair<double, double> get_mass_eff_scalar_density(double degeneracy, double nucleon_mass,double a,double number_density,double mass_scalar, bool print);
int T_mass_eff_mu_eff_scalar_density_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_T_mass_eff_mu_eff_scalar_density (size_t iter, gsl_multiroot_fsolver * s);
tuple<double, double,double,bool> get_T_mass_eff_mu_eff_scalar_density(double a,double b,double nucleon_mass,double number_density, double temperature, double degeneracy, double mass_scalar, bool print);
int coeff_root(const gsl_vector * x, void *params, gsl_vector * f);
int print_state_coeff (size_t iter, gsl_multiroot_fsolver * s);
pair<double, double> get_coeff(double nucleon_mass, double saturation_density, double binding_energy, double degeneracy, double mass_scalar, double mass_vector, bool print);
int crit_root(const gsl_vector * x, void *params, gsl_vector * f);
pair<double, double> get_crit(double nucleon_mass,double saturation_density, double binding_energy, double mass_scalar, double mass_vector, double degeneracy, bool print);
double energy_pp_minus_mass_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b);
double energy_pp_minus_mass_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b);
double press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b);
double incsolv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b);
double T_press_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double b, double temperature);
double  T_press_dmu_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double b, double temperature, double delta=1.);
double  T_press_dT_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a, double b, double temperature, double delta=1e-2);
double  T_press_dn_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double b, double temperature, double delta=1e-2);
double  T_press_dn2_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar, double a,double b, double temperature, double delta=100);
double T_eps_solv(double degeneracy, double nucleon_mass,double number_density, double mass_vector, double mass_scalar,double a, double b, double temperature);
