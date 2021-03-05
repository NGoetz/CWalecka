#include <vector>
using namespace std;
double fmom(double degeneracy, double number_density);
double inv_p(double degeneracy,double p);
double E_eff(double p, double mass_eff);
double scalar_density(double degeneracy, double mass_eff,double number_density);
double scalar_density_dn(double degeneracy, double mass_eff,double number_density,vector<double> scalar_coeff,vector<double> scalar_exp);
double mass_eff_from_scalar_density(double nucleon_mass,vector<double> scalar_coeff,vector<double> scalar_exp, double scalar_density);
double mass_eff_from_scalar_density_dn(double mass_eff,  double number_density,  vector<double> scalar_coeff,vector<double> scalar_exp, double degeneracy);
double eps_over_n(double degeneracy, double mass_eff,double number_density,vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density);
double eps_over_n_dn(double degeneracy, double mass_eff,double number_density,vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density);
double press(double degeneracy, double mass_eff, double number_density, vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density);
double fermiplus(double temperature, double mass_eff, double p, double mu_eff);
double fermiminus(double temperature, double mass_eff, double p, double mu_eff);
double muB_from_number_density(double degeneracy, double mass_eff, double number_density, vector<double> vec_coeff, vector<double> vec_exp);
double T_mu_eff(double number_density, double muB, vector<double> vec_coeff, vector<double> vec_exp);
double T_pressur_integrand (double p, void * params);
double T_scalar_density_integrand (double p, void * params);
double T_number_density_integrand (double p, void * params);
double T_pressure(double degeneracy, double temperature,double mass_eff,double mu_eff,vector<double> scalar_coeff, vector<double> scalar_exp, vector<double> vec_coeff, vector<double> vec_exp, double scalar_density, double number_density);
double T_scalar_density(double degeneracy,double temperature,double mass_eff,double mu_eff);
double T_number_density(double degeneracy,double temperature,double mass_eff,double mu_eff);


