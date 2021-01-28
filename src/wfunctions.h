//a and b refer to the normalized coupling constants square of the walecka model
double fmom(double degeneracy, double number_density);
double inv_p(double degeneracy,double p);
double E_eff(double p, double mass_eff);
double scalar_density(double degeneracy, double mass_eff,double number_density);
double scalar_density_dn(double degeneracy, double mass_eff,double number_density,double a, double mass_scalar);
double mass_eff_from_scalar_density(double nucleon_mass,double a,double scalar_density,double mass_scalar);
double mass_eff_from_scalar_density_dn(double mass_eff,  double number_density,  double a, double mass_scalar, double degeneracy);
double eps_over_n(double degeneracy, double mass_eff,double number_density,double a,double b,double scalar_density, double mass_vector,double mass_scalar);
double eps_over_n_dn(double degeneracy, double mass_eff,double number_density,double a,double b,double scalar_density, double mass_vector,double mass_scalar);
double eps_over_n_d2n(double degeneracy, double mass_eff,double number_density,double a,double b,double scalar_density, double mass_vector,double mass_scalar);
double incompress(double degeneracy, double mass_eff,double number_density,double a,double b,double scalar_density, double mass_vector,double mass_scalar);
double press(double degeneracy, double mass_eff, double number_density, double a, double b, double scalar_density, double mass_vector, double mass_scalar);
double fermiplus(double temperature, double mass_eff, double p, double mu_eff);
double fermiminus(double temperature, double mass_eff, double p, double mu_eff);
double muB_from_number_density(double degeneracy, double mass_eff, double number_density,double mass_vector,double b);
double T_mu_eff(double number_density, double muB, double b,double mass_vector);
double T_pressur_integrand (double p, void * params);
double T_scalar_density_integrand (double p, void * params);
double T_number_density_integrand (double p, void * params);
double T_pressure(double degeneracy, double temperature,double mass_eff,double mu_eff,double a, double b,double scalar_density, double number_density,double mass_vector, double mass_scalar);
double T_scalar_density(double degeneracy,double temperature,double mass_eff,double mu_eff);
double T_number_density(double degeneracy,double temperature,double mass_eff,double mu_eff);


