#include<stdlib.h>
#include<vector>
using namespace std;

struct plot_params
{
    double binding_energy;
    double saturation_density;
    double nucleon_mass;
    double degeneracy;
    vector<double> scalar_exp;
    vector<double> vec_exp;
    vector<double> scalar_coeff;
    vector<double> vec_coeff;
    double temperature;
    
};

void energy_pp_plot(double start, double step, double nmax, void * params, int num);
void energy_pp_dn_plot(double start, double step, double nmax, void * params);
void press_plot(double start, double step, double nmax, void * params);
void scalar_density_plot(double start, double step, double nmax, void * params);
void eff_mass_plot_p(double start, double step, double nmax,void * params, int num );
void T_press_plot(double start, double step, double nmax, void * params, int num);
void T_energy_pp_plot(double start, double step, double nmax, void * params);
void T_mass_eff_plot(double start, double step, double nmax, void * params);
void interactions_plot(double * lowerbound, double * upperbound, double step, void * params);
