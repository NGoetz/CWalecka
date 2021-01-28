#include <../../src/wplot.h>
#include <../../src/wsolvers.h>
#include <../../src/wfunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include <iostream>
#include <matplot/matplot.h>
#include <vector>

//the convergence factor from fm to MeV
double const conv=197.3;

// convenience for Python-like behaviour
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

// what follows are convenience plotting functions
//energy density at T=0; input in multiples of saturation density
void energy_pp_plot(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=energy_pp_minus_mass_solv(degeneracy, nucleon_mass,t[i]*saturation_density, mass_vector, mass_scalar, res.first, res.second);
    }
    auto ax = matplot::gca();
    matplot::plot(t,s)->line_width(3);
    matplot::hold(matplot::on);
    matplot::title("Energy/nucleon");
    matplot::xlabel("n/n_0");
    matplot::ylabel("epsilon/n -m_N");
    ax->y_axis().label_font_size(12);
    ax->y_axis().label_weight("bold");
    ax->x_axis().label_font_size(12);
    ax->x_axis().label_weight("bold");
   
    matplot::hold(matplot::off);
    matplot::show();

    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
//derivative of energy density at T=0
void energy_pp_dn_plot(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=energy_pp_minus_mass_dn_solv(degeneracy, nucleon_mass,t[i]*saturation_density, mass_vector, mass_scalar, res.first, res.second);
    }

    auto ax = matplot::gca();
    
    matplot::plot(ax,t,s);
    
    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
//pressure at T=0
void press_plot(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=press_solv(degeneracy, nucleon_mass,t[i]*saturation_density, mass_vector, mass_scalar, res.first, res.second);
    }

    auto ax = matplot::gca();
    matplot::ylim(matplot::manual);
    matplot::ylim({*min_element(s.begin(), s.end()),*max_element(s.begin(), s.end())});
    
    matplot::plot(ax,t,s);
    
    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
//scalar density at T=0
void scalar_density_plot(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=get_mass_eff_scalar_density(degeneracy, nucleon_mass,res.first,t[i]*saturation_density, mass_scalar, false).second;
    }

    auto ax = matplot::gca();
    matplot::plot(ax,t,s);
    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
//effective mass at T=0, as function of p instead of n
void eff_mass_plot_p(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=get_mass_eff_scalar_density(degeneracy, nucleon_mass,res.first,t[i]*saturation_density, mass_scalar, false).first/nucleon_mass;
        t[i]=fmom(degeneracy,t[i]*saturation_density);
    }

    auto ax = matplot::gca();
    matplot::plot(ax,t,s);
    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
//pressure at finite T, as function of n at T=0
void T_press_plot(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    double temperature=(( struct plot_params *) params )->temperature;
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=1e10*T_press_solv(degeneracy, nucleon_mass, t[i]*saturation_density, mass_vector, mass_scalar, res.first, res.second, temperature)/pow((conv),3);
    }
       auto ax = matplot::gca();
   
    matplot::ylim(matplot::manual);
    matplot::ylim({*min_element(s.begin(), s.end()),*max_element(s.begin(), s.end())});
    
    matplot::plot(ax,t,s)->line_width(3);
    matplot::hold(matplot::on);
    matplot::title("Pressure");
    matplot::xlabel("n/n_0");
    matplot::ylabel("P");
    ax->y_axis().label_font_size(12);
    ax->y_axis().label_weight("bold");
    ax->x_axis().label_font_size(12);
    ax->x_axis().label_weight("bold");
    matplot::hold(matplot::off);
    matplot::show();
    
    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}

//energy density at finite T, as function of n at T=0
void T_energy_pp_plot(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    double temperature=(( struct plot_params *) params )->temperature;
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=T_eps_solv(degeneracy, nucleon_mass, t[i]*saturation_density, mass_vector, mass_scalar, res.first, res.second,temperature)/(t[i]*saturation_density)-nucleon_mass;
    }
    auto ax = matplot::gca();
   
    matplot::ylim(matplot::manual);
    matplot::ylim({*min_element(s.begin(), s.end()),*max_element(s.begin(), s.end())});
    
    matplot::plot(ax,t,s)->line_width(3);
    matplot::hold(matplot::on);
    matplot::title("Energy density");
    matplot::xlabel("n/n_0");
    matplot::ylabel("epsilon/n-m");
    ax->y_axis().label_font_size(12);
    ax->y_axis().label_weight("bold");
    ax->x_axis().label_font_size(12);
    ax->x_axis().label_weight("bold");
    matplot::hold(matplot::off);
    matplot::show();
    
    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}

void T_mass_eff_plot(double start, double step, double nmax, void * params){
    double binding_energy = ((struct plot_params *) params)->binding_energy;
    double saturation_density=((struct plot_params *) params)->saturation_density;
    double nucleon_mass=((struct plot_params *) params)->nucleon_mass;
    double degeneracy=(( struct plot_params *) params )->degeneracy;
    double mass_scalar=((struct plot_params *) params)->mass_scalar;
    double mass_vector=(( struct plot_params *) params )->mass_vector;
    double temperature=(( struct plot_params *) params )->temperature;
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    pair<double, double> res=get_coeff(nucleon_mass, saturation_density, binding_energy, degeneracy, mass_scalar, mass_vector, false);
    for (int i=0; i<t.size(); i=i+1){
        tuple<double, double,double,bool> res2=get_T_mass_eff_mu_eff_scalar_density(res.first,res.second, nucleon_mass, t[i]*saturation_density,  temperature, degeneracy, mass_scalar, false);
        s[i]=get<0>(res2);
        t[i] = fmom(degeneracy,t[i]*saturation_density);
    }
   
    auto ax = matplot::gca();
    matplot::ylim(matplot::manual);
    matplot::ylim({*min_element(s.begin(), s.end()),*max_element(s.begin(), s.end())});
    
    matplot::plot(ax,t,s);
    
    do 
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}


