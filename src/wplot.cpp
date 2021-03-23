#include <wplot.h>
#include <wsolvers.h>
#include <wfunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include <iostream>
#include <matplot/matplot.h>
#include <vector>
#include <cassert>

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
void energy_pp_plot(double start, double step, double nmax, void * params, int num){
    auto ax = matplot::gca();
    for (int i=0; i<num; i++){
        double binding_energy =(((struct plot_params *) params))[i].binding_energy;
        double saturation_density=(((struct plot_params *) params))[i].saturation_density;
        double nucleon_mass=(((struct plot_params *) params))[i].nucleon_mass;
        double degeneracy=(( (struct plot_params *) params))[i].degeneracy;
        vector<double> scalar_exp=(((struct plot_params *) params))[i].scalar_exp;
        vector<double> vec_exp=(( (struct plot_params *) params)) [i].vec_exp;
        vector<double> scalar_coeff=(((struct plot_params *) params))[i].scalar_coeff;
        vector<double> vec_coeff=(( (struct plot_params *) params)) [i].vec_coeff;
        auto t = arange<double>(start, nmax,step);
        auto s = arange<double>(start, nmax, step);
        for (int i=0; i<t.size(); i=i+1){
            s[i]=energy_pp_minus_mass_solv(degeneracy, nucleon_mass,t[i]*saturation_density, scalar_coeff, scalar_exp, vec_coeff, vec_exp);
        }

        matplot::plot(t,s)->line_width(3);
        matplot::hold(matplot::on);
        matplot::title("Energy/nucleon at T=0");
        matplot::xlabel("n/n_0");
        matplot::ylabel("epsilon/n -m_N (MeV)");
        ax->y_axis().label_font_size(12);
        ax->y_axis().label_weight("bold");
        ax->x_axis().label_font_size(12);
        ax->x_axis().label_weight("bold");


    }
    auto l=matplot::legend(ax, {"V1S1", "V2S0", "V0S2"});
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
    vector<double> scalar_exp=(((struct plot_params *) params))->scalar_exp;
    vector<double> vec_exp=(( (struct plot_params *) params))->vec_exp;
    vector<double> scalar_coeff=(((struct plot_params *) params))->scalar_coeff;
    vector<double> vec_coeff=(( (struct plot_params *) params))->vec_coeff;
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=energy_pp_minus_mass_dn_solv(degeneracy, nucleon_mass,t[i]*saturation_density, scalar_coeff, scalar_exp, vec_coeff, vec_exp);
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
    vector<double> scalar_exp=(((struct plot_params *) params))->scalar_exp;
    vector<double> vec_exp=(( (struct plot_params *) params))->vec_exp;
    vector<double> scalar_coeff=(((struct plot_params *) params))->scalar_coeff;
    vector<double> vec_coeff=(( (struct plot_params *) params))->vec_coeff;
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    auto sdn = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=press_solv(degeneracy, nucleon_mass,t[i]*saturation_density, scalar_coeff, scalar_exp, vec_coeff, vec_exp);
        sdn[i]=press_dn_solv(degeneracy, nucleon_mass,t[i]*saturation_density, scalar_coeff, scalar_exp, vec_coeff, vec_exp)*1e5;
    }
    matplot::hold(matplot::on);
    auto ax = matplot::gca();

    matplot::plot(ax,t,s)->line_width(3);;
    matplot::plot(ax,t,sdn)->line_width(3);;
    auto l=matplot::legend(ax, {"pressure", "derivative"});
    matplot::hold(matplot::off);
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
    vector<double> scalar_exp=(((struct plot_params *) params))->scalar_exp;
    vector<double> vec_exp=(( (struct plot_params *) params))->vec_exp;
    vector<double> scalar_coeff=(((struct plot_params *) params))->scalar_coeff;
    vector<double> vec_coeff=(( (struct plot_params *) params))->vec_coeff;
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=get_mass_eff_scalar_density(degeneracy, nucleon_mass,scalar_coeff , scalar_exp, t[i]*saturation_density, false).second;
    }

    auto ax = matplot::gca();
    matplot::plot(ax,t,s);
    do
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
//effective mass at T=0, as function of p instead of n
void eff_mass_plot_p(double start, double step, double nmax, void * params, int num){
    auto ax = matplot::gca();
    for (int i=0; i<num; i++){
        double binding_energy =(((struct plot_params *) params))[i].binding_energy;
        double saturation_density=(((struct plot_params *) params))[i].saturation_density;
        double nucleon_mass=(((struct plot_params *) params))[i].nucleon_mass;
        double degeneracy=(( (struct plot_params *) params))[i].degeneracy;
        vector<double> scalar_exp=(((struct plot_params *) params))[i].scalar_exp;
        vector<double> vec_exp=(( (struct plot_params *) params))[i].vec_exp;
        vector<double> scalar_coeff=(((struct plot_params *) params))[i].scalar_coeff;
        vector<double> vec_coeff=(( (struct plot_params *) params))[i].vec_coeff;
        auto t = arange<double>(start, nmax,step);
        auto s = arange<double>(start, nmax, step);
        for (int i=0; i<t.size(); i=i+1){
            s[i]=get_mass_eff_scalar_density(degeneracy, nucleon_mass,scalar_coeff,scalar_exp,t[i]*saturation_density, false).first/nucleon_mass;
            t[i]=fmom(degeneracy,t[i]*saturation_density);
        }
        matplot::plot(t,s)->line_width(3);
        matplot::hold(matplot::on);
        matplot::title("Effective mass at T=0");
        matplot::xlabel("p_F in fm-1");
        matplot::ylabel("m_eff/m_N");
        ax->y_axis().label_font_size(12);
        ax->y_axis().label_weight("bold");
        ax->x_axis().label_font_size(12);
        ax->x_axis().label_weight("bold");


    }
    matplot::legend(ax, {"V0S2"});
    matplot::hold(matplot::off);
    matplot::show();

    do
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
//pressure at finite T, as function of n at T=0
void T_press_plot(double start, double step, double nmax, void * params, int num){
    auto ax = matplot::gca();
    for (int i=0; i<num; i++){
        cout<<i<<endl;
        double binding_energy =(((struct plot_params *) params))[i].binding_energy;
        double saturation_density=(((struct plot_params *) params))[i].saturation_density;
        double nucleon_mass=(((struct plot_params *) params))[i].nucleon_mass;
        double degeneracy=(( (struct plot_params *) params))[i].degeneracy;
        vector<double> scalar_exp=(((struct plot_params *) params))[i].scalar_exp;
        vector<double> vec_exp=(( (struct plot_params *) params))[i].vec_exp;
        vector<double> scalar_coeff=(((struct plot_params *) params))[i].scalar_coeff;
        vector<double> vec_coeff=(( (struct plot_params *) params))[i].vec_coeff;
        double temperature=(( (struct plot_params *) params)) [i].temperature;
        assert(temperature>0.2);
        auto t = arange<double>(start, nmax,step);
        auto s = arange<double>(start, nmax, step);
        auto v = arange<double>(start, nmax, step);
        auto w = arange<double>(start, nmax, step);

        for (int i=0; i<t.size(); i=i+1){
            //cout<<t[i]<<endl;
            //cout<<" "<<vec_coeff[0]<<" "<<vec_exp[0]<<endl;
            s[i]=T_press_solv(degeneracy, nucleon_mass, t[i]*saturation_density, scalar_coeff, scalar_exp, vec_coeff, vec_exp, temperature)/pow((conv),3);
           // v[i]=1e6*T_press_dn_solv(degeneracy, nucleon_mass, t[i]*saturation_density,scalar_coeff, scalar_exp, vec_coeff, vec_exp, temperature)/pow((conv),3);
            //w[i]=1e11*T_press_dn2_solv(degeneracy, nucleon_mass, t[i]*saturation_density, scalar_coeff, scalar_exp, vec_coeff, vec_exp, temperature)/pow((conv),3);

        }

       //matplot::ylim(matplot::manual);
        //matplot::ylim({0.4,1});

        matplot::hold(matplot::on);
        matplot::plot(ax,t,s)->line_width(2);
        //matplot::plot(ax,t,v)->line_width(3);
        //matplot::plot(ax,t,w)->line_width(3);
        matplot::title("Critical pressure for all 2 term models");
        auto l=matplot::legend(ax, {"V1S1", "V2S0", "V0S2"});
        l->location(matplot::legend::general_alignment::topleft);
        matplot::xlabel("n/n_0");
        matplot::ylabel("P (MeV fm-3)");
        ax->y_axis().label_font_size(12);
        ax->y_axis().label_weight("bold");
        ax->x_axis().label_font_size(12);
        ax->x_axis().label_weight("bold");
    }
    matplot::hold(matplot::off);
    //matplot::legend(ax, {"Walecka solution", "New solution"});
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
    double temperature=(( struct plot_params *) params )->temperature;
    vector<double> scalar_exp=(((struct plot_params *) params))->scalar_exp;
    vector<double> vec_exp=(( (struct plot_params *) params))->vec_exp;
    vector<double> scalar_coeff=(((struct plot_params *) params))->scalar_coeff;
    vector<double> vec_coeff=(( (struct plot_params *) params))->vec_coeff;
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);
    for (int i=0; i<t.size(); i=i+1){
        s[i]=T_eps_solv(degeneracy, nucleon_mass, t[i]*saturation_density,scalar_coeff, scalar_exp, vec_coeff, vec_exp,temperature)/(t[i]*saturation_density)-nucleon_mass;
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
    double temperature=(( struct plot_params *) params )->temperature;
    vector<double> scalar_exp=(((struct plot_params *) params))->scalar_exp;
    vector<double> vec_exp=(( (struct plot_params *) params))->vec_exp;
    vector<double> scalar_coeff=(((struct plot_params *) params))->scalar_coeff;
    vector<double> vec_coeff=(( (struct plot_params *) params))->vec_coeff;
    auto t = arange<double>(start, nmax,step);
    auto s = arange<double>(start, nmax, step);

    for (int i=0; i<t.size(); i=i+1){
        tuple<double, double,double,bool> res2=get_T_mass_eff_mu_eff_scalar_density(scalar_coeff, scalar_exp, vec_coeff, vec_exp,nucleon_mass, t[i]*saturation_density,  temperature, degeneracy, false);
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



void interactions_plot(double * lowerbound, double * upperbound, double step, void * params){
    gsl_set_error_handler(&gsl_handler);
    double nucleon_mass=((struct interaction_params *) params)->nucleon_mass;
    double degeneracy=((struct interaction_params *) params)->degeneracy;
    double saturation_density=((struct interaction_params *) params)->saturation_density;
    double binding_energy=((struct interaction_params *) params)->binding_energy;
    double critical_temperature=((struct interaction_params *) params)->critical_temperature;
    double critical_density=((struct interaction_params *) params)->critical_density;
    auto x = arange<double>(lowerbound[0], upperbound[0],step);
    auto y = arange<double>(lowerbound[1], upperbound[1], step);
    tuple<double, double,double,double,double, bool>res;
    vector<vector<double>> data = {};
    unsigned int terms [2]={1,1};
    for (int i=0; i<x.size(); i=i+1){
        vector<double> storage;
        for (int j=0; j<y.size(); j=j+1){

            try{
                struct interaction_params p={nucleon_mass, degeneracy, saturation_density,binding_energy, critical_temperature, critical_density, terms };
                res= get_interaction_4D(&p, false,{x[i],y[j]},{1e-5,1e-5});
                double x0=get<2>(res);//scalar_coeff
                double x3=get<1>(res);//vec_exp
                double x1=get<3>(res);//vec_coeff
                double x2=get<0>(res);//scalar_exp

                pair<double, double> output=get_mass_eff_scalar_density(degeneracy, nucleon_mass,{x0}, {x2}, saturation_density, false);
                double y0 =binding_energy-nucleon_mass+eps_over_n(degeneracy, output.first, saturation_density, {x0}, {x2}, {x1}, {x3}, output.second);
                double y1 =eps_over_n_dn(degeneracy, output.first, saturation_density,{x0}, {x2}, {x1}, {x3}, output.second);
                double y2 =T_press_dn_solv(degeneracy, nucleon_mass, critical_density,  {x0}, {x2}, {x1}, {x3}, critical_temperature);
                double y3 =T_press_dn2_solv(degeneracy, nucleon_mass, critical_density,  {x0}, {x2}, {x1}, {x3}, critical_temperature);
                storage.push_back(log(abs(y0)+abs(y1)+abs(y2)+abs(y3)));
            }catch(int exception)
        {
            storage.push_back(1e6);
        }


        }
        data.push_back(storage);

    }
    auto ax = matplot::gca();
    vector<string>xs,ys;
    transform(x.begin(),x.end(),back_inserter(xs),[](double const& val){return to_string(val);});
    transform(y.begin(),y.end(),back_inserter(ys),[](double const& val){return to_string(val);});
    ax->x_axis().ticklabels(xs);
    ax->y_axis().ticklabels(ys);
    matplot::heatmap(data);
    matplot::colorbar()
        .limits({-11, 6})
        .tick_values({-14,-11,-8,-5, -3, 0, 3,6})
        .ticklabels({"0","1e-11", "1e-8", "1e-5", "1e-3", "1","1e3","1e6"});
    matplot::show();

    do
    {
        cout << '\n' << "Press a key to continue...";
    } while (cin.get() != '\n');
}
