#include <tubex.h>
#include <tubex-rob.h>
#include "ibex.h"
#include <chrono>

#include "SoundSimulation.h"
#include "DelayComputation.h"
#include "TdoaLocalizationPaving.h"

using namespace std;
using namespace tubex;
using namespace ibex;

int main()
{

    double dt = 0.2;
    Interval tdomain(-40.,40.);
    TFunction signal ("(2/(sqrt(3*1.5)*3.14)^(1/4))*(1-(t^2)/(1.5^2))*exp(-(t^2)/(2*(1.5^2)))");

    // Environment
    
    IntervalVector sea({{-7.,7.},{-7.,-1.}});

    // Robots

    Vector pa{5.,-2.}, pb{-5.,-5.};

    double velocity = 0.2;
    double attenuation_coefficient = 0.1;
    bool only_tdoa = false;

    // Setting up the simulation

    SoundSimulation my_simulation (dt, tdomain, sea, pa, pb, signal, velocity, attenuation_coefficient, only_tdoa);

    vibes::beginDrawing();
    VIBesFigMap fig_map("Map");
    fig_map.set_properties(100, 10, 1200, 600);
    fig_map.axis_limits(sea[0].lb(),sea[0].ub(),sea[1].lb(),sea[1].ub());
    VIBesFigTube fig_signals("Signals");
    fig_signals.set_properties(100, 700, 1200, 400);

    my_simulation.draw_map(fig_map);
    my_simulation.draw_signals(fig_signals);

    Interval init_delay = (Interval(my_simulation.true_delay[0]) | my_simulation.true_delay[1] | my_simulation.true_delay[2]) + Interval(-10,10);
    IntervalVector delay (3, init_delay);

    vector<IntervalVector> subpaving_delay;
    vector<IntervalVector> subpaving_tdoa;
    IntervalVector tdoa(3);

    double sivia_precision = 5.0;
    compute_delay_subpaving(delay, my_simulation.y, my_simulation.e, sivia_precision, subpaving_delay, subpaving_tdoa, tdoa);

//    tdoa = my_simulation.true_tdoa + IntervalVector(3,Interval(-1,1));

    cout << "delay: " << delay << endl;
    cout << "tdoa: " << tdoa << endl;

    cout << "subpaving: " << endl;
    for (std::size_t i = 0; i < subpaving_delay.size(); ++i) {
        cout << "delay: " << subpaving_delay[i] << endl;
        cout << "tdoa: " << subpaving_tdoa[i] << endl;
    }

    TdoaLocalizationPaving pav(sea); // paver
    pav.compute(pa, sea[1], Interval(velocity), 0.2, fig_map, tdoa);

    vibes::endDrawing();
    return EXIT_SUCCESS;
}
