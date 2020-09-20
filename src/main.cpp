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

    // Sampling time of the signal
    double dt = 0.2;

    // Time domain over which the signal is simulated. Outside this domain the signal is assumed to be zero (to allow the CtcDelay to work)
    Interval tdomain(-10,10);

    // Simulated signal
    TFunction signal ("(2/(sqrt(3*1.5)*3.14)^(1/4))*(1-(t^2)/(1.5^2))*exp(-(t^2)/(2*(1.5^2)))");
//    TFunction signal ("cos(t)*sin(t+1)*exp(-(t^2)/t)");

    // Environment
    IntervalVector sea({{-7.,7.},{-14.,-1.}});

    // Robots
    Vector pa{6.,-10.}, pb{-6.,-8.};

    // speed of signal
    double velocity = 0.2;

    // attenuation of signal
    double attenuation_coefficient = 0.1;

    // Setting up the simulation
    SoundSimulation my_simulation (dt, tdomain, sea, pa, pb, signal, velocity, attenuation_coefficient);

    vibes::beginDrawing();
    VIBesFigMap fig_map("Map");
    fig_map.set_properties(100, 10, 1200, 600);
    VIBesFigTube fig_signals("Signals");
    fig_signals.set_properties(100, 700, 1200, 400);

    my_simulation.draw_map(fig_map);
    fig_map.axis_limits(sea[0].lb(),sea[0].ub(),sea[1].lb(),sea[1].ub());
    my_simulation.draw_signals(fig_signals);

//    return EXIT_SUCCESS;

    vector<IntervalVector> subpaving_delay;
    vector<IntervalVector> subpaving_tdoa;
    IntervalVector tdoa(3);
    IntervalVector delay (3, my_simulation.getInit_delay());
//    delay[0] = Interval(118.4130802819196, 121.0410646866565);
//    delay[1] = Interval(123.6690490913932, 126.2970334961301);
//    delay[2] = Interval(121.0410646866564, 123.6690490913933);

    double sivia_precision_delay = 5.0;
    compute_delay_subpaving(delay, my_simulation.getY(), my_simulation.getE(), sivia_precision_delay, subpaving_delay, subpaving_tdoa, tdoa);

    if(subpaving_delay.empty()){
        cout << "no valid delay found!" << endl;
        return EXIT_SUCCESS;
    }

//    tdoa = my_simulation.true_tdoa + IntervalVector(3,Interval(-1,1));
//    subpaving_tdoa.push_back(tdoa);

    cout << "delay: " << delay << endl;
    cout << "tdoa: " << tdoa << endl;

    cout << "subpaving: " << endl;
    for (std::size_t i = 0; i < subpaving_delay.size(); ++i) {
        cout << "delay: " << subpaving_delay[i] << endl;
        cout << "tdoa: " << subpaving_tdoa[i] << endl;
    }

    bool use_tdoa_subpaving_for_paving = false;
    double sivia_precision_localization = 0.2;
    TdoaLocalizationPaving pav(sea); // paver
    pav.compute(pa, sea[1], Interval(velocity), sivia_precision_localization, fig_map, subpaving_tdoa, tdoa, use_tdoa_subpaving_for_paving);

    std::vector<ConnectedSubset> subsets = pav.get_connected_subsets();

    for (size_t i = 0; i < subsets.size(); ++i) {
        double area_covered = 0.0;
        std::vector<IntervalVector> boxes = subsets[i].get_boxes();
        for (size_t j = 0; j < boxes.size(); ++j) {
            area_covered += boxes[i][0].diam()*boxes[i][1].diam();
        }
        cout << "subset " << i << ": " << area_covered << "m²" << endl;

    }

    vibes::endDrawing();
    return EXIT_SUCCESS;
}
