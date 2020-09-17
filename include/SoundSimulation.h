#ifndef SOUNDSIMULATION_H
#define SOUNDSIMULATION_H

#include <tubex.h>
#include <tubex-rob.h>
#include "ibex.h"

class SoundSimulation
{
public:

    ibex::Vector pa = ibex::Vector(2);
    ibex::Vector pb = ibex::Vector(2);

    ibex::IntervalVector sea;

    ibex::Vector pa_surface = ibex::Vector(2);
    ibex::Vector pa_seabed = ibex::Vector(2);

    ibex::Vector true_distances = ibex::Vector(3);
    ibex::Vector true_delay = ibex::Vector(3);
    ibex::Vector true_tdoa = ibex::Vector(3);
    ibex::Vector true_attenuation = ibex::Vector(3);

    tubex::Tube y = tubex::Tube(ibex::Interval(0,1));
    tubex::Tube e = tubex::Tube(ibex::Interval(0,1));

    tubex::Trajectory e_;
    tubex::TrajectoryVector v_r = tubex::TrajectoryVector(3);

    SoundSimulation(double dt, ibex::Interval tdomain, ibex::IntervalVector sea_, ibex::Vector pa_, ibex::Vector pb_,
                    tubex::TFunction signal, double velocity, double attenuation_coefficient, bool only_tdoa);

    void draw_map(tubex::VIBesFigMap& fig_map);

    void draw_signals(tubex::VIBesFigTube& fig_signals);

protected:
    SoundSimulation();

};

#endif // SOUNDSIMULATION_H
