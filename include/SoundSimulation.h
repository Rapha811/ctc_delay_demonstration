#ifndef SOUNDSIMULATION_H
#define SOUNDSIMULATION_H

#include <tubex.h>
#include <tubex-rob.h>
#include "ibex.h"

class SoundSimulation
{
public:

    SoundSimulation(double dt, ibex::Interval tdomain, ibex::IntervalVector sea_, ibex::Vector pa_, ibex::Vector pb_,
                    tubex::TFunction signal, double velocity, double attenuation_coefficient);

    void draw_map(tubex::VIBesFigMap& fig_map);

    void draw_signals(tubex::VIBesFigTube& fig_signals);

    tubex::Tube get_reception_tube_y() const;

    tubex::Tube get_emission_tube_e() const;

    ibex::Interval get_init_delay() const;

protected:
    SoundSimulation();

private:

    // positions of the sound source and sound receiver
    ibex::Vector position_receiver = ibex::Vector(2);
    ibex::Vector position_emitter = ibex::Vector(2);

    // bounds of the sea
    ibex::IntervalVector sea;

    // reflected positions of the sound receiver on the sea surface and seabed
    ibex::Vector position_receiver_surface = ibex::Vector(2);
    ibex::Vector position_receiver_seabed = ibex::Vector(2);

    // the true simulation values are saved
    ibex::Vector true_distances = ibex::Vector(3);
    ibex::Vector true_delay = ibex::Vector(3);
    ibex::Vector true_tdoa = ibex::Vector(3);
    ibex::Vector true_attenuation = ibex::Vector(3);

    // an initial range for the delay must be specified for CtcDelay to work (otherwise the behaviour would be undefined since the signals are not defined infinitely)
    // we select init_delay such that it enclosed the true delay +/- some uncertainty which will be contracted. this is just a very rough assumption.
    ibex::Interval init_delay;

    // reception tube of the three superimposed signals
    tubex::Tube y = tubex::Tube(ibex::Interval(0,1));

    // tube of the emitted signals which we assume to know
    tubex::Tube e = tubex::Tube(ibex::Interval(0,1));

    // emitted signal without uncertainty
    tubex::Trajectory e_;

    // three components of the received signal without uncertainty
    tubex::TrajectoryVector v_r = tubex::TrajectoryVector(3);

};

#endif // SOUNDSIMULATION_H
