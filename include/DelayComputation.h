#ifndef DELAYCOMPUTATION_H
#define DELAYCOMPUTATION_H

#include <codac.h>
#include <codac-rob.h>
#include "ibex.h"
#include <chrono>

void compute_delay_subpaving(ibex::IntervalVector& delay_in_out, const codac::Tube& y_, const codac::Tube& e_, float precision, std::vector<ibex::IntervalVector>& subpaving_delay,
                             std::vector<ibex::IntervalVector>& subpaving_tdoa, ibex::IntervalVector& tdoa_out, const codac::TrajectoryVector& v_r_traj);

#endif // DELAYCOMPUTATION_H
