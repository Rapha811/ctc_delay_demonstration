#ifndef DELAYCOMPUTATION_H
#define DELAYCOMPUTATION_H

#include <tubex.h>
#include <tubex-rob.h>
#include "ibex.h"
#include <chrono>

void compute_delay_subpaving(ibex::IntervalVector& delay_in_out, const tubex::Tube& y_, const tubex::Tube& e_, float precision,
                                       std::vector<ibex::IntervalVector>& subpaving_delay, std::vector<ibex::IntervalVector>& subpaving_tdoa, ibex::IntervalVector& tdoa_out);

#endif // DELAYCOMPUTATION_H
