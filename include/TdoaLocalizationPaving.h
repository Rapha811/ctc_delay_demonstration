#ifndef TDOALOCALIZATIONPAVING_H
#define TDOALOCALIZATIONPAVING_H

#include "tubex.h"
#include "ibex.h"


class TdoaLocalizationPaving : public tubex::Paving
{
public:
    TdoaLocalizationPaving(const IntervalVector& init_box) : tubex::Paving(init_box)
    {

    }

    void compute(const ibex::IntervalVector& box_a, const ibex::Interval& sea_bounds, const ibex::Interval speed, float precision, tubex::VIBesFigMap& fig_map,
                 const std::vector<ibex::IntervalVector> &tdoa_subpaving_in, const ibex::IntervalVector &tdoa_in, bool use_subpaving);
};

#endif // TDOALOCALIZATIONPAVING_H
