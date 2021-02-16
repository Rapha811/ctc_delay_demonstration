#include "TdoaLocalizationPaving.h"

using namespace ibex;
using namespace tubex;
using namespace std;


bool build_cn_and_contract(const IntervalVector& box_a, // emission position
                           IntervalVector& box_b, // emission position
                           const Interval& sea_bounds, // sea levels (surface+seabed)
                           const std::vector<IntervalVector>& tdoa_subpaving_in,
                           const IntervalVector& tdoa_in,
                           const Interval speed,
                           bool use_subpaving)
{

    // Domains

    Interval sea_bounds_ub = sea_bounds.ub();
    Interval sea_bounds_lb = sea_bounds.lb();

    std::vector<IntervalVector> tdoa_vec;
    if(use_subpaving){
        tdoa_vec = tdoa_subpaving_in;
    } else {
        tdoa_vec.push_back(tdoa_in);
    }

    // Defining custom-built contractors:

    Variable a(2);
    Variable b(2);
    Variable surface, surface1, surface2;
    Variable tdoa_var;


    Function direct_distance(a, b, sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) ));
    Function reflected_distance(a, b, surface, sqrt(sqr(a[0]-b[0]) + sqr(a[1]-(-b[1]+2*surface))));

    Function relative_distance1(a, b, surface, tdoa_var, reflected_distance(a,b,surface) - direct_distance(a,b) - tdoa_var*speed);
    Function relative_distance2(a, b, surface1, surface2, tdoa_var, reflected_distance(a,b,surface1) - reflected_distance(a,b,surface2) - tdoa_var*speed);

    IntervalVector box_b_out(2, Interval::EMPTY_SET);

    for (size_t i = 0; i < tdoa_vec.size(); ++i) {

        IntervalVector tdoa_cur (tdoa_vec[i]);

        Function f1_1(b, relative_distance1(box_a,b,sea_bounds_ub,tdoa_cur[0]));
        CtcFunction ctc1_1(f1_1);
        Function f2_1(b, relative_distance1(box_a,b,sea_bounds_lb,tdoa_cur[1]));
        CtcFunction ctc2_2(f2_1);
        Function f3_1(b, relative_distance2(box_a,b,sea_bounds_lb,sea_bounds_ub,tdoa_cur[2]));
        CtcFunction ctc3_1(f3_1);
        CtcCompo ctc_compo1 (ctc1_1, ctc2_2, ctc3_1);

        Function f1_2(b, relative_distance1(box_a,b,sea_bounds_ub,tdoa_cur[1]));
        CtcFunction ctc1_2(f1_2);
        Function f2_2(b, relative_distance1(box_a,b,sea_bounds_lb,tdoa_cur[0]));
        CtcFunction ctc2_1(f2_2);
        Function f3_2(b, relative_distance2(box_a,b,sea_bounds_ub,sea_bounds_lb,tdoa_cur[2]));
        CtcFunction ctc3_2(f3_2);
        CtcCompo ctc_compo2 (ctc1_2, ctc2_1, ctc3_2);

        CtcUnion ctc_union(ctc_compo1, ctc_compo2);
        CtcFixPoint ctc_fix (ctc_union, 0.1);

        IntervalVector box_b_ (box_b);
        ctc_fix.contract(box_b_);

        box_b_out |= box_b_;

    }

    box_b = box_b_out;
    return box_b.is_empty();
}


void TdoaLocalizationPaving::compute(const IntervalVector& box_a, // emission position
                                     const Interval& sea_bounds, // sea levels (surface+seabed)
                                     const Interval speed,
                                     float precision, // SIVIA precision
                                     VIBesFigMap& fig_map, // for graphics
                                     const std::vector<IntervalVector>& tdoa_subpaving_in,
                                     const IntervalVector& tdoa_in,
                                     bool use_subpaving)
{

    IntervalVector box_b = box();

    bool emptiness = build_cn_and_contract(box_a, box_b, sea_bounds, tdoa_subpaving_in, tdoa_in, speed, use_subpaving);

    // SIVIA

    if(emptiness)
    {
        set_value(SetValue::OUT);
        fig_map.draw_box(box(), "#666666ff[#cccccc55]");
    }

    else if(box().max_diam() < precision)
    {
        set_value(SetValue::MAYBE);
        fig_map.draw_box(box(), "#009b00ff[#009b00b4]");
    }

    else
    {
        bisect();
        ((TdoaLocalizationPaving*)m_first_subpaving)->compute(box_a, sea_bounds, speed, precision, fig_map, tdoa_subpaving_in, tdoa_in, use_subpaving);
        ((TdoaLocalizationPaving*)m_second_subpaving)->compute(box_a, sea_bounds, speed, precision, fig_map, tdoa_subpaving_in, tdoa_in, use_subpaving);
    }

    // todo: inner test?
}
