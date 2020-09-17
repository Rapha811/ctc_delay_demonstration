#include "TdoaLocalizationPaving.h"

using namespace ibex;
using namespace tubex;
using namespace std;


bool build_cn_and_contract(const IntervalVector& box_a, // emission position
                           IntervalVector& box_b, // emission position
                           const Interval& sea_bounds, // sea levels (surface+seabed)
                           const IntervalVector& tdoa_in,
                           const Interval speed)
{

    // Domains

    IntervalVector tdoa (tdoa_in);
    IntervalVector box_a_ (box_a);

    Interval sea_bounds_ub = sea_bounds.ub();
    Interval sea_bounds_lb = sea_bounds.lb();


    // Defining custom-built contractors:

    Variable a(2);
    Variable b(2);
    Variable surface, surface1, surface2;
    Variable tdoa_var;


    Function direct_distance(a, b, sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) ));
    Function reflected_distance(a, b, surface, sqrt(sqr(a[0]-b[0]) + sqr(a[1]-(-b[1]+2*surface))));

    Function relative_distance1(a, b, surface, tdoa_var, reflected_distance(a,b,surface) - direct_distance(a,b) - tdoa_var*speed);
    Function relative_distance2(a, b, surface1, surface2, tdoa_var, reflected_distance(a,b,surface1) - reflected_distance(a,b,surface2) - tdoa_var*speed);

    Variable tdoa_var_array(3);
    Function f1_1(a, b, tdoa_var_array, relative_distance1(a,b,sea_bounds_ub,tdoa_var_array[0]));
    Function f1_2(a, b, tdoa_var_array, relative_distance1(a,b,sea_bounds_ub,tdoa_var_array[1]));
    CtcFunction ctc1_1(f1_1);
    CtcFunction ctc1_2(f1_2);
    CtcUnion ctc1 (ctc1_1,ctc1_2);

    Function f2_1(a, b, tdoa_var_array, relative_distance1(a,b,sea_bounds_lb,tdoa_var_array[0]));
    Function f2_2(a, b, tdoa_var_array, relative_distance1(a,b,sea_bounds_lb,tdoa_var_array[1]));
    CtcFunction ctc2_1(f2_1);
    CtcFunction ctc2_2(f2_2);
    CtcUnion ctc2 (ctc2_1,ctc2_2);

    Function f3_1(a, b, tdoa_var_array, relative_distance2(a,b,sea_bounds_lb,sea_bounds_ub,tdoa_var_array[2]));
    Function f3_2(a, b, tdoa_var_array, relative_distance2(a,b,sea_bounds_ub,sea_bounds_lb,tdoa_var_array[2]));
    CtcFunction ctc3_1(f3_1);
    CtcFunction ctc3_2(f3_2);
    CtcUnion ctc3 (ctc3_1,ctc3_2);


    // Solver

    ContractorNetwork cn;
    cn.set_fixedpoint_ratio(0.1);

    cn.add(ctc1, {box_a_, box_b, tdoa});
    cn.add(ctc2, {box_a_, box_b, tdoa});
    cn.add(ctc3, {box_a_, box_b, tdoa});

    cn.contract(false); // use 'true' for verbose mode

    return cn.emptiness();

}


void TdoaLocalizationPaving::compute(const IntervalVector& box_a, // emission position
                                     const Interval& sea_bounds, // sea levels (surface+seabed)
                                     const Interval speed,
                                     float precision, // SIVIA precision
                                     VIBesFigMap& fig_map,
                                     const IntervalVector &tdoa_in) // for graphics
{

    IntervalVector box_b = box();

    bool emptiness = build_cn_and_contract(box_a, box_b, sea_bounds, tdoa_in, speed);

    // SIVIA

    if(emptiness)
    {
        set_value(SetValue::OUT);
        fig_map.draw_box(box(), "#0059AF[#0059AF55]");
    }

    else if(box().max_diam() < precision)
    {
        set_value(SetValue::MAYBE);
        fig_map.draw_box(box_b, "#FFE207[#FFE20755]");
    }

    else
    {
        bisect();
        ((TdoaLocalizationPaving*)m_first_subpaving)->compute(box_a, sea_bounds, speed, precision, fig_map, tdoa_in);
        ((TdoaLocalizationPaving*)m_second_subpaving)->compute(box_a, sea_bounds, speed, precision, fig_map, tdoa_in);
    }

    // todo: inner test?
}
