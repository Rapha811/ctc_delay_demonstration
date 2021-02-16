#include "DelayComputation.h"

using namespace std;
using namespace ibex;
using namespace tubex;

bool build_cn_and_contract_delay(IntervalVector& delay_in_out, IntervalVector& tdoa, const Tube& y_, const Tube& e_){

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    cout << "in: " << delay_in_out << endl;

    tdoa[0] = delay_in_out[1] - delay_in_out[0];
    tdoa[1] = delay_in_out[2] - delay_in_out[0];
    tdoa[2] = delay_in_out[2] - delay_in_out[1];

    // Defining custom-built contractors:

    // Sum of signals (y: reception, xi: independant signals, ai: attenuation)
    Function f_y ("a1", "a2", "a3", "x1", "x2", "x3", "y", "a1*x1+a2*x2+a3*x3-y");
    CtcFunction ctc_y(f_y, Interval(0.0)); // f=0

    Function f_minus("a1", "a2", "a1-a2");
    CtcFunction ctc_greater(f_minus, Interval::POS_REALS);

    // Defining domains

    // Copy of tubes, that may be contracted
    Tube y(y_);
    Tube e(e_);
    TubeVector v_r(y_.tdomain(), 3); // 3 received signals, independant
    v_r.sample(y_); // same sampling as y for tubes computations

    // Delays

    // [ Assumption on the delay ]
    tdoa[0] &= Interval::POS_REALS;
    tdoa[1] &= Interval::POS_REALS;

    // Attenuation
    IntervalVector attenuation (3, Interval(0,1));

    // Solver

    ContractorNetwork cn;
    cn.set_fixedpoint_ratio(0.0);

    // Constraint: decomposition of the received signal
    cn.add(ctc_y, {attenuation[0], attenuation[1], attenuation[2], v_r[0],v_r[1],v_r[2],y});

    // Constraint: attenuation of first signal must be the smallest, i.e. a0 must be larger than a1 and a2
    cn.add(ctc_greater, {attenuation[0], attenuation[1]});
    cn.add(ctc_greater, {attenuation[0], attenuation[2]});

    cn.add(ctc_greater, {delay_in_out[1], delay_in_out[0]});
    cn.add(ctc_greater, {delay_in_out[2], delay_in_out[0]});

    // WLOG: one signal MUST arrive before the other or they arrive at the same time. we just have to remember this later during the localization.
    cn.add(ctc_greater, {delay_in_out[2], delay_in_out[1]});
    cn.add(ctc_greater, {attenuation[1], attenuation[2]});

    cn.contract(false);

    // The CtcDelay is not added to the CN to improve computation time (the CN does not know WHEN to call CtcDelay and calls it too often)
    CtcDelay ctc_delay;
    bool continue_contraction = true;
    double old_volume_delay = delay_in_out.volume();
    double old_volume_attenuation = attenuation.volume();
    double old_volume_vr = v_r.volume();

    int i = 0;
    while(continue_contraction){

        continue_contraction = false;
        old_volume_delay = delay_in_out.volume();
        old_volume_attenuation = attenuation.volume();
        old_volume_vr = v_r.volume();

        // Constraints: delays between the emitted signal and the three signal components
        ctc_delay.contract(delay_in_out[0], e, v_r[0]);
        ctc_delay.contract(delay_in_out[1], e, v_r[1]);
        ctc_delay.contract(delay_in_out[2], e, v_r[2]);

        // If CtcDelay does not contract (or results in the empty set), we do not need to call the CN
        if(!delay_in_out.is_empty() && !attenuation.is_empty() && !v_r.is_empty() &&
                (old_volume_delay > delay_in_out.volume() || old_volume_attenuation > attenuation.volume() || old_volume_vr > v_r.volume())){

            // If there is a contraction, we trigger all contractors in the CN and contract
            cn.trigger_all_contractors();
            cn.contract(false);

            // We continue contraction until no more contraction is achieved
            if(!delay_in_out.is_empty() && !attenuation.is_empty() && !v_r.is_empty()){
                continue_contraction = true;
            }
        }

        i++;
    }

    cout << i << " iterations for contraction!" << endl;

    tdoa[0] = delay_in_out[1] - delay_in_out[0];
    tdoa[1] = delay_in_out[2] - delay_in_out[0];
    tdoa[2] = delay_in_out[2] - delay_in_out[1];

    // Constraints: delays between three respective receptions
    ctc_delay.contract(tdoa[0], v_r[0], v_r[1]);
    ctc_delay.contract(tdoa[1], v_r[0], v_r[2]);
    ctc_delay.contract(tdoa[2], v_r[1], v_r[2]);

    cout << "contracted: " << endl;
    cout << "delay: " << delay_in_out << endl;
    cout << "tdoa: " << tdoa << endl;
    cout << "attenuation: " << attenuation << endl;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    cout << endl;

    return cn.emptiness() || tdoa.is_empty();
}


void compute_delay_subpaving(IntervalVector& delay_in_out, const Tube& y_, const Tube& e_, float precision, vector<IntervalVector>& subpaving_delay,
                             vector<IntervalVector>& subpaving_tdoa, IntervalVector& tdoa_out){


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    IntervalVector hull_delay(3, Interval::EMPTY_SET);
    IntervalVector hull_tdoa(3, Interval::EMPTY_SET);

    stack<IntervalVector> s;
    s.push(delay_in_out);

    subpaving_delay.clear();
    subpaving_tdoa.clear();

    while (!s.empty()) {

        IntervalVector box=s.top();
        s.pop();

        IntervalVector tdoa(3);
        bool emptiness = build_cn_and_contract_delay(box, tdoa, y_, e_);

        if(emptiness){
            continue;
        }
        if (box.max_diam()<precision) {

            hull_delay |= box;
            hull_tdoa |= tdoa;

            subpaving_delay.push_back(box);
            subpaving_tdoa.push_back(tdoa);

        } else{
            pair<IntervalVector,IntervalVector> p=box.bisect(box.extr_diam_index(false));
            s.push(p.first);
            s.push(p.second);
        }
    }

    delay_in_out = hull_delay;
    tdoa_out = hull_tdoa;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time for delay computation = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

}
