#include "DelayComputation.h"

using namespace std;
using namespace ibex;
using namespace tubex;

bool build_cn_and_contract_delay(IntervalVector& delay_in_out, IntervalVector& tdoa, const Tube& y_, const Tube& e_, const tubex::TrajectoryVector& v_r_traj){

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

    // Constraints: delays between emission and reception
//    cn.add(ctc::delay, {delay_in_out[0], e, v_r[0]});
//    cn.add(ctc::delay, {delay_in_out[1], e, v_r[1]});
//    cn.add(ctc::delay, {delay_in_out[2], e, v_r[2]});

    // Constraint: decomposition of the received signal
    cn.add(ctc_y, {attenuation[0], attenuation[1], attenuation[2], v_r[0],v_r[1],v_r[2],y});

    // Constraint: attenuation of first signal must be the smallest, i.e. a0 must be larger than a1 and a2
    cn.add(ctc_greater, {attenuation[0], attenuation[1]});
    cn.add(ctc_greater, {attenuation[0], attenuation[2]});

    cn.add(ctc_greater, {delay_in_out[1], delay_in_out[0]});
    cn.add(ctc_greater, {delay_in_out[2], delay_in_out[0]});

    // WLOG: one signal MUST arrive before the other or they arrive at the same time. we just have to remember this later during the localization.
    cn.add(ctc_greater, {delay_in_out[1], delay_in_out[2]});

    cn.contract(false);

    CtcDelay ctc_delay;
    bool continue_contraction = true;
    double old_volume_delay = delay_in_out.volume();
    double old_volume_attenuation = attenuation.volume();
    double old_volume_vr = v_r.volume();

    cn.contract(false);

    int i = 0;
    while(continue_contraction){

        continue_contraction = false;
        old_volume_delay = delay_in_out.volume();
        old_volume_attenuation = attenuation.volume();
        old_volume_vr = v_r.volume();

        ctc_delay.contract(delay_in_out[0], e, v_r[0]);
        ctc_delay.contract(delay_in_out[1], e, v_r[1]);
        ctc_delay.contract(delay_in_out[2], e, v_r[2]);

        VIBesFigTube fig_signals_after_delay(std::to_string(i*2));
        fig_signals_after_delay.add_trajectories(&v_r_traj, "x", "black");
        fig_signals_after_delay.add_tube(&v_r[0], "vr0", "[#ff000064]");
        fig_signals_after_delay.add_tube(&v_r[1], "vr1", "[#00800064]");
        fig_signals_after_delay.add_tube(&v_r[2], "vr2", "[#0000ff64]");
        fig_signals_after_delay.show();
        fig_signals_after_delay.axis_limits(25.2, 61.2,-0.555562881137756, 1.244933115499084);
        fig_signals_after_delay.save_image("", "png", "/home/raphael/Arbeit/Automatica_Time_Delay/latex/figures/contraction_figures_temp");

        if(!delay_in_out.is_empty() && !attenuation.is_empty() && !v_r.is_empty() &&
                (old_volume_delay > delay_in_out.volume() || old_volume_attenuation > attenuation.volume() || old_volume_vr > v_r.volume())){

            cn.trigger_all_contractors();
            cn.contract(false);

            if(!delay_in_out.is_empty() && !attenuation.is_empty() && !v_r.is_empty()){
                continue_contraction = true;
            }
        }

        VIBesFigTube fig_signals_after_cn(std::to_string(i*2+1));
        fig_signals_after_cn.add_trajectories(&v_r_traj, "x", "black");
        fig_signals_after_cn.add_tube(&v_r[0], "vr0", "[#ff000064]");
        fig_signals_after_cn.add_tube(&v_r[1], "vr1", "[#00800064]");
        fig_signals_after_cn.add_tube(&v_r[2], "vr2", "[#0000ff64]");
        fig_signals_after_cn.show();
        fig_signals_after_cn.axis_limits(25.2, 61.2,-0.555562881137756, 1.244933115499084);
        fig_signals_after_cn.save_image("", "png", "/home/raphael/Arbeit/Automatica_Time_Delay/latex/figures/contraction_figures_temp");

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
                             vector<IntervalVector>& subpaving_tdoa, IntervalVector& tdoa_out,const tubex::TrajectoryVector& v_r_traj){


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
        bool emptiness = build_cn_and_contract_delay(box, tdoa, y_, e_, v_r_traj);

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
