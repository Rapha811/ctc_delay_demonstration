#include "SoundSimulation.h"

using namespace ibex;
using namespace tubex;
using namespace std;

// Computes the boxed refraction point between signal emission/reception (for multipaths)
IntervalVector refraction_point(const IntervalVector& pa, // emission position (box)
                                const IntervalVector& pb, // receiver position (box)
                                double y_surface) // coordinate of the surface
{
    // Computing the reflection of the receiver position on the surface
    IntervalVector reflected_pb{pb[0],-pb[1]+2.*y_surface};

    // Computing refraction point on the surface
    Interval r = (y_surface-reflected_pb[1])/(pa[1]-reflected_pb[1]);
    IntervalVector pi{r*(pa[0]-reflected_pb[0])+reflected_pb[0],
                y_surface};
    return pi;
}

// Computes the actual refracted ray: displays the result and returns the propagated distance
void draw_ray_truth(VIBesFigMap& fig_map, // graphical view
                 const Vector& pa, // emission position (box)
                 const Vector& pb, // receiver position (box)
                 double y_surface) // coordinate of the surface
{
    // Refraction point
    Vector pi = refraction_point(pa, pb, y_surface).mid();

    fig_map.draw_circle(pi[0], pi[1], 0.05, "black[black]");
    fig_map.draw_edge(Edge(pa,pi));

    fig_map.draw_circle(pb[0], pb[1], 0.05, "black[black]");
    fig_map.draw_edge(Edge(pb,pi));
}

// Computes the actual refracted ray: displays the result and returns the propagated distance
double compute_true_distances(const Vector& pa, // emission position (box)
                 const Vector& pb, // receiver position (box)
                 double y_surface) // coordinate of the surface
{
    // Refraction point
    Vector pi = refraction_point(pa, pb, y_surface).mid();

    return Edge(pa,pi).length().mid() + Edge(pb,pi).length().mid();
}

tubex::Tube add_signals_and_build_reception_tube(const TrajectoryVector &v_r, const Interval &common_t, double dt){
    vector<Interval> v_tdomains;
    vector<Interval> v_codomains;
    double lb, ub = common_t.lb();

    do
    {
        lb = ub; // we guarantee all slices are adjacent
        ub = min(lb + dt, common_t.ub()); // the tdomain of the last slice may be smaller

        Interval tdomain (lb,ub);
        v_tdomains.push_back(tdomain);
        Interval codomain(0.0);
        for (int i = 0; i < v_r.size(); ++i) {
            if(tdomain.is_subset(v_r[i].tdomain())){
                codomain += v_r[i](tdomain);
            }
            else if(tdomain.intersects(v_r[i].tdomain())){
                Interval intersecting_tdomain = tdomain & v_r[i].tdomain();
                codomain += Interval(0.0) | v_r[i](intersecting_tdomain);
            }
        }
        v_codomains.push_back(codomain);
    } while(ub < common_t.ub());

    return Tube(v_tdomains, v_codomains);
}

tubex::Tube SoundSimulation::get_reception_tube_y() const
{
    return y;
}

tubex::Tube SoundSimulation::get_emission_tube_e() const
{
    return e;
}

ibex::Interval SoundSimulation::get_init_delay() const
{
    return init_delay;
}

SoundSimulation::SoundSimulation(double dt, ibex::Interval tdomain, ibex::IntervalVector sea_, ibex::Vector pa_, ibex::Vector pb_,
                                 tubex::TFunction signal, double velocity, double attenuation_coefficient)
{

    position_receiver = pa_;
    position_emitter = pb_;
    sea = sea_;

    // reflected positions of the receiver
    position_receiver_surface = {position_receiver[0],-position_receiver[1]+2.*sea[1].ub()};
    position_receiver_seabed = {position_receiver[0],-position_receiver[1]+2.*sea[1].lb()};

    // Computing truth
    true_distances[0] = Edge(position_receiver,position_emitter).length().mid(); // direct path
    true_distances[1] = compute_true_distances(position_receiver, position_emitter, sea[1].ub()); // surface
    true_distances[2] = compute_true_distances(position_emitter, position_receiver, sea[1].lb()); // seabed

    // Reference signal ([e]mission)
    e_ = Trajectory(tdomain, signal, dt);

    // Delayed [r]eceived signals
    v_r = TrajectoryVector(3, e_);

    // Memo: shift_tdomain: \f$[t_0,t_f]:=[t_0+a,t_f+a]\f$
    // Memo: delay constraint \f$\mathbf{x}(t)=\mathbf{y}(t+\tau)\f$

    true_delay[0] = true_distances[0]/velocity;
    true_delay[1] = true_distances[1]/velocity;
    true_delay[2] = true_distances[2]/velocity;

    // an initial range for the delay must be specified for CtcDelay to work (otherwise the behaviour would be undefined since the signals are not defined infinitely)
    // we select init_delay such that it enclosed the true delay +/- some uncertainty which will be contracted. this is just a very rough assumption.
    init_delay = (Interval(true_delay[0]) | true_delay[1] | true_delay[2]) + Interval(-10,10);

    cout << "init delay: " << init_delay << endl;

    // Computing shifted signals
    v_r[0].shift_tdomain(true_delay[0]);
    v_r[1].shift_tdomain(true_delay[1]);
    v_r[2].shift_tdomain(true_delay[2]);

    true_tdoa[0] = true_delay[1] - true_delay[0];
    true_tdoa[1] = true_delay[2] - true_delay[0];
    true_tdoa[2] = true_delay[2] - true_delay[1];

    cout << "true delay: " << true_delay << endl;
    cout << "true tdoa: " << true_tdoa << endl;

    true_attenuation[0] = 1/pow(true_distances[0],attenuation_coefficient);
    true_attenuation[1] = 1/pow(true_distances[1],attenuation_coefficient);
    true_attenuation[2] = 1/pow(true_distances[2],attenuation_coefficient);
    v_r[0] *= true_attenuation[0];
    v_r[1] *= true_attenuation[1];
    v_r[2] *= true_attenuation[2];

    cout << "true attenuation: " << true_attenuation << endl;


    // generation of the tubes enclosing the emitted and received signal
    // the signals are padded with zeros inside the computed time domains to allow the CtcDelay to work
    Interval y_tdomain = v_r[0].tdomain() | v_r[1].tdomain() | v_r[2].tdomain();

    Interval y_tdomain_padded = y_tdomain | (tdomain+init_delay);
    y = add_signals_and_build_reception_tube(v_r, y_tdomain_padded, dt);

    Interval e_tdomain_padded = tdomain | (y_tdomain-init_delay);
    e = add_signals_and_build_reception_tube(TrajectoryVector(1, e_), e_tdomain_padded, dt);

    // merge slices in the beginning and end that are very much similar
    double distance_threshold = 1e-3;
    e.merge_similar_slices(distance_threshold);
    y.merge_similar_slices(distance_threshold);

    cout << "v_r0: " << v_r[0] << endl;
    cout << "v_r1: " << v_r[1] << endl;
    cout << "v_r2: " << v_r[2] << endl;
    cout << "e: " << e << endl;
    cout << "y: " << y << endl;

}

void draw_hyperbola(tubex::VIBesFigMap& fig_map, const IntervalVector &sea, double dt, double diff, const Vector &b, const Vector &c, string color){

    double xb = b[0];
    double yb = b[1];

    double xc = c[0];
    double yc = c[1];

    for (double x = sea[0].lb(); x < sea[0].ub(); x+=dt) {

        // solved using MATLAB
        double y1 = (diff*sqrt((xb*xc*-2.0-yb*yc*2.0-diff*diff+xb*xb+xc*xc+yb*yb+yc*yc)*(x*xb*-4.0-x*xc*4.0+xb*xc*2.0-yb*yc*2.0-diff*diff+(x*x)*4.0+xb*xb+xc*xc+yb*yb+yc*yc))+(diff*diff)*yb+(diff*diff)*yc-(xb*xb)*yb+(xb*xb)*yc+(xc*xc)*yb-(xc*xc)*yc+yb*(yc*yc)+(yb*yb)*yc-yb*yb*yb-yc*yc*yc+x*xb*yb*2.0-x*xb*yc*2.0-x*xc*yb*2.0+x*xc*yc*2.0)/(yb*yc*4.0+(diff*diff)*2.0-(yb*yb)*2.0-(yc*yc)*2.0);
        double y2 = (-diff*sqrt((xb*xc*-2.0-yb*yc*2.0-diff*diff+xb*xb+xc*xc+yb*yb+yc*yc)*(x*xb*-4.0-x*xc*4.0+xb*xc*2.0-yb*yc*2.0-diff*diff+(x*x)*4.0+xb*xb+xc*xc+yb*yb+yc*yc))+(diff*diff)*yb+(diff*diff)*yc-(xb*xb)*yb+(xb*xb)*yc+(xc*xc)*yb-(xc*xc)*yc+yb*(yc*yc)+(yb*yb)*yc-yb*yb*yb-yc*yc*yc+x*xb*yb*2.0-x*xb*yc*2.0-x*xc*yb*2.0+x*xc*yc*2.0)/(yb*yc*4.0+(diff*diff)*2.0-(yb*yb)*2.0-(yc*yc)*2.0);

        if(sea[1].contains(y1)) fig_map.draw_point(Point(x,y1), color);
        if(sea[1].contains(y2)) fig_map.draw_point(Point(x,y2), color);
    }
}

void SoundSimulation::draw_map(tubex::VIBesFigMap& fig_map){

    fig_map.draw_box(sea, "#CEEEFF[#CEEEFF]");
    fig_map.draw_vehicle({position_receiver[0],position_receiver[1],M_PI}, 0.8);
    fig_map.draw_vehicle({position_emitter[0],position_emitter[1],0.}, 0.8);

    fig_map.draw_edge(Edge(position_receiver,position_emitter));
    draw_ray_truth(fig_map, position_receiver, position_emitter, sea[1].ub()); // surface
    draw_ray_truth(fig_map, position_emitter, position_receiver, sea[1].lb()); // seabed

    /// can be used to draw three circles corresponding to the distance between emitter and (reflected) receiver positions
//    fig_map.draw_circle(pa[0],pa[1],true_distances[0],"red");
//    fig_map.draw_circle(pa_surface[0],pa_surface[1],true_distances[1],"red");
//    fig_map.draw_circle(pa_seabed[0],pa_seabed[1],true_distances[2],"red");

    /// can be used to draw the hyperbolas corresponding to the pseudo ranges
    /// in practice, we do not know whether true_distances[1] and true_distances[2] correspond to reflections on the surface OR seabed

    // correct correspondences between distances and surface/seabead
    draw_hyperbola(fig_map, sea, 0.01, true_distances[1]-true_distances[0], position_receiver_surface, position_receiver, "blue");
    draw_hyperbola(fig_map, sea, 0.01, true_distances[2]-true_distances[0], position_receiver_seabed, position_receiver, "blue");
    // the third hyperbola is redundant
//    draw_hyperbola(fig_map, sea, 0.01, true_distances[2]-true_distances[1], position_receiver_seabed, position_receiver_surface, "blue");

    // incorrect correspondences between distances and surface/seabead
    draw_hyperbola(fig_map, sea, 0.01, true_distances[2]-true_distances[0], position_receiver_surface, position_receiver, "orange");
    draw_hyperbola(fig_map, sea, 0.01, true_distances[1]-true_distances[0], position_receiver_seabed, position_receiver, "orange");
    // the third hyperbola is redundant
//    draw_hyperbola(fig_map, sea, 0.01, true_distances[1]-true_distances[2], position_receiver_seabed, position_receiver_surface, "orange");

    fig_map.show();

}

void SoundSimulation::draw_signals(tubex::VIBesFigTube& fig_signals){
    fig_signals.add_trajectory(&e_, "e", "#AF0800");
    fig_signals.add_trajectories(&v_r, "e", "blue");
    fig_signals.add_tube(&y, "y", "[magenta]");
    fig_signals.add_tube(&e, "e_", "[orange]");
    Interval x_range(e.first_slice()->tdomain().ub()-5.0,y.last_slice()->tdomain().lb()+5.0);
    Interval y_range = y.codomain() | e.codomain();
    fig_signals.show();
    fig_signals.axis_limits(x_range.lb(),x_range.ub(),y_range.lb(),y_range.ub());
}


SoundSimulation::SoundSimulation(){}
