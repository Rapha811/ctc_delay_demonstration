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

SoundSimulation::SoundSimulation(double dt, ibex::Interval tdomain, ibex::IntervalVector sea_, ibex::Vector pa_, ibex::Vector pb_,
                                 tubex::TFunction signal, double velocity, double attenuation_coefficient, bool only_tdoa)
{

    pa = pa_;
    pb = pb_;
    sea = sea_;

    pa_surface = {pa[0],-pa[1]+2.*sea[1].ub()};
    pa_seabed = {pa[0],-pa[1]+2.*sea[1].lb()};

    // Computing truth

    true_distances[0] = Edge(pa,pb).length().mid(); // direct path
    true_distances[1] = compute_true_distances(pa, pb, sea[1].ub()); // surface
    true_distances[2] = compute_true_distances(pb, pa, sea[1].lb()); // seabed

    // Computing shifted signals

    // Reference signal ([e]mission), inspired from Ricker wavelet
    e_ = Trajectory(tdomain, signal, dt);

    // Delayed [r]eceived signals
    v_r = TrajectoryVector(3, e_);

    // Memo: shift_tdomain: \f$[t_0,t_f]:=[t_0+a,t_f+a]\f$
    // Memo: delay constraint \f$\mathbf{x}(t)=\mathbf{y}(t+\tau)\f$

    true_delay[0] = true_distances[0]/velocity;
    true_delay[1] = true_distances[1]/velocity;
    true_delay[2] = true_distances[2]/velocity;

    if(only_tdoa){
        true_delay = true_delay - Vector(3,true_delay[0]);
    }

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

    // Defining all signals on the same t-domain
    Interval common_t = v_r[0].tdomain() & v_r[1].tdomain() & v_r[2].tdomain();
//    Interval common_t (30, 75);
    v_r.truncate_tdomain(common_t);

    // Building reception tube
    y = Tube(v_r[0], dt) + v_r[1] + v_r[2];
    e = Tube(e_,dt);

    // merge slices in the beginning and end that are very much similar
    double distance_threshold = 1e-5;
    e.merge_similar_slices(distance_threshold);
    y.merge_similar_slices(distance_threshold);

}

void SoundSimulation::draw_map(tubex::VIBesFigMap& fig_map){

    fig_map.draw_box(sea, "#CEEEFF[#CEEEFF]");
    fig_map.draw_vehicle({pa[0],pa[1],M_PI}, 0.8);
    fig_map.draw_vehicle({pb[0],pb[1],0.}, 0.8);

    fig_map.draw_edge(Edge(pa,pb));
    draw_ray_truth(fig_map, pa, pb, sea[1].ub()); // surface
    draw_ray_truth(fig_map, pb, pa, sea[1].lb()); // seabed

    fig_map.draw_circle(pa[0],pa[1],true_distances[0],"red");
    fig_map.draw_circle(pa_surface[0],pa_surface[1],true_distances[1],"red");
    fig_map.draw_circle(pa_seabed[0],pa_seabed[1],true_distances[2],"red");

    for (double x = sea[0].lb(); x < sea[0].ub(); x+=0.01) {
        for (double y = sea[1].lb(); y < sea[1].ub(); y+=0.01) {

            double distance = sqrt(pow(x-pa_surface[0],2) + pow(y-pa_surface[1],2)) - sqrt(pow(x-pa[0],2) + pow(y-pa[1],2));
            if(abs(distance - (true_distances[1]-true_distances[0])) < 0.01){
                fig_map.draw_point(Point(x,y), "blue");
            }
            if(abs(distance - (true_distances[2]-true_distances[0])) < 0.01){
                fig_map.draw_point(Point(x,y), "orange");
            }

            double distance1 = sqrt(pow(x-pa_seabed[0],2) + pow(y-pa_seabed[1],2)) - sqrt(pow(x-pa[0],2) + pow(y-pa[1],2));
            if(abs(distance1 - (true_distances[2]-true_distances[0])) < 0.01){
                fig_map.draw_point(Point(x,y), "blue");
            }
            if(abs(distance1 - (true_distances[1]-true_distances[0])) < 0.01){
                fig_map.draw_point(Point(x,y), "orange");
            }

            double distance2 = sqrt(pow(x-pa_surface[0],2) + pow(y-pa_surface[1],2)) - sqrt(pow(x-pa_seabed[0],2) + pow(y-pa_seabed[1],2));
            if(abs(distance2 - (true_distances[1]-true_distances[2])) < 0.01){
                fig_map.draw_point(Point(x,y), "blue");
            }
            if(abs(distance2 - (true_distances[2]-true_distances[1])) < 0.01){
                fig_map.draw_point(Point(x,y), "orange");
            }

        }

    }
    fig_map.show();

}

void SoundSimulation::draw_signals(tubex::VIBesFigTube& fig_signals){
    fig_signals.add_trajectory(&e_, "e", "#AF0800");
    fig_signals.add_trajectories(&v_r, "e", "blue");
    fig_signals.add_tube(&y, "y", "[magenta]");
    fig_signals.add_tube(&e, "e_", "[orange]");
    fig_signals.show();
}


SoundSimulation::SoundSimulation(){}
