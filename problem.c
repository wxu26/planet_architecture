/**
 * Velocity dependent drag force
 *
 * This is a very simple example on how to implement a velocity 
 * dependent drag force. The example uses the IAS15 integrator, which 
 * is ideally suited to handle non-conservative forces.
 * No gravitational forces or collisions are present.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rebound.h"

// global variables: to be set by input later
const double tmax = 1.e5*2.*M_PI;
const double dt_output = tmax/1000.;
const char output_name[] = "sim.out";

const double initial_planet_mass = 3e-6;
const int initial_n_planet = 4;
const double initial_a_ratio = 1.2; // initial semimajor axes ratios: a_{i+1} / a_i
const double initial_a = 1.; // initial semimajor axis of innermost planet

const double Mdot_0_Msunyr = 1.e-8;
const double alpha = 0.01;
const double tau_d_yr = 1.e5;
const double B_kG = 1.;
const double h0 = 0.025;

const double Msun_in_g = 1.9891e33;
const double au_in_cm = 1.49597871e13;
const double G = 1.;

const double M = 1.;

double mass_to_radius(const double m){
    // just scale by earth values for now...
    double r = pow(m/3.e-6, 1./3.)*4.26e-04;
    return r;
}

int collision(struct reb_simulation* const sim, struct reb_collision c){
    printf("collision happened!\n");
    struct reb_particle* particles = sim->particles;
    double m = particles[c.p1].m+particles[c.p2].m;
    double r3 = pow(particles[c.p1].r,3)+pow(particles[c.p2].r,3);
    double px = particles[c.p1].m*particles[c.p1].vx+particles[c.p2].m*particles[c.p2].vx;
    double py = particles[c.p1].m*particles[c.p1].vy+particles[c.p2].m*particles[c.p2].vy;
    double pz = particles[c.p1].m*particles[c.p1].vz+particles[c.p2].m*particles[c.p2].vz;
    double mx = particles[c.p1].m*particles[c.p1].x+particles[c.p2].m*particles[c.p2].x;
    double my = particles[c.p1].m*particles[c.p1].y+particles[c.p2].m*particles[c.p2].y;
    double mz = particles[c.p1].m*particles[c.p1].z+particles[c.p2].m*particles[c.p2].z;
    particles[c.p1].m = m;
    particles[c.p1].r = mass_to_radius(m);
    particles[c.p1].x = mx/m;
    particles[c.p1].y = my/m;
    particles[c.p1].z = mz/m;
    particles[c.p1].vx = px/m;
    particles[c.p1].vy = py/m;
    particles[c.p1].vz = pz/m;
    return 2; // means removing particle 2
}

void migration(struct reb_simulation* const sim){
    const int N = sim->N;
    struct reb_particle* p0 = &sim->particles[0];
    double t = sim->t;
    for (int i=1;i<N;i++){ // only loop through planets
        // orbital elements
        struct reb_particle* p = &sim->particles[i];
        struct reb_orbit o= reb_orbit_from_particle(sim->G, sim->particles[i], sim->particles[0]);
        double a = o.a;
        double e = o.e;
        double I = o.inc;
        double m = p->m;

        // disk properties
        double tau_d = tau_d_yr*2.*M_PI;
        double Mdot_Msunyr = Mdot_0_Msunyr * exp(-t/tau_d);
        double Sigma_gcm2 = 75. * (Mdot_Msunyr/1.e-9) * 1./(alpha/1.e-2) * 1./sqrt(a/0.1);
        double Sigma = Sigma_gcm2 / (Msun_in_g/au_in_cm/au_in_cm);
        double r_c = 0.1 * pow(Mdot_Msunyr/1.e-9, -2./7.) * pow(B_kG,4./7.); // truncation radius
        if (a<r_c) Sigma = 1.e-40;
        double beta = 0.5; // surface density slope, Sigma ~ r^-beta
        double h = h0;

        // migration following Cresswell & Nelson
        double Omega = sqrt(G*M/(a*a*a));
        double t_wave = M/m * M/(Sigma*a*a) * pow(h,4) / Omega;
        double t_e = t_wave/0.78;
        double t_I = t_wave/0.544;
        double t_m = 2*t_wave/(2.7+1.1*beta)/(h*h);
        double t_m_inv=1./t_m;
        double t_e_inv=1./t_e;
        double t_I_inv=1./t_I;

        // update t_m_inv following Liu et al.
        
        double rsq = p->x*p->x + p->y*p->y + p->z*p->z;
        double v_dot_r = p->x*p->vx + p->y*p->vy + p->z*p->vz;
        double q = p->m/p0->m;
        double u = v_dot_r/rsq;
        p->ax += -p->vx*t_m_inv;
        p->ay += -p->vy*t_m_inv;
        p->az += -p->vz*t_m_inv;
        p0->ax += p->vx*t_m_inv*q;
        p0->ay += p->vy*t_m_inv*q;
        p0->az += p->vz*t_m_inv*q;
        p->ax += -2*u*p->x*t_e_inv;
        p->ay += -2*u*p->y*t_e_inv;
        p->az += -2*u*p->z*t_e_inv;
        p0->ax += 2*u*p->x*t_e_inv*q;
        p0->ay += 2*u*p->y*t_e_inv*q;
        p0->az += 2*u*p->z*t_e_inv*q;
    }
}

int main(){
    // initialize simulation
    struct reb_simulation* sim = reb_simulation_create();
    // integrator
    sim->dt             = 1e-2;        // initial timestep.
    sim->integrator     = REB_INTEGRATOR_WHFAST;
    // add particles
    struct reb_particle p; 
    p.m = 1.;    
    reb_simulation_add(sim, p);
    struct reb_particle primary = sim->particles[0];
    double m = initial_planet_mass;
    double a = initial_a;
    double e = 0.0;
    double inc = 0.0;
    double Omega = 0.0;
    double omega = 0.0;
    double f = 0.;
    for (int i=1;i<=initial_n_planet;++i) {
        struct reb_particle plt = reb_particle_from_orbit(sim->G, primary, m, a, e, inc, Omega, omega, f);
        plt.r = mass_to_radius(plt.m);
        plt.hash = i;
        reb_simulation_add(sim,plt);
        // for next planet
        a *= initial_a_ratio;
        f += 1.; // just separate adjacent planets by 1 rad initially
    }
    // add migration force
    sim->additional_forces = migration;
    sim->force_is_velocity_dependent = 1;
    // add collision
    sim->collision         = REB_COLLISION_LINE;
    sim->collision_resolve = collision;
    // initialize output
    FILE *output_file;
    output_file = fopen(output_name, "w");
    fprintf(output_file, "# pid, tid, t, m, r_p, a, e, inc, Omega, omega, f\n");
    fclose(output_file);
    // run simulation and save outputs
    double t_start = clock();
    int tid=0; // time id
    while (sim->t<tmax) {
        // integrate
        reb_simulation_integrate(sim, sim->t + dt_output);
        tid++;
        double t = sim->t;
        // open output file and save data
        output_file = fopen(output_name, "a");
        const int N = sim->N;
        for (int i=1;i<N;i++){
            struct reb_orbit o= reb_orbit_from_particle(sim->G, sim->particles[i], sim->particles[0]);
            int pid = sim->particles[i].hash; // planet id
            double m = sim->particles[i].m;
            double r_p = sim->particles[i].r;
            double a = o.a;
            double e = o.e;
            double inc = o.inc;
            double Omega = o.Omega;
            double omega = o.omega;
            double f = o.f;
            fprintf(output_file, "%d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", pid, tid, t, m, r_p, a, e, inc, Omega, omega, f);
        }
        fclose(output_file);
    }
    double t_end = clock();
    printf("simulation took %f seconds\n",(t_end-t_start)/CLOCKS_PER_SEC);
}