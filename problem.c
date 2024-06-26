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
#include <string.h>
#include <math.h>
#include <time.h>
#include "rebound.h"

const char input_name[] = "sim.in";
const char output_name[] = "sim.out";

// these variables can be updated by the input file
double tmax_yr = 1.e5;
double dt_output_yr = 1.e2;
double initial_planet_mass = 3e-6; // initial mass of innermost planet
double initial_planet_mass_ratio = 1.; // initial planet mass ratios: m_{i+1} / m_i
int initial_n_planet = 4;
double initial_a = 1.; // initial semimajor axis of innermost planet
double initial_a_ratio = 1.2; // initial semimajor axes ratios: a_{i+1} / a_i

// disk properties
// these numbers follow Liu et al.
const double Mdot_0_Msunyr = 1.e-8;
const double alpha = 0.01;
const double tau_d_yr = 1.e5;
const double B_kG = 1.;
const double h0 = 0.025;

const double Msun_in_g = 1.9891e33;
const double au_in_cm = 1.49597871e13;
const double G = 1.;

const double M = 1.; // stellar mass

double mass_to_radius(const double m){
    // CAVEAT: we just scale by earth values for now...
    double r = pow(m/3.e-6, 1./3.)*4.26e-04;
    return r;
}

int collision(struct reb_simulation* const sim, struct reb_collision c){
    printf("\rcollision happened at t = %.2f                                            \n", sim->t);
    struct reb_particle* particles = sim->particles;
    double m = particles[c.p1].m+particles[c.p2].m;
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
        double beta_T = 1; // temperature slope, T ~ r^-beta_T
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
        // CAVEAT: e and I dependence are not implemented

        // update t_m_inv following Liu et al.
        // normalized torque = torque / m r^2 Omega^2
        // 2 * normalized torque * Omega = dadt/a = 2/tm
        // 1/tm = normalized torque * Omega
        double qp = m/M;
        double qd = Sigma*a*a/M;
        double normalized_torque_2s = (-(2.5+0.5*beta_T+0.1*beta)+1.4*beta_T+1.1*(1.5+beta))*qd*qp/(h*h);
        double normalized_torque_1s = -0.65*qd*qp/(h*h*h) + 2.46*qd*sqrt(qp/(h*h*h))*exp(-e/(.5*h+0.01));
        double x_hs = 1.7*sqrt(qp/h)*a;
        double normalized_torque = 0.;
        if (a>r_c) {
            double f = exp(-(a-r_c)/x_hs);
            normalized_torque = f*normalized_torque_1s + (1-f)*normalized_torque_2s;
        }
        t_m_inv = normalized_torque * Omega;
        
        // apply force to planet
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
        p->ax += -2.*u*p->x*t_e_inv;
        p->ay += -2.*u*p->y*t_e_inv;
        p->az += -2.*u*p->z*t_e_inv;
        p0->ax += 2.*u*p->x*t_e_inv*q;
        p0->ay += 2.*u*p->y*t_e_inv*q;
        p0->az += 2.*u*p->z*t_e_inv*q;
        // CAVEAT: no inclination damping now
    }
}

int main(int argc, char *argv[]){
    // parse input file
    FILE *input_file = fopen(input_name, "r");
    if (input_file == NULL) {
        perror("Could not open input file! Check if sim.in exists.");
        exit(EXIT_FAILURE);
    }
    char line[256];
    char key[256];
    double value;
    while (fgets(line, sizeof(line), input_file)) {
        if (sscanf(line, "%s = %lf", key, &value) == 2 && strcmp(key, "tmax_yr") == 0) tmax_yr = value;
        if (sscanf(line, "%s = %lf", key, &value) == 2 && strcmp(key, "dt_output_yr") == 0) dt_output_yr = value;
        if (sscanf(line, "%s = %lf", key, &value) == 2 && strcmp(key, "initial_planet_mass") == 0) initial_planet_mass = value;
        if (sscanf(line, "%s = %lf", key, &value) == 2 && strcmp(key, "initial_planet_mass_ratio") == 0) initial_planet_mass_ratio = value;
        if (sscanf(line, "%s = %lf", key, &value) == 2 && strcmp(key, "initial_n_planet") == 0) initial_n_planet = value;
        if (sscanf(line, "%s = %lf", key, &value) == 2 && strcmp(key, "initial_a_ratio") == 0) initial_a_ratio = value;
        if (sscanf(line, "%s = %lf", key, &value) == 2 && strcmp(key, "initial_a") == 0) initial_a = value;
    }
    fclose(input_file);
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
        m *= initial_planet_mass_ratio;
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
    double tmax = tmax_yr*2.*M_PI;
    double dt_output = dt_output_yr*2.*M_PI;
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
        // print progress
        double t_current = clock();
        printf("\rt = %.2f, tmax = %.2f, progress = %.2f%%, %.2f seconds taken", sim->t, tmax, sim->t/tmax*100., (t_current-t_start)/CLOCKS_PER_SEC);
        fflush(stdout);
    }
    double t_end = clock();
    printf("\nsimulation took %f seconds\n",(t_end-t_start)/CLOCKS_PER_SEC);
}