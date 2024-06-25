import numpy as np
import rebound as rb

"""
units:
use rebound units (Msun - au - 1/2pi yr) unless otherwise noted
"""
Msun_in_g = 1.9891e33
au_in_cm = 1.49597871e13
G = 1


def migration_rate_CN(a,e,I,M,m,Sigma,h,beta):
    """
    migration rate following Cresswell & Nelson
    input:
        a,e,I: planet orbital elements
        M: stellar mass
        m: planet mass
        Sigma: surface density
        h: aspect ratio (H/R)
        beta: surface density slope
    output:
        inverse of t_m, t_e, t_I
    """
    Omega = np.sqrt(G*M/a**3)
    t_wave = M/m * M/(Sigma*a**2) * h**4 / Omega
    # TODO: add the e,i dependences
    t_e = t_wave/0.78
    t_I = t_wave/0.544
    t_m = 2*t_wave/(2.7+1.1*beta)*h**-2
    return 1/t_m, 1/t_e, 1/t_I


def migration_rate_Liu():
    """
    migration rate folloing Lui et al.
    """
    # TODO: implement this
    t_m_inv = 0.
    return t_m_inv


def disk(r,t,Mdot_0_Msunyr,alpha,tau_d_yr,B_kG,h0):
    """
    disk surface density
    input:
        Mdot_0_Msunyr: initial Mddot in Msun/yr
        alpha: viscosity
        tau_d_yr: disk dispersal time in yr
        B_kG: stellar magnetic field in kilo Gauss
        h0: aspect ratio (H/R)
    """
    tau_d = tau_d_yr*2*np.pi
    Mdot_Msunyr = Mdot_0_Msunyr * np.exp(-t/tau_d)
    Sigma_gcm2 = 75 * (Mdot_Msunyr/1e-9) * (alpha/1e-2)**-1 * (r/0.1)**-0.5
    Sigma = Sigma_gcm2 / (Msun_in_g/au_in_cm**2)
    r_c = 0.1 * (Mdot_Msunyr/1e-9)**(-2/7) * (B_kG)**(4/7) # truncation radius
    if (r<r_c): Sigma = 1e-40
    beta = 0.5 # surface density slope, Sigma ~ r^-beta
    h = h0
    return Sigma, h, beta


def simulation(
    Mdot_0_Msunyr=1e-9,alpha=1e-2,tau_d_yr=1e4,B_kG=1,h0=0.025, # TODO revisit fiducial values
    duration_yr = 1e6, n_output = 1000,
    ):
    """
    output:
    data: dictionary that stores time evolution of the system
    for example, data[a] is a (n_t, n_p) array where n_p = number of planets and n_t = number of output times
    """

    # 1. set up the simulation
    sim = rb.Simulation()
    ps = sim.particles # setting ps as global variable to be called in force routine - this is intentional
    # 1.1 initialize simulation and add force
    def migraton(reb_sim):
        """
        force perscription used for rebound
        """
        for i in range(1,len(ps)):
            # CAVEAT: we assume coplanar here - no incination damping
            a = ps[i].a
            e = ps[i].e
            I = ps[i].inc
            m = ps[i].m
            M = ps[0].m
            t = sim.t
            Sigma, h, beta = disk(a,t,Mdot_0_Msunyr,alpha,tau_d_yr,B_kG,h0) # use r=a for disk properties
            t_m_inv, t_e_inv, t_i_inv = migration_rate_CN(a,e,I,M,m,Sigma,h,beta)
            v_dot_r = (ps[i].vx*ps[i].x) + (ps[i].vy*ps[i].y) + (ps[i].vz*ps[i].z)
            r = np.sqrt(((ps[i].x)**2)+((ps[i].y)**2)+((ps[i].z)**2))
            u = v_dot_r / r**2
            # force on planet
            ps[i].ax -= 2*u*ps[i].x*t_e_inv
            ps[i].ay -= 2*u*ps[i].y*t_e_inv
            ps[i].az -= 2*u*ps[i].z*t_e_inv
            ps[i].ax -= ps[i].vx*t_m_inv
            ps[i].ay -= ps[i].vy*t_m_inv
            ps[i].az -= ps[i].vz*t_m_inv
            # force on star - keep center of mass fixed
            ps[0].ax += 2*u*ps[i].x*t_e_inv*(m/M)
            ps[0].ay += 2*u*ps[i].y*t_e_inv*(m/M)
            ps[0].az += 2*u*ps[i].z*t_e_inv*(m/M)
            ps[0].ax += ps[i].vx*t_m_inv*(m/M)
            ps[0].ay += ps[i].vy*t_m_inv*(m/M)
            ps[0].az += ps[i].vz*t_m_inv*(m/M)
        return
    sim.additional_forces = migraton
    sim.force_is_velocity_dependent = 1
    # 1.2 collision routine
    # TODO
    # 1.3 initial conditions & time integrator
    # star
    sim.add(m=1)
    # planets
    sim.add(m = 3e-5 ,r=5.2e-5, a = 1, e = 0.1, f = np.pi)
    sim.add(m = 3e-5 ,r=5.2e-5, a = 2, e = 0.1, f = 0)
    # move frame
    sim.move_to_com()
    # integrator setup
    sim.integrator = "whfast"
    sim.dt = 1e-2

    # 2. run the simulation and save data
    # initialize data
    data_key_sim = ["t","t_yr","n_p"]
    data_key_particle = [
        "hash", # unique id for each particle; string in rebound
        "m","r", # mass and size of particle
        "a","e","inc","f","pomega","Omega" # orbital elements
    ]
    data = {}
    for k in data_key_sim:
        data[k] = np.zeros(n_output+1)
    for k in data_key_particle:
        data[k] = np.zeros((n_output+1, len(ps)-1))
    data["hash"] = np.zeros((n_output+1, len(ps)-1), dtype='str')
    duration = duration_yr*2*np.pi
    t = np.linspace(0,duration,n_output+1)
    for i, t1 in enumerate(t):
        if t1>0: sim.integrate(t1)
        # save data
        data["t"][i] = sim.t
        data["t_yr"][i] = sim.t*2*np.pi
        data["n_p"][i] = len(ps)
        for ip,p in enumerate(ps[1:]):
            for k in data_key_particle:
                data[k][i,ip] = getattr(p, k)
    return data