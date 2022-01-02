/**
 * @file lennard.c
 * @brief 2D simulation of particles interacting with Lennard-Jones potential and boundary conditions.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define box_l 20. // Square's side
#define nb_particles 180
#define rho (box_a / nb_particles)  // 2D density of particles
#define Rc_2 6.25                   // Truncation radius
#define Ep_Rc -0.014512279230439821 // Potential energy as truncation radius

typedef struct Particle Particle;
struct Particle
{
    double x;
    double y;
    double vx;
    double vy;
};

typedef struct SumForces SumForces;
struct SumForces
{
    double sf_x;
    double sf_y;
};

/**
 * @brief Compute the total energy of the system (kinetic energies + potential energies)
 * 
 * @param N size of particles array
 * @param particles array containing the particles in the system
 * 
 * @return double, Total energy
 */
double total_energy(int N, Particle particles[])
{
    double E, Ec, Ep, dx, dy, r_2i, r_6i;
    E = Ec = Ep = 0;
    for (int i = 0; i < N; i++)
    {
        Ec += 0.5 * (particles[i].vx * particles[i].vx +
                     particles[i].vy * particles[i].vy);
        for (int j = i + 1; j < N; j++)
        {
            dx = particles[j].x - particles[i].x;
            if (dx > box_l / 2)
                dx -= box_l;
            else if (dx < -box_l / 2)
                dx += box_l;

            dy = particles[j].y - particles[i].y;
            if (dy > box_l / 2)
                dy -= box_l;
            else if (dy < -box_l / 2)
                dy += box_l;

            r_2i = 1 / (dx * dx + dy * dy);

            if (r_2i > 1 / Rc_2)
            {
                r_6i = r_2i * r_2i * r_2i;
                Ep += 4 * (r_6i * r_6i) - 4 * r_6i - Ep_Rc;
            }
        }
    }
    E = Ec + Ep;
    return E;
}

/**
 * @brief  Compute lennard-jones forces between 2 particles in 2D space
 * 
 * @param dx distance on axis x between the particles
 * @param dy distance on axis y between the particles
 * 
 * @return double, Force between the particles
 */
static inline double lennard_jones_forces(double dx, double dy)
{
    double r_2i = 1 / (dx * dx + dy * dy); // Squared radius
    if (r_2i < 1 / Rc_2)
        return 0;
    double r_6i = r_2i * r_2i * r_2i;
    double r_8i = r_6i * r_2i;
    return 48. * dx * r_8i * r_6i - 24. * dx * r_8i;
}

/**
 * @brief Velocity verlet algorithm
 * 
 * @param dt Step precision of time
 * @param tmax Computation time
 * @param N Row size of particles array
 * @param particles Particles container
 */
void v_verlet(double dt, double tmax, int N, Particle particles[])
{
    // Forces' sum at time t and t + 1
    double force_x, force_y; // Value of force x/y between 2 particles
    SumForces *forces;       // Forces at times t
    SumForces *forces_1;     // Forces at times t + dt

    // Initialising array containing forces of size N
    forces = calloc(N, sizeof(SumForces));
    forces_1 = calloc(N, sizeof(SumForces));

    double dx, dy; // Computation, distance x/y between 2 particles
    double t = 0;  // Time
    unsigned int nt = 0;
    double E_0 = total_energy(N, particles);
    // Compute forces at t = 0
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            dx = particles[i].x - particles[j].x;
            // Periodic conditions
            if (dx > box_l / 2)
                dx -= box_l;
            else if (dx < -box_l / 2)
                dx += box_l;

            dy = particles[i].y - particles[j].y;
            if (dy > box_l / 2)
                dy -= box_l;
            else if (dy < -box_l / 2)
                dy += box_l;

            force_x = lennard_jones_forces(dx, dy);
            force_y = lennard_jones_forces(dy, dx);
            forces[i].sf_x += force_x;
            forces[j].sf_x -= force_x;
            forces[i].sf_y += force_y;
            forces[j].sf_y -= force_y;
        }
    }
    while (t < tmax)
    {
        // Moving particles with forces at time t computed at t - dt
        for (int i = 0; i < N; i++)
        {
            particles[i].x += particles[i].vx * dt + 0.5 * dt * dt * forces[i].sf_x;
            particles[i].y += particles[i].vy * dt + 0.5 * dt * dt * forces[i].sf_y;
            if (particles[i].x < 0)
                particles[i].x += box_l;

            if (particles[i].x > box_l)
                particles[i].x -= box_l;

            if (particles[i].y < 0)
                particles[i].y += box_l;

            if (particles[i].y > box_l)
                particles[i].y -= box_l;
            
            if (nt % 500 == 0)
			{
				printf("%lf;%lf;", particles[i].x, particles[i].y);
			}
        }

        for (int i = 0; i < N; i++)
        {
            // Compute forces at time t + dt
            for (int j = i + 1; j < N; j++)
            {
                dx = particles[i].x - particles[j].x;
                // Periodic conditions
                if (dx > box_l / 2)
                    dx -= box_l;
                else if (dx < -box_l / 2)
                    dx += box_l;

                dy = particles[i].y - particles[j].y;
                if (dy > box_l / 2)
                    dy -= box_l;
                else if (dy < -box_l / 2)
                    dy += box_l;

                force_x = lennard_jones_forces(dx, dy);
                force_y = lennard_jones_forces(dy, dx);
                forces_1[i].sf_x += force_x;
                forces_1[j].sf_x -= force_x;
                forces_1[i].sf_y += force_y;
                forces_1[j].sf_y -= force_y;
            }

            // Computing and adding velocity
            particles[i].vx += 0.5 * (forces_1[i].sf_x + forces[i].sf_x) * dt;
            particles[i].vy += 0.5 * (forces_1[i].sf_y + forces[i].sf_y) * dt;

            // t + dt -> t
            forces[i].sf_x = forces_1[i].sf_x;
            forces[i].sf_y = forces_1[i].sf_y;
            // Reseting force t + dt
            forces_1[i].sf_x = 0;
            forces_1[i].sf_y = 0;
        }
        t += dt;
        if (nt % 500 == 0)
		{
			printf("%lf\n", t);
		}
        nt++;
    }
    free(forces);
    free(forces_1);
}

int main()
{
    double dt = 3e-4; // Precision
    double tmax;
    double N_line = (ceil(sqrt(nb_particles))); // Nb particles per line
    char u;

    // N-particles in 2D space, 4 columns (x, y, Vx, Vy) for each particle
    Particle *particles;
    particles = calloc(nb_particles, sizeof(Particle));
    tmax = dt * 1000000;
    // Disposing
    u = -1;
    for (int i = 0; i < nb_particles; i++)
    {
        if (i % (int)N_line == 0)
        {
            u++;
        }
        particles[i].x = (double)(i % (int)N_line) * (box_l / N_line);
        particles[i].y = (double)u * (box_l / N_line);
    }
    v_verlet(dt, tmax, nb_particles, particles);
    free(particles);

    return 0;
}
