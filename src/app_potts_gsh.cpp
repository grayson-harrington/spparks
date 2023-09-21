/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "app_potts_gsh.h"

#include <fstream>
#include <iostream>

#include "error.h"
#include "math.h"
#include "random_park.h"
#include "solve.h"
#include "stdlib.h"
#include "string.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsGSH::AppPottsGSH(SPPARKS* spk, int narg, char** arg) : AppPotts(spk, narg, arg) {
    ninteger = 1;
    ndouble = 0;

    delpropensity = 1;
    delevent = 0;

    allow_kmc = 0;        // TODO: do we want kmc or rejection?
    allow_rejection = 1;  // if both kmc and rejection are allowed, what is used?
    allow_masking = 1;
    numrandom = 1;

    create_arrays();
    sites = unique = NULL;

    // parse arguments for PottsGSH class only, not children

    if (strcmp(style, "potts_gsh") != 0) return;

    if (narg != 2) error->all(FLERR, "Illegal app_style command");

    // get the map which is used to convert spin to euler angles and gsh coef
    char* map_path = arg[1];
    SpinMaps result = read_spin2angle_map(map_path, nspins, n_euler_angles, n_gsh_coef);
    spin2euler = result.spin2euler;
    spin2gsh = result.spin2gsh;

    // TODO: create spin2quat

    if (nspins <= 0) error->all(FLERR, "Illegal app_style command");
    dt_sweep = 1.0 / nspins;  // TODO: Because nspins is so large, this really slows down the simulation
                              // where is this taking effect? Why does it slow the program so much?
                              // when I make it smaller, the microstructure doesn't evolve the same
                              // I think this will most likely be in one of the AppPotts functions
                              // propensity, rejection, etc...
}

/* ---------------------------------------------------------------------- */

AppPottsGSH::~AppPottsGSH() {
    delete[] sites;
    delete[] unique;

    // Don't forget to free the memory for the custom arrays being used
    for (int i = 0; i < nspins; ++i) {
        delete[] spin2euler[i];
        delete[] spin2gsh[i];
    }
    delete[] spin2euler;
    delete[] spin2gsh;
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsGSH::grow_app() {
    spin = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsGSH::init_app() {
    delete[] sites;
    delete[] unique;
    sites = new int[1 + maxneigh];
    unique = new int[1 + maxneigh];

    // Initialize the spins based on file inputs.
    // TODO: This will eventually be EBSD data rather than random spins
    for (int i = 0; i < nlocal; i++) {
        spin[i] = rand() % nspins;
    }

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
        if (spin[i] < 0 || spin[i] > nspins - 1) flag = 1;
    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
    if (flagall) error->all(FLERR, "One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
    user defined optional parameters
 ------------------------------------------------------------------------- */
void AppPottsGSH::input_app(char* command, int narg, char** arg) {
    if (strcmp(command, "site_energy_type") == 0) {
        if (narg == 0) error->all(FLERR, "Illegal humph_n command\n");

        site_energy_type = atoi(arg[0]);

        double input;
        char* units;

        switch (site_energy_type) {
            case 1:
                // spin count
                // no additional parameters required
                break;
            case 2:
                // euler misorientation
                // theta_m required
                if (narg != 3) error->all(FLERR, "Illegal type 2 parameter count\n");

                input = atof(arg[1]);
                units = arg[2];

                if (strcmp(units, "d") == 0 || strcmp(units, "deg") == 0 || strcmp(units, "degree") == 0 || strcmp(units, "degrees") == 0) {
                    theta_m = input / 180.0 * M_PI;
                } else if (strcmp(units, "r") == 0 || strcmp(units, "rad") == 0 || strcmp(units, "radian") == 0 || strcmp(units, "radians") == 0) {
                    theta_m = input;
                } else {
                    error->all(FLERR, "Illegal angle unit specifier\n");
                }

                if (theta_m < 0) {
                    error->all(FLERR, "theta_m input must be positive\n");
                } else if (theta_m > M_PI) {
                    error->all(FLERR, "theta_m must be between 0-180 degrees or 0-PI radians\n");
                }

                if (logfile)
                    fprintf(logfile, "  site_energy_type set to %d\n", site_energy_type);
                if (logfile)
                    fprintf(logfile, "  theta_m set to %f radians\n", theta_m);
                if (screen && me == 0)
                    fprintf(screen, "  theta_m set to %f radians\n", theta_m);

                break;
            case 3:
                // gsh distance
                // no additional parameters required, yet
                break;
            default:
                error->all(FLERR, "Invalid site_energy_type\n");
        }
    }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPottsGSH::site_energy(int i) {
    // here is where we will use our custom distance function.

    switch (site_energy_type) {
        case 1:
            return site_energy_spin(i);
        case 2:
            return site_energy_read_shockley(i);
        case 3:
            return site_energy_gsh(i);
        default:
            error->all(FLERR, "Invalid site_energy_type\n");
            return INFINITY;
    }

    
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   flip to random spin from 1 to nspins
------------------------------------------------------------------------- */

void AppPottsGSH::site_event_rejection(int i, RandomPark* random) {
    // save old state
    int oldstate = spin[i];

    // get initial energy
    double einitial = site_energy(i);

    // get new state and update spin[i]
    // TODO: should this be completely random, or from neighbors?
    int iran = rand() % nspins;
    spin[i] = iran;

    // get new energy
    double efinal = site_energy(i);

    // accept or reject via Boltzmann criterion
    // null bin extends to nspins

    if (efinal <= einitial) {
    } else if (temperature == 0.0) {
        spin[i] = oldstate;
    } else if (random->uniform() > exp((einitial - efinal) * t_inverse)) {
        spin[i] = oldstate;
    }

    if (spin[i] != oldstate) naccept++;

    // set mask if site could not have changed
    // if site changed, unset mask of sites with affected propensity
    // OK to change mask of ghost sites since never used
    if (Lmask) {
        if (einitial == 0) mask[i] = 1;
        if (spin[i] != oldstate)
            for (int j = 0; j < numneigh[i]; j++)
                mask[neighbor[i][j]] = 0;
    }
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

// TODO: Unchanged so far
double AppPottsGSH::site_propensity(int i) {
    // events = spin flips to neighboring site different than self
    // disallow wild flips = flips to value different than all neighs

    int j, m, value;
    int nevent = 0;

    for (j = 0; j < numneigh[i]; j++) {
        value = spin[neighbor[i][j]];
        if (value == spin[i]) continue;
        for (m = 0; m < nevent; m++)
            if (value == unique[m]) break;
        if (m < nevent) continue;
        unique[nevent++] = value;
    }

    // for each flip:
    // compute energy difference between initial and final state
    // if downhill or no energy change, propensity = 1
    // if uphill energy change, propensity = Boltzmann factor

    int oldstate = spin[i];
    double einitial = site_energy(i);
    double efinal;
    double prob = 0.0;

    for (m = 0; m < nevent; m++) {
        spin[i] = unique[m];
        efinal = site_energy(i);
        if (efinal <= einitial)
            prob += 1.0;
        else if (temperature > 0.0)
            prob += exp((einitial - efinal) * t_inverse);
    }

    spin[i] = oldstate;
    return prob;
}

/*
-------------------------------------------------------------------------
-------------------------------------------------------------------------
----------------- Custom AppPottsGSH Functions --------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
*/

double AppPottsGSH::site_energy_spin(int i) {
    int isite = spin[i];
    int eng = 0;
    for (int j = 0; j < numneigh[i]; j++)
        if (isite != spin[neighbor[i][j]]) eng++;
    return (double)eng;
}

double AppPottsGSH::site_energy_read_shockley(int i) {
    // double quat_i[4] = {quat_a[i], quat_b[i], quat_c[i], quat_d[i]};

    // double energy = 0.0;
    // int nei;
    // for (int j = 0; j < numneigh[i]; j++) {
    //     nei = neighbor[i][j];
    //     double quat_nei[4] = {quat_a[nei], quat_b[nei], quat_c[nei], quat_d[nei]};
    //     double miso = angle_between(quat_i, quat_nei);
    //     energy += read_shockley(miso);
    // }
    // return (double)energy;

    // TODO: have this change based on site_energy_type

    return 0;
}

double AppPottsGSH::site_energy_gsh(int i) {
    // TODO: This will be a custom learned function in the future
    double gsh_dist;
    double eng = 0.0;

    double* gsh_isite = spin2gsh[spin[i]];
    double* gsh_jsite;

    int nei;
    for (int j = 0; j < numneigh[i]; j++) {
        nei = neighbor[i][j];
        gsh_jsite = spin2gsh[spin[nei]];
        eng += euclideanDistance(gsh_isite, gsh_jsite, n_gsh_coef);
    }

    return eng;
}

// Read-Shockley misorientation energy
//  https://www.desmos.com/calculator/gleyqfqseq
double AppPottsGSH::read_shockley(const double theta) {
    double energy = 0;

    if (theta >= theta_m) {
        energy = 1;
    } else if (theta == 0) {
        energy = 0;
    } else {
        double ratio = theta / theta_m;
        energy = ratio * (1 - log(ratio));
    }

    return energy;
}

/*
-------------------------------------------------------------------------
-------------------------------------------------------------------------
------------------------ Custom Functions -------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
*/

SpinMaps AppPottsGSH::read_spin2angle_map(const char* filePath, int& n_lines, int& n_eul, int& n_gsh) {
    std::ifstream inputFile(filePath);

    if (!inputFile) {
        std::cerr << "Error: Cannot open the file." << std::endl;
        exit(1);  // Exit the program with an error code (e.g., 1)
    }

    inputFile >> n_lines >> n_eul >> n_gsh;

    // Create dynamic 2D arrays for Euler angles and GSH coefficients
    double** euler_angles = new double*[n_lines];
    double** gsh_coefs = new double*[n_lines];

    for (int i = 0; i < n_lines; ++i) {
        int id;
        inputFile >> id;

        euler_angles[id] = new double[n_eul];
        gsh_coefs[id] = new double[n_gsh];

        // Read Euler angles
        for (int j = 0; j < n_eul; ++j) {
            inputFile >> euler_angles[id][j];
        }

        // Read GSH coefficients
        for (int j = 0; j < n_gsh; ++j) {
            inputFile >> gsh_coefs[id][j];
        }
    }

    inputFile.close();

    SpinMaps data;
    data.spin2euler = euler_angles;
    data.spin2gsh = gsh_coefs;

    return data;
}

double AppPottsGSH::euclideanDistance(const double* array1, const double* array2, const int size) {
    double sum = 0.0;

    for (int i = 0; i < size; ++i) {
        double diff = array1[i] - array2[i];
        sum += diff * diff;
    }

    return sqrt(sum);
}

void AppPottsGSH::euler2quaternion(const double euler_in[3], double quat_out[4]) {
    double alpha = euler_in[0];
    double beta = euler_in[1];
    double gamma = euler_in[2];

    double sigma = 0.5 * (alpha + gamma);
    double delta = 0.5 * (alpha - gamma);

    double c = cos(beta / 2);
    double s = sin(beta / 2);

    quat_out[0] = c * cos(sigma);
    quat_out[1] = -s * cos(delta);
    quat_out[2] = -s * sin(delta);
    quat_out[3] = -c * sin(sigma);

    if (quat_out[0] < 0) {
        quat_out[0] *= -1;
        quat_out[1] *= -1;
        quat_out[2] *= -1;
        quat_out[3] *= -1;
    }
}

double AppPottsGSH::angle_between(const double quat1[4], const double quat2[4]) {
    double dot = quat1[0] * quat2[0] + quat1[1] * quat2[1] + quat1[2] * quat2[2] + quat1[3] * quat2[3];
    double angle = acos(2 * pow(dot, 2) - 1);
    if (isnan(angle)) {
        angle = 0;
    }
    return angle;
}