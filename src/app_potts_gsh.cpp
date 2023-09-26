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
    create_arrays();

    allow_kmc = 0;
    allow_rejection = 1;  
    allow_masking = 1;

    // parse arguments for PottsGSH class only, not children

    if (strcmp(style, "potts_gsh") != 0) return;

    if (narg != 2) error->all(FLERR, "Illegal app_style command");

    // get the map which is used to convert spin to euler angles and gsh coef
    char* map_path = arg[1];
    SpinMaps result = read_spin2angle_map(map_path, nspins, n_euler_angles, n_gsh_coef);
    spin2euler = result.spin2euler;
    spin2gsh = result.spin2gsh;

    if (nspins <= 0) error->all(FLERR, "Illegal app_style command");

    // dt_sweep = 1.0 / nspins;  // TODO: Because nspins is so large, this really slows down the simulation
    //                           // where is this taking effect? Why does it slow the program so much?
    //                           // when I make it smaller, the microstructure doesn't evolve the same
    //                           // I think this will most likely be in one of the AppPotts functions
    //                           // propensity, rejection, etc...

    // TODO: this is arbitrary, so we can speed up by doing only neighbors in site_event_rejection
}

/* ---------------------------------------------------------------------- */

AppPottsGSH::~AppPottsGSH() {
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
    // set dt_sweep based on number of neighbors
    dt_sweep = 1.0 / maxneigh;

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
    // TODO: ?
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPottsGSH::site_energy(int i) {

    // TODO: This will be a custom learned function in the future
    //       will most likely take as input two gshs and output energy

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

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   flip to random spin from 1 to nspins
------------------------------------------------------------------------- */

void AppPottsGSH::site_event_rejection(int i, RandomPark* random) {
    
    // save old state and get initial energy
    int oldstate = spin[i];
    double einitial = site_energy(i);

    // get new state, update spin[i], get new energy
    int j = (int) (numneigh[i]*random->uniform());
    spin[i] = spin[neighbor[i][j]];
    double efinal = site_energy(i);

    // accept or reject via Boltzmann criterion
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

/*
-------------------------------------------------------------------------
-------------------------------------------------------------------------
----------------- Custom AppPottsGSH Functions --------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
*/


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
