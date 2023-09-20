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

    allow_kmc = 1;
    allow_rejection = 1;
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
   compute energy of site
------------------------------------------------------------------------- */

double AppPottsGSH::site_energy(int i) {
    // here is where we will use our custom distance function.

    // instead of the traditional (1-delta(spin1, spin2))
    // we will use a distance function to represent energy

    double eng = 0.0;

    double* gsh_isite = spin2gsh[spin[i]];
    double* gsh_jsite;

    int nei;
    for (int j = 0; j < numneigh[i]; j++) {
        nei = neighbor[i][j];

        gsh_jsite = spin2gsh[spin[nei]];

        // std::cout << i << ": ";
        // for (int k = 0; k < 10; ++k) {
        //     std::cout << gsh_isite[k] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << neighbor[i][j] << ": ";
        // for (int k = 0; k < 10; ++k) {
        //     std::cout << gsh_jsite[k] << " ";
        // }
        // std::cout << std::endl;

        // exit(0);

        eng += euclideanDistance(gsh_isite, gsh_jsite, n_gsh_coef);
    }

    // if (eng == 0) {
    //     std::cout << "energy: " << eng << std::endl;
    //     exit(0);
    // }

    // exit(0);

    return eng;
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
    int iran = rand() % nspins;
    spin[i] = iran;

    // get new energy
    double efinal = site_energy(i);

    // std::cout << "initial energy: " << einitial << std::endl;
    // std::cout << "new energy: " << efinal << std::endl;
    // exit(0);

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

    // update this masking option. If we can mask, it will run faster.
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

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppPottsGSH::site_event(int i, RandomPark* random) {
    int j, m, value;

    // pick one event from total propensity by accumulating its probability
    // compare prob to threshhold, break when reach it to select event
    // perform event

    double threshhold = random->uniform() * propensity[i2site[i]];
    double efinal;

    int oldstate = spin[i];
    double einitial = site_energy(i);
    double prob = 0.0;
    int nevent = 0;

    for (j = 0; j < numneigh[i]; j++) {
        value = spin[neighbor[i][j]];
        if (value == oldstate) continue;
        for (m = 0; m < nevent; m++)
            if (value == unique[m]) break;
        if (m < nevent) continue;
        unique[nevent++] = value;

        spin[i] = value;
        efinal = site_energy(i);
        if (efinal <= einitial)
            prob += 1.0;
        else if (temperature > 0.0)
            prob += exp((einitial - efinal) * t_inverse);
        if (prob >= threshhold) break;
    }

    // compute propensity changes for self and neighbor sites
    // ignore update of neighbor sites with isite < 0

    int nsites = 0;
    int isite = i2site[i];
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(i);

    for (j = 0; j < numneigh[i]; j++) {
        m = neighbor[i][j];
        isite = i2site[m];
        if (isite < 0) continue;
        sites[nsites++] = isite;
        propensity[isite] = site_propensity(m);
    }

    solve->update(nsites, sites, propensity);
}

/*
-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------- Custom Class Functions -----------------------------
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

double AppPottsGSH::euclideanDistance(const double* array1, const double* array2, int size) {
    double sum = 0.0;

    for (int i = 0; i < size; ++i) {
        double diff = array1[i] - array2[i];
        sum += diff * diff;
    }

    return sqrt(sum);
}