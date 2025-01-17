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

#include "app_potts_neighonly.h"

#include "error.h"
#include "math.h"
#include "random_park.h"
#include "stdlib.h"
#include "string.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsNeighOnly::AppPottsNeighOnly(SPPARKS *spk, int narg, char **arg) : AppPotts(spk, narg, arg) {
    // parse arguments for PottsNeighOnly class only, not children

    if (strcmp(style, "potts/neighonly") != 0) return;

    if (narg != 2) error->all(FLERR, "Illegal app_style command");

    nspins = atoi(arg[1]);
    if (nspins <= 0) error->all(FLERR, "Illegal app_style command");
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsNeighOnly::init_app() {
    delete[] sites;
    delete[] unique;
    sites = new int[1 + maxneigh];
    unique = new int[1 + maxneigh];

    dt_sweep = 1.0 / maxneigh;

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
        if (spin[i] < 1 || spin[i] > nspins) flag = 1;
    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
    if (flagall) error->all(FLERR, "One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */

void AppPottsNeighOnly::site_event_rejection(int i, RandomPark *random) {
    int oldstate = spin[i];
    double einitial = site_energy(i);

    // events = spin flips to neighboring site different than self

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

    if (nevent == 0) return;
    int iran = (int)(nevent * random->uniform());
    if (iran >= nevent) iran = nevent - 1;
    spin[i] = unique[iran];
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
        if (einitial < 0.5 * numneigh[i]) mask[i] = 1;
        if (spin[i] != oldstate)
            for (int j = 0; j < numneigh[i]; j++)
                mask[neighbor[i][j]] = 0;
    }
}
