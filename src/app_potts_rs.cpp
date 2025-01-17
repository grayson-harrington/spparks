/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
   ----------------------------------------------------------------------
   Read-Shockley implementation into SPPARKS AppPotts using quaternions.
   Developed by Efrain Hernandez-Rivera
   US Army Research Laboratory
   Version: v2.4 (06Jul2017)
   --
   Use of maps for computationally efficient misorientation calculation
   Edited for clarity
            v3.0 (28Apr2018)
   Fixed bug with ghost/maps
            v3.1 (10May2018)
   Fixed bug with creation of new boudaries if not in map
            v3.2 (28Jun2018)
------------------------------------------------------------------------- */

#include "app_potts_rs.h"

#include "comm_lattice.h"
#include "domain.h"
#include "error.h"
#include "math.h"
#include "random_park.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

using namespace SPPARKS_NS;

#define MY_PI 3.14159265358979323846   // pi
#define MY_2PI 6.28318530717958647692  // 2pi

/* ---------------------------------------------------------------------- */

AppPottsRS::AppPottsRS(SPPARKS *spk, int narg, char **arg) : AppPotts(spk, narg, arg) {
    // allow for definition of double arrays
    // 1 per Euler angle
    ninteger = 1;
    ndouble = 3;

    // add the extra arrays
    recreate_arrays();

    // only error check for this class, not derived classes
    if (strcmp(arg[0], "potts/rs") == 0 && narg < 2)
        error->all(FLERR, "Illegal app_style command");

    // cutoff misorientation angle between low and high angle GB
    // probably on the high end (ie typically 15 deg)
    thetam = 25.0 / 180.0 * MY_PI;

    // interaction (interfacial) energy
    Jij = 1.0;

    // Mobility parameters
    nmob = 4.0;
    bmob = 5.0;

    // Symmetry operator
    Osym = 24;
    if (narg == 3)
        Osym = atoi(arg[2]);
}

/* ----------------------------------------------------------------------
   destructor
------------------------------------------------------------------------- */
AppPottsRS::~AppPottsRS() {
    // free up memory from quaternion symmetry operator
    for (int i = 0; i < Osym; i++)
        delete[] symquat[i];
    delete[] symquat;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */
void AppPottsRS::init_app() {
    delete[] sites;
    delete[] unique;
    sites = new int[1 + maxneigh];
    unique = new int[1 + maxneigh];

    dt_sweep = 1.0 / maxneigh;

    int flag = 0;
    // Check that all three angles are within range
    for (int i = 0; i < nlocal; i++) {
        // Randonly distribute Euler angles
        phi1[i] = MY_2PI * phi1[i];
        Phi[i] = acos(2 * Phi[i] - 1);
        phi2[i] = MY_2PI * phi2[i];
        if (phi1[i] < 0 || phi1[i] >= MY_2PI) flag = 1;
        if (phi2[i] < 0 || phi2[i] >= MY_2PI) flag = 1;
        if (Phi[i] < 0 || Phi[i] >= MY_PI) flag = 1;
    }

    // Initialize symmetry operator as quaternion vectors
    // Osym = 24 (cubic), 12 (hexagonal)
    symmat(&symquat);

    comm->all();

    int spi, spn, nei;
    double qi[4], qj[4];

    for (int i = 0; i < nlocal + nghost; i++) {
        spi = spin[i];
        euler2quat(i, qi);

        for (int n = 0; n < numneigh[i]; n++) {
            nei = neighbor[i][n];
            spn = spin[nei];

            // order by min/max to avoid duplicate pairs
            int smin = MIN(spi, spn);
            int smax = MAX(spi, spn);

            std::pair<int, int> spins = std::make_pair(smin, smax);

            if (spi != spn && misos.count(spins) == 0) {
                euler2quat(nei, qj);

                // insert misorientation angle between spin ID pair (thetar)
                misos[spins] = quaternions(qi, qj) / thetam;
            }
        }
    }

    if (logfile)
        fprintf(logfile, "  Pairs misorientation map created\n");
    if (screen && me == 0)
        fprintf(screen, "  Pairs misorientation map created\n");

    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
    if (flagall) error->all(FLERR, "One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
 set site value ptrs each time iarray/darray are reallocated
 ------------------------------------------------------------------------- */

void AppPottsRS::grow_app() {
    // set pointers
    // to define these, use command create_sites box iN and set iN
    spin = iarray[0];
    phi1 = darray[0];
    Phi = darray[1];
    phi2 = darray[2];
}

/* ----------------------------------------------------------------------
 *    user defined optional parameters
 *    ---------------------------------------------------------------------- */

void AppPottsRS::input_app(char *command, int narg, char **arg) {
    if (narg < 1) {
        error->all(FLERR, "Invalid command for app_style");
    }

    // Redefine mobility parameters (n,b)
    if (strcmp(command, "mobility") == 0) {
        if (narg != 2)
            error->all(FLERR,
                       "Illegal mobility flag requires two arguments,"
                       " parameter-flag and parameter-value, (e.g. mobility expo"
                       " 3.0)\n");

        if (strcmp(arg[0], "expo") == 0) {
            nmob = atof(arg[1]);

            if (logfile)
                fprintf(logfile, "  Mobility exponent reset to %g\n", nmob);
            if (screen && me == 0)
                fprintf(screen, "  Mobility exponent reset to %g\n", nmob);

        } else if (strcmp(arg[0], "scale") == 0) {
            bmob = atof(arg[1]);

            if (logfile)
                fprintf(logfile, "  Mobility scaling reset to %g\n", bmob);
        } else
            error->all(FLERR, "Mobility parameter not recognized\n");
    }
    // Cutoff angle for Read-Shockley
    else if (strcmp(command, "cutoff") == 0) {
        if (narg < 1)
            error->all(FLERR, "Illegal cutoff angle command\n");
        thetam = fabs(atof(arg[0])) / 180.0 * MY_PI;
        if (thetam > MY_2PI)
            error->all(FLERR,
                       "Cutoff angle must be defined in terms of degrees "
                       "(0,360)\n");

        if (logfile)
            fprintf(logfile, "  Low-to-high angle cutoff reset to %s deg\n", arg[0]);
        if (screen && me == 0)
            fprintf(screen, "  Low-to-high angle cutoff reset to %s deg\n", arg[0]);
    }
    // Potts interfacial energy scaler
    else if (strcmp(command, "energy_scaling") == 0) {
        if (narg < 1)
            error->all(FLERR, "Illegal scaling energy command\n");
        Jij = atof(arg[0]);
        if (Jij < 0)
            error->all(FLERR, "Illegal energy value (>0)\n");

        if (logfile)
            fprintf(logfile, "  PMC energy scaling by %g.\n", Jij);
        if (screen && me == 0)
            fprintf(screen, "  PMC energy scaling by %g.\n", Jij);
    } else
        error->all(FLERR, "Input command not recognized by app\n");
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPottsRS::site_energy(int i) {
    int nei;
    double eng = 0.0, qi[4], qj[4], thetar;

    euler2quat(i, qi);

    for (int j = 0; j < numneigh[i]; j++) {
        nei = neighbor[i][j];
        if (spin[i] == spin[nei]) continue;

        int smin = MIN(spin[i], spin[nei]);
        int smax = MAX(spin[i], spin[nei]);

        std::pair<int, int> spins = std::make_pair(smin, smax);

        // ratio of theta/theta_m
        if (misos.count(spins) == 1)
            thetar = misos[spins];
        else {
            euler2quat(nei, qj);
            thetar = quaternions(qi, qj) / thetam;
            misos[spins] = thetar;
        }

        if (thetar >= 1.0 || thetam < 1e-8)
            eng += 1;
        else if (thetar > 0.0)
            eng += thetar * (1.0 - log(thetar));
    }

    return Jij * eng;
}

/* ----------------------------------------------------------------------
  Convert symmetry matrix to quaternion form
------------------------------------------------------------------------- */

void AppPottsRS::mat2quat(const double O[3][3], double q[4]) {
    double q4 = 0;
    if ((1 + O[0][0] + O[1][1] + O[2][2]) > 0) {
        q4 = sqrt(1 + O[0][0] + O[1][1] + O[2][2]) / 2;
        q[0] = q4;
        q[1] = (O[2][1] - O[1][2]) / (4 * q4);
        q[2] = (O[0][2] - O[2][0]) / (4 * q4);
        q[3] = (O[1][0] - O[0][1]) / (4 * q4);
    } else if ((1 + O[0][0] - O[1][1] - O[2][2]) > 0) {
        q4 = sqrt(1 + O[0][0] - O[1][1] - O[2][2]) / 2;
        q[0] = (O[2][1] - O[1][2]) / (4 * q4);
        q[1] = q4;
        q[2] = (O[1][0] + O[0][1]) / (4 * q4);
        q[3] = (O[0][2] + O[2][0]) / (4 * q4);
    } else if ((1 - O[0][0] + O[1][1] - O[2][2]) > 0) {
        q4 = sqrt(1 - O[0][0] + O[1][1] - O[2][2]) / 2;
        q[0] = (O[0][2] - O[2][0]) / (4 * q4);
        q[1] = (O[1][0] + O[0][1]) / (4 * q4);
        q[2] = q4;
        q[3] = (O[2][1] + O[1][2]) / (4 * q4);
    } else if ((1 - O[0][0] - O[1][1] + O[2][2]) > 0) {
        q4 = sqrt(1 - O[0][0] - O[1][1] + O[2][2]) / 2;
        q[0] = (O[1][0] - O[0][1]) / (4 * q4);
        q[1] = (O[0][2] + O[2][0]) / (4 * q4);
        q[2] = (O[2][1] + O[1][2]) / (4 * q4);
        q[3] = q4;
    }
}

/* ----------------------------------------------------------------------
   define the symmetry operator for cubic or hexagonal
------------------------------------------------------------------------- */
void AppPottsRS::symmat(double ***sym) {
    // grow by number of symmetric operators
    (*sym) = new double *[Osym];

    // grow for symmetry quaternion vectors
    for (int o = 0; o < Osym; o++)
        (*sym)[o] = new double[4];

    // buffer for quaternion
    double q[4];

    if (Osym == 24) {
        // cubic symmetry
        double SYM[24][3][3] =
            {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
             {{1, 0, 0}, {0, -1, 0}, {0, 0, -1}},
             {{1, 0, 0}, {0, 0, -1}, {0, 1, 0}},
             {{1, 0, 0}, {0, 0, 1}, {0, -1, 0}},
             {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}},
             {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
             {{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}},
             {{-1, 0, 0}, {0, 0, 1}, {0, 1, 0}},
             {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
             {{0, 1, 0}, {0, 0, -1}, {-1, 0, 0}},
             {{0, 1, 0}, {1, 0, 0}, {0, 0, -1}},
             {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
             {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
             {{0, -1, 0}, {0, 0, -1}, {1, 0, 0}},
             {{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}},
             {{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}},
             {{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}},
             {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
             {{0, 0, 1}, {0, -1, 0}, {1, 0, 0}},
             {{0, 0, 1}, {-1, 0, 0}, {0, -1, 0}},
             {{0, 0, -1}, {0, 1, 0}, {1, 0, 0}},
             {{0, 0, -1}, {-1, 0, 0}, {0, 1, 0}},
             {{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}},
             {{0, 0, -1}, {1, 0, 0}, {0, -1, 0}}};

        // initialize global operator
        for (int o = 0; o < Osym; o++) {
            mat2quat(SYM[o], q);
            for (int i = 0; i < 4; i++)
                (*sym)[o][i] = q[i];
        }
    } else if (Osym == 12) {
        // hexagonal symmetry
        double a = sqrt(3) / 2;
        double SYM[12][3][3] =
            {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
             {{-0.5, a, 0}, {-a, -0.5, 0}, {0, 0, 1}},
             {{-0.5, -a, 0}, {a, -0.5, 0}, {0, 0, 1}},
             {{0.5, a, 0}, {-a, 0.5, 0}, {0, 0, 1}},
             {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
             {{0.5, -a, 0}, {a, 0.5, 0}, {0, 0, 1}},
             {{-0.5, -a, 0}, {-a, 0.5, 0}, {0, 0, -1}},
             {{1, 0, 0}, {0, -1, 0}, {0, 0, -1}},
             {{-0.5, a, 0}, {a, 0.5, 0}, {0, 0, -1}},
             {{0.5, a, 0}, {a, -0.5, 0}, {0, 0, -1}},
             {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}},
             {{0.5, -a, 0}, {-a, -0.5, 0}, {0, 0, -1}}};

        // initialize global operator
        for (int o = 0; o < Osym; o++) {
            mat2quat(SYM[o], q);
            for (int i = 0; i < 4; i++)
                (*sym)[o][i] = q[i];
        }
    }
}

double AppPottsRS::quaternions(const double qi[4], const double qj[4]) {
    double miso0, misom = MY_2PI;

    double q[4], qib[4], qjb[4], qmin[4] = {0, 0, 0, 0};
    for (int o1 = 0; o1 < Osym; o1++) {
        for (int o2 = 0; o2 < Osym; o2++) {
            quat_mult(symquat[o1], qi, qib);
            quat_mult(symquat[o2], qj, qjb);

            // j-grain conjugate quaternion
            qjb[1] = -qjb[1];
            qjb[2] = -qjb[2];
            qjb[3] = -qjb[3];

            quat_mult(qib, qjb, q);
            miso0 = 2 * acos(q[0]);

            if (miso0 > MY_PI)
                miso0 = miso0 - MY_2PI;
            if (fabs(miso0) < misom) {
                misom = fabs(miso0);
                qmin[0] = q[0];
                qmin[1] = q[1];
                qmin[2] = q[2];
                qmin[3] = q[3];
            }
        }
    }

    miso0 = 2 * acos(qmin[0]);
    if (miso0 > MY_PI)
        miso0 = miso0 - MY_2PI;

    return fabs(miso0);
}

void AppPottsRS::quat_mult(const double qi[4], const double qj[4], double q[4]) {
    // Hamilton multiplication/product
    // multiplying quaternions and update
    q[0] = qi[0] * qj[0] - qi[1] * qj[1] - qi[2] * qj[2] - qi[3] * qj[3];
    q[1] = qi[0] * qj[1] + qi[1] * qj[0] + qi[2] * qj[3] - qi[3] * qj[2];
    q[2] = qi[0] * qj[2] - qi[1] * qj[3] + qi[2] * qj[0] + qi[3] * qj[1];
    q[3] = qi[0] * qj[3] + qi[1] * qj[2] - qi[2] * qj[1] + qi[3] * qj[0];
}

void AppPottsRS::euler2quat(int i, double q[4]) {
    // Convert grain Euler angles to quaternion vector
    double p1 = phi1[i], P = Phi[i], p2 = phi2[i];

    q[0] = cos(P / 2.) * cos((p1 + p2) / 2.);
    q[1] = sin(P / 2.) * cos((p1 - p2) / 2.);
    q[2] = sin(P / 2.) * sin((p1 - p2) / 2.);
    q[3] = cos(P / 2.) * sin((p1 + p2) / 2.);
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */

void AppPottsRS::site_event_rejection(int i, RandomPark *random) {
    int oldstate = spin[i];
    double iphi[3] = {phi1[i], Phi[i], phi2[i]};

    // events = spin flips to neighboring site different than self

    int j, nei;
    int nevent = 0;

    // Nearest-neighbor sampling
    for (j = 0; j < numneigh[i]; j++) {
        nei = neighbor[i][j];
        if (spin[i] == spin[nei])
            continue;

        sites[nevent++] = nei;
    }

    if (nevent == 0) return;

    int iran = (int)(nevent * random->uniform());
    if (iran >= nevent) iran = nevent - 1;

    double einitial = site_energy(i), qold[4];
    euler2quat(i, qold);

    spin[i] = spin[sites[iran]];
    phi1[i] = phi1[sites[iran]];
    phi2[i] = phi2[sites[iran]];
    Phi[i] = Phi[sites[iran]];

    double efinal = site_energy(i), qnew[4];

    // Determing misorientation between ij states to calculate mobility
    double thetar;
    int smin = MIN(oldstate, spin[i]);
    int smax = MAX(oldstate, spin[i]);

    std::pair<int, int> spins = std::make_pair(smin, smax);

    // ratio of theta/theta_m
    if (misos.count(spins) == 1)
        thetar = misos[spins];
    else {
        euler2quat(i, qnew);
        thetar = quaternions(qold, qnew) / thetam;
        misos[spins] = thetar;
    }

    double p0 = (1.0 - exp(-bmob * pow(thetar, nmob)));

    // Check for isotropic case
    if (thetam < 1e-8) p0 = 1;

    // accept or reject via Boltzmann criterion
    if (efinal <= einitial) {
        if (thetar < 1e-8) {
            // Check if two instances of random Euler angle assignment were initiated
            // highly unlikely?
            if (fabs(phi1[i] - iphi[0]) < 1e-8 && fabs(Phi[i] - iphi[1]) < 1e-8 &&
                fabs(phi2[i] - iphi[2]) < 1e-8) {
                spin[i] = oldstate = MIN(spin[i], oldstate);
            }
        } else if (random->uniform() < p0) {
        } else {
            spin[i] = oldstate;
            phi1[i] = iphi[0];
            phi2[i] = iphi[2];
            Phi[i] = iphi[1];
        }
    } else if (temperature == 0.0) {
        spin[i] = oldstate;
        phi1[i] = iphi[0];
        phi2[i] = iphi[2];
        Phi[i] = iphi[1];
    } else if (random->uniform() > p0 * exp((einitial - efinal) * t_inverse)) {
        spin[i] = oldstate;
        phi1[i] = iphi[0];
        phi2[i] = iphi[2];
        Phi[i] = iphi[1];
    }

    if (spin[i] != oldstate) naccept++;

    // TODO: masking?
}
