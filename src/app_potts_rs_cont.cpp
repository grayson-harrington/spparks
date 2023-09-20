/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
-------------------------------------------------------------------------
    Continuous orientation state space Read-Shockley SPPARKS app_potts
        app_potts_rs_cont uses euler angles and quaternions to measure
        misorientation angles which are then used to calculate
        misorientation energy between neighboring voxels and to determine
        state change acceptance probability. This method can only be used
        via rKMC and MMC. Traditional KMC is not available.
    Developed by Grayson Harrington, PhD Student
    MINED group, Dr. Kalidindi
    Georgia Institute of Technology
------------------------------------------------------------------------- */

#include "app_potts_rs_cont.h"

#include <iostream>

#include "error.h"
#include "math.h"
#include "random_park.h"
#include "solve.h"
#include "stdlib.h"
#include "string.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsRSCont::AppPottsRSCont(SPPARKS *spk, int narg, char **arg) : AppPotts(spk, narg, arg) {
    ndouble = 7;
    allow_rejection = 1;
    allow_kmc = 0;
    allow_masking = 0;

    recreate_arrays();

    // error checking for AppPottsRSCont class
    if (strcmp(style, "potts_rs_cont") != 0)
        return;

    if (narg > 1)
        error->all(FLERR, "Illegal app_style command");
}

/* ---------------------------------------------------------------------- */

AppPottsRSCont::~AppPottsRSCont() {
    // make sure to delete variables with allocated memory
    //      none within this subclass
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */
void AppPottsRSCont::grow_app() {
    // set pointers to the arrays which contain euler angle and quaternion data
    phi1 = darray[0];
    phi = darray[1];
    phi2 = darray[2];
    quat_a = darray[3];
    quat_b = darray[4];
    quat_c = darray[5];
    quat_d = darray[6];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */
void AppPottsRSCont::init_app() {
    dt_sweep = 1.0;  // 1.0  because I want one sweep to be one MCS

    // make sure that all sites have been initialized
    double euler[3];
    double quat[4];
    for (int i = 0; i < nlocal; i++) {
        // generate euler angles and assign to sites
        generate_euler(euler);
        phi1[i] = euler[0];
        phi[i] = euler[1];
        phi2[i] = euler[2];

        // prepare quaternion and assign to sites
        euler2quaternion(euler, quat);
        quat_a[i] = quat[0];
        quat_b[i] = quat[1];
        quat_c[i] = quat[2];
        quat_d[i] = quat[3];
    }

    int flag = 0;
    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
    if (flagall)
        error->all(FLERR, "One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
    user defined optional parameters
 ------------------------------------------------------------------------- */
void AppPottsRSCont::input_app(char *command, int narg, char **arg) {
    // high-angle cutoff for Read-Shockley mmisorientation calculations
    if (strcmp(command, "theta_m") == 0) {
        if (narg != 2) error->all(FLERR, "Illegal high-angle cutoff angle command\n");

        double input = atof(arg[0]);
        char *units = arg[1];

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
            fprintf(logfile, "  theta_m set to %f radians\n", theta_m);
        if (screen && me == 0)
            fprintf(screen, "  theta_m set to %f radians\n", theta_m);
    }
    // boltzmann equation denominator, kbT (not physical boltzmann constant and temperature)
    else if (strcmp(command, "kbT") == 0) {
        if (narg != 1) error->all(FLERR, "Illegal kbT command\n");

        kbT = atof(arg[0]);
        if (kbT < 0) {
            error->all(FLERR, "kbT input must be positive\n");
        }
        kbT_inverse = 1 / kbT;

        if (logfile)
            fprintf(logfile, "  kbT set to %f\n", kbT);
        if (screen && me == 0)
            fprintf(screen, "  kbT set to %f\n", kbT);
    }
    // scalar component in Humphrey's general mobility formulation
    else if (strcmp(command, "humph_B") == 0) {
        if (narg != 1) error->all(FLERR, "Illegal humph_B command\n");

        humph_B = atof(arg[0]);

        if (logfile)
            fprintf(logfile, "  humph_B set to %f\n", humph_B);
        if (screen && me == 0)
            fprintf(screen, "  humph_B set to %f\n", humph_B);
    }
    // exponent component in Humphrey's general mobility formulation
    else if (strcmp(command, "humph_n") == 0) {
        if (narg != 1) error->all(FLERR, "Illegal humph_n command\n");

        humph_n = atof(arg[0]);

        if (logfile)
            fprintf(logfile, "  humph_n set to %f\n", humph_n);
        if (screen && me == 0)
            fprintf(screen, "  humph_n set to %f\n", humph_n);
    }
}

/* ----------------------------------------------------------------------
   update state variables with provided euler angles
        calculates quaternion and updates
------------------------------------------------------------------------- */
void AppPottsRSCont::update_state(int i, double euler_new[3], double quat_new[4]) {
    phi1[i] = euler_new[0];
    phi[i] = euler_new[1];
    phi2[i] = euler_new[2];

    quat_a[i] = quat_new[0];
    quat_b[i] = quat_new[1];
    quat_c[i] = quat_new[2];
    quat_d[i] = quat_new[3];
}

/* ----------------------------------------------------------------------
   compute energy of site
        energy of each neighbor contribute as the probability of two
        neighbors having the same orientation is negligible
------------------------------------------------------------------------- */
double AppPottsRSCont::site_energy(int i) {
    double quat_i[4] = {quat_a[i], quat_b[i], quat_c[i], quat_d[i]};

    double energy = 0.0;
    int nei;
    for (int j = 0; j < numneigh[i]; j++) {
        nei = neighbor[i][j];
        double quat_nei[4] = {quat_a[nei], quat_b[nei], quat_c[nei], quat_d[nei]};
        double miso = angle_between(quat_i, quat_nei);
        energy += misorientation_energy(miso);
    }
    return (double)energy;
}

/* ----------------------------------------------------------------------
   rKMC method
        perform a site event with probabilistic rejection
        flip to random orientation and test for energy change
------------------------------------------------------------------------- */
void AppPottsRSCont::site_event_rejection(int i, RandomPark *random) {
    // Save old state
    double euler_old[3] = {phi1[i], phi[i], phi2[i]};
    double quat_old[4] = {quat_a[i], quat_b[i], quat_c[i], quat_d[i]};

    // get initial energy
    double energy_old = site_energy(i);

    // generate new state
    double euler_new[3];
    double quat_new[4];
    generate_euler(euler_new);
    euler2quaternion(euler_new, quat_new);

    // change state, get change in energy
    update_state(i, euler_new, quat_new);
    double energy_new = site_energy(i);

    // accept or reject via Boltzmann criterion
    //      get mobility prefactor
    //          Humphrey's mobility: M ~ 1-exp(B(theta/theta_m)^n)
    //          theta, in this case, is the difference between old and new theta values
    double theta = angle_between(quat_old, quat_new);
    double prob_0 = 1 - exp(-humph_B * pow(theta / theta_m, humph_n));

    //      get probability of state change occuring
    double prob = 0.0;
    if (energy_new <= energy_old) {
        prob = 1.0;
    } else if (kbT == 0.0) {
        prob = 0.0;  // boltzmann case
    } else {
        prob = exp((energy_old - energy_new) * kbT_inverse);  // boltzmann general
    }
    prob = prob_0*prob;

    if (random->uniform() < prob) {
        // state change accepted
        naccept++;
    } else {
        // state change rejected
        update_state(i, euler_old, quat_old);
    }
}

/* ----------------------------------------------------------------------
   FUNCTIONS FOR ORIENTATION ANALYSIS
        generate_euler --> generate euler angle within fundamental zone
        euler2quaternion --> convert euler to quaternion
        angle_between --> angle between two quaternions
        misorientation_energy --> Read-Shockley misorientation energy
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   generate an euler angle within the allowable range of angles
        RADIANS
------------------------------------------------------------------------- */
void AppPottsRSCont::generate_euler(double euler_out[3]) {
    double phi1 = ((double)rand() / (RAND_MAX)) * 2 * M_PI;
    double phi = ((double)rand() / (RAND_MAX)) * M_PI;
    double phi2 = ((double)rand() / (RAND_MAX)) * 2 * M_PI;

    double euler_in[3] = {phi1, phi, phi2};
    cubic_euler2fz(euler_in, euler_out);
}

/* ----------------------------------------------------------------------
   calculate the angle between two quaternions
        https://iopscience.iop.org/article/10.1088/0965-0393/23/8/083501
            A.5 and A.6
------------------------------------------------------------------------- */
void AppPottsRSCont::euler2quaternion(const double euler_in[3], double quat_out[4]) {
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

/* ----------------------------------------------------------------------
   calculate the angle between two quaternions
        https://math.stackexchange.com/questions/90081/quaternion-distance
        https://link.springer.com/article/10.1007/s10851-009-0161-2
------------------------------------------------------------------------- */
double AppPottsRSCont::angle_between(const double quat1[4], const double quat2[4]) {
    double dot = quat1[0] * quat2[0] + quat1[1] * quat2[1] + quat1[2] * quat2[2] + quat1[3] * quat2[3];
    double angle = acos(2 * pow(dot, 2) - 1);
    if (isnan(angle)) {
        angle = 0;
    }
    return angle;
}

/* ----------------------------------------------------------------------
   calculate Read-Shockley misorientation energy
        https://journals.aps.org/pr/abstract/10.1103/PhysRev.78.275
------------------------------------------------------------------------- */
double AppPottsRSCont::misorientation_energy(const double theta) {
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

/* ----------------------------------------------------------------------
   given an input euler angle, find an equivalent fundamental zone euler
   angle based on 24-fold cubic symmetry
------------------------------------------------------------------------- */
void AppPottsRSCont::cubic_euler2fz(const double euler_in[3], double euler_out[3]) {
    double phi1 = euler_in[0];
    double phi = euler_in[1];
    double phi2 = euler_in[2];

    double temp_expr[48][3];
    double r = ((double)rand() / (RAND_MAX));

    phi = phi + r * 1e-7;  // why is this here???

    double t1 = sin(phi2);
    double t2 = cos(phi1);
    double t4 = cos(phi2);
    double t5 = cos(phi);
    double t6 = t4 * t5;
    double t7 = sin(phi1);
    double t10 = sin(phi);
    double t11 = pow(t10, 2);
    double t12 = pow(t4, 2);
    double t15 = sqrt(1 - t11 * t12);
    double t16 = 1 / t15;
    double t17 = (t1 * t2 + t6 * t7) * t16;
    double t21 = (t1 * t7 - t6 * t2) * t16;
    double t22 = atan2(t17, -t21);
    double t23 = t22 - phi1;
    double t24 = t4 * t10;
    double t25 = atan2(t15, -t24);
    double t26 = t25 - phi;
    double t27 = t1 * t10;
    double t28 = t27 * t16;
    double t29 = t5 * t16;
    double t30 = atan2(t28, t29);
    double t31 = t30 - phi2;
    double t32 = atan2(-t17, t21);
    double t33 = t32 - phi1;
    double t34 = atan2(-t15, -t24);
    double t35 = t34 - phi;
    double t36 = atan2(-t28, -t29);
    double t37 = t36 - phi2;
    double t38 = 2 * phi;
    double t39 = M_PI - t38;
    double t40 = 2 * phi2;
    double t41 = M_PI - t40;
    double t43 = t1 * t5;
    double t46 = pow(t1, 2);
    double t49 = sqrt(-t11 * t46 + 1);
    double t50 = 1 / t49;
    double t51 = (-t4 * t2 + t43 * t7) * t50;
    double t55 = (t4 * t7 + t43 * t2) * t50;
    double t56 = atan2(-t51, -t55);
    double t57 = t56 - phi1;
    double t58 = atan2(t49, t27);
    double t59 = t58 - phi;
    double t60 = t24 * t50;
    double t61 = t5 * t50;
    double t62 = atan2(t60, t61);
    double t63 = t62 - phi2;
    double t64 = atan2(t51, t55);
    double t65 = t64 - phi1;
    double t66 = atan2(-t49, t27);
    double t67 = t66 - phi;
    double t68 = atan2(-t60, -t61);
    double t69 = t68 - phi2;
    double t70 = atan2(t29, -t28);
    double t71 = t70 - phi2;
    double t72 = atan2(-t29, t28);
    double t73 = t72 - phi2;
    double t74 = atan2(t15, t24);
    double t75 = t74 - phi;
    double t76 = atan2(-t28, t29);
    double t77 = t76 - phi2;
    double t78 = atan2(-t15, t24);
    double t79 = t78 - phi;
    double t80 = atan2(t28, -t29);
    double t81 = t80 - phi2;
    double t82 = M_PI / 2;
    double t83 = atan2(t49, -t27);
    double t84 = t83 - phi;
    double t85 = atan2(-t60, t61);
    double t86 = t85 - phi2;
    double t87 = atan2(-t49, -t27);
    double t88 = t87 - phi;
    double t89 = atan2(t60, -t61);
    double t90 = t89 - phi2;
    double t91 = -t82 - t40;
    double t92 = t82 - t40;
    double t93 = atan2(-t61, t60);
    double t94 = t93 - phi2;
    double t95 = atan2(t61, -t60);
    double t96 = t95 - phi2;
    double t97 = atan2(-t61, -t60);
    double t98 = t97 - phi2;
    double t99 = atan2(t61, t60);
    double t100 = t99 - phi2;
    double t101 = atan2(-t29, -t28);
    double t102 = t101 - phi2;
    double t103 = atan2(t29, t28);
    double t104 = t103 - phi2;

    temp_expr[0][0] = t23;
    temp_expr[0][1] = t26;
    temp_expr[0][2] = t31;
    temp_expr[1][0] = t33;
    temp_expr[1][1] = t35;
    temp_expr[1][2] = t37;
    temp_expr[2][0] = -M_PI;
    temp_expr[2][1] = t39;
    temp_expr[2][2] = t41;
    temp_expr[3][0] = 0;
    temp_expr[3][1] = -M_PI;
    temp_expr[3][2] = -t40;
    temp_expr[4][0] = t57;
    temp_expr[4][1] = t59;
    temp_expr[4][2] = t63;
    temp_expr[5][0] = t65;
    temp_expr[5][1] = t67;
    temp_expr[5][2] = t69;
    temp_expr[6][0] = t23;
    temp_expr[6][1] = t26;
    temp_expr[6][2] = t71;
    temp_expr[7][0] = t33;
    temp_expr[7][1] = t35;
    temp_expr[7][2] = t73;
    temp_expr[8][0] = t33;
    temp_expr[8][1] = t75;
    temp_expr[8][2] = t77;
    temp_expr[9][0] = t23;
    temp_expr[9][1] = t79;
    temp_expr[9][2] = t81;
    temp_expr[10][0] = 0;
    temp_expr[10][1] = 0;
    temp_expr[10][2] = -M_PI;
    temp_expr[11][0] = -M_PI;
    temp_expr[11][1] = -t38;
    temp_expr[11][2] = 0;
    temp_expr[12][0] = t23;
    temp_expr[12][1] = t26;
    temp_expr[12][2] = t37;
    temp_expr[13][0] = t33;
    temp_expr[13][1] = t35;
    temp_expr[13][2] = t31;
    temp_expr[14][0] = -M_PI;
    temp_expr[14][1] = t39;
    temp_expr[14][2] = -t40;
    temp_expr[15][0] = 0;
    temp_expr[15][1] = -M_PI;
    temp_expr[15][2] = t41;
    temp_expr[16][0] = t57;
    temp_expr[16][1] = t59;
    temp_expr[16][2] = t69;
    temp_expr[17][0] = t65;
    temp_expr[17][1] = t67;
    temp_expr[17][2] = t63;
    temp_expr[18][0] = 0;
    temp_expr[18][1] = 0;
    temp_expr[18][2] = -t82;
    temp_expr[19][0] = -M_PI;
    temp_expr[19][1] = -t38;
    temp_expr[19][2] = t82;
    temp_expr[20][0] = t65;
    temp_expr[20][1] = t84;
    temp_expr[20][2] = t86;
    temp_expr[21][0] = t57;
    temp_expr[21][1] = t88;
    temp_expr[21][2] = t90;
    temp_expr[22][0] = t23;
    temp_expr[22][1] = t26;
    temp_expr[22][2] = t73;
    temp_expr[23][0] = t33;
    temp_expr[23][1] = t35;
    temp_expr[23][2] = t71;
    temp_expr[24][0] = -M_PI;
    temp_expr[24][1] = t39;
    temp_expr[24][2] = t91;
    temp_expr[25][0] = 0;
    temp_expr[25][1] = -M_PI;
    temp_expr[25][2] = t92;
    temp_expr[26][0] = t57;
    temp_expr[26][1] = t59;
    temp_expr[26][2] = t94;
    temp_expr[27][0] = t65;
    temp_expr[27][1] = t67;
    temp_expr[27][2] = t96;
    temp_expr[28][0] = t65;
    temp_expr[28][1] = t84;
    temp_expr[28][2] = t98;
    temp_expr[29][0] = t57;
    temp_expr[29][1] = t88;
    temp_expr[29][2] = t100;
    temp_expr[30][0] = t33;
    temp_expr[30][1] = t75;
    temp_expr[30][2] = t81;
    temp_expr[31][0] = t23;
    temp_expr[31][1] = t79;
    temp_expr[31][2] = t77;
    temp_expr[32][0] = 0;
    temp_expr[32][1] = 0;
    temp_expr[32][2] = 0;
    temp_expr[33][0] = -M_PI;
    temp_expr[33][1] = -t38;
    temp_expr[33][2] = -M_PI;
    temp_expr[34][0] = t33;
    temp_expr[34][1] = t75;
    temp_expr[34][2] = t102;
    temp_expr[35][0] = t23;
    temp_expr[35][1] = t79;
    temp_expr[35][2] = t104;
    temp_expr[36][0] = t57;
    temp_expr[36][1] = t59;
    temp_expr[36][2] = t96;
    temp_expr[37][0] = t65;
    temp_expr[37][1] = t67;
    temp_expr[37][2] = t94;
    temp_expr[38][0] = t33;
    temp_expr[38][1] = t75;
    temp_expr[38][2] = t104;
    temp_expr[39][0] = t23;
    temp_expr[39][1] = t79;
    temp_expr[39][2] = t102;
    temp_expr[40][0] = t65;
    temp_expr[40][1] = t84;
    temp_expr[40][2] = t100;
    temp_expr[41][0] = t57;
    temp_expr[41][1] = t88;
    temp_expr[41][2] = t98;
    temp_expr[42][0] = 0;
    temp_expr[42][1] = 0;
    temp_expr[42][2] = t82;
    temp_expr[43][0] = -M_PI;
    temp_expr[43][1] = -t38;
    temp_expr[43][2] = -t82;
    temp_expr[44][0] = t65;
    temp_expr[44][1] = t84;
    temp_expr[44][2] = t90;
    temp_expr[45][0] = t57;
    temp_expr[45][1] = t88;
    temp_expr[45][2] = t86;
    temp_expr[46][0] = -M_PI;
    temp_expr[46][1] = t39;
    temp_expr[46][2] = t92;
    temp_expr[47][0] = 0;
    temp_expr[47][1] = -M_PI;
    temp_expr[47][2] = t91;

    double euler_equivs[48][3];
    for (int i = 0; i < 48; i++) {
        euler_equivs[i][0] = fmod(temp_expr[i][0] + phi1, 2 * M_PI);
        euler_equivs[i][0] = euler_equivs[i][0] >= 0 ? euler_equivs[i][0] : euler_equivs[i][0] + 2 * M_PI;
        euler_equivs[i][1] = fmod(temp_expr[i][1] + phi, 2 * M_PI);
        euler_equivs[i][1] = euler_equivs[i][1] >= 0 ? euler_equivs[i][1] : euler_equivs[i][1] + 2 * M_PI;
        euler_equivs[i][2] = fmod(temp_expr[i][2] + phi2, 2 * M_PI);
        euler_equivs[i][2] = euler_equivs[i][2] >= 0 ? euler_equivs[i][2] : euler_equivs[i][2] + 2 * M_PI;

        // check to see if this euler equiv is within the fundamental zone.
        // if it is, then return it.
        double phi2_cutoff_high = M_PI / 4;
        double Phi_cutoff_high = M_PI / 2;
        double Phi_cutoff_low = acos(cos(euler_equivs[i][2]) / sqrt(1 + pow(cos(euler_equivs[i][2]), 2)));
        if ((euler_equivs[i][2] <= phi2_cutoff_high) && (euler_equivs[i][1] <= Phi_cutoff_high) && (euler_equivs[i][1] > Phi_cutoff_low)) {
            euler_out[0] = euler_equivs[i][0];
            euler_out[1] = euler_equivs[i][1];
            euler_out[2] = euler_equivs[i][2];
            break;
        }
    }
}
