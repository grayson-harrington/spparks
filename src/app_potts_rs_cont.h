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
    Developed by Grayson Harrington
        Georgia Institute of Technology, PhD Student
            MINED Research Group, Advisor: Dr. Surya Kalidindi
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(potts_rs_cont, AppPottsRSCont)
#else

#ifndef SPK_APP_POTTS_RS_CONT_H
#define SPK_APP_POTTS_RS_CONT_H

#include "app_potts.h"

#include <map>
#include <string>
#include <array>

namespace SPPARKS_NS
{

    class AppPottsRSCont : public AppPotts
    {
    public:
        AppPottsRSCont(class SPPARKS *, int, char **);
        ~AppPottsRSCont();

        // define app initialization methods
        void grow_app();
        void init_app();
        void input_app(char *, int, char **);

        // define required rKMC/MMC methods
        double site_energy(int);
        virtual void site_event_rejection(int, class RandomPark *);

    protected:
        // define all custom variables here...
        // euler angles, quats, etc...
        double *phi1, *phi, *phi2;
        double *quat_a, *quat_b, *quat_c, *quat_d;
        void update_state(int i, double euler_new[3], double quat_new[4]);

        // orientation variables and functions
        double theta_m;
        double humph_B, humph_n;
        double kbT, kbT_inverse;

        //      create euler angles which lie in the fundamental zone for cubic systems
        void generate_euler(double euler_out[3]);

        //      convert euler angles to the cubic fundamental zone
        void cubic_euler2fz(const double euler_in[3], double euler_out[3]);

        //      convert euler angles to quaternion
        void euler2quaternion(const double euler_in[3], double quat_out[4]);

        //      angle between two quaternions
        double angle_between(const double quat1[4], const double quat2[4]);

        //      Read-Shockley misorientation energy
        //          calculate the misorientation energy based on
        //          RS provided the misorientation angle
        double misorientation_energy(const double theta);
    };
}

#endif
#endif
