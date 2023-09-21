/* -----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(potts_gsh, AppPottsGSH)
#else

#ifndef SPK_APP_POTTS_GSH_H
#define SPK_APP_POTTS_GSH_H

#include "app_potts.h"

namespace SPPARKS_NS {

    struct SpinMaps {
        double **spin2euler;
        double **spin2gsh;
    };

    class AppPottsGSH : public AppPotts {
    public:
        AppPottsGSH(class SPPARKS *, int, char **);
        virtual ~AppPottsGSH();
        virtual void grow_app();
        virtual void init_app();
        virtual void input_app(char *, int, char **);

        virtual double site_energy(int);

        //  TODO: I might be able to remove these from the header and cpp files if they are just inherited from AppPotts
        // wait to see if you will be modifying any of them
        virtual void site_event_rejection(int, class RandomPark *);
        virtual double site_propensity(int);

    protected:
        int nspins;
        int *spin;
        int *sites, *unique;

        // spin maps 
        int n_euler_angles, n_gsh_coef;
        double **spin2euler, **spin2quat, **spin2gsh;
    
        // site energy functions and variables
        //      type:
        //          1 - spin count
        //          2 - euler misorientation
        //          3 - gsh distance
        int site_energy_type; 

        double site_energy_spin(int i);

        double theta_m;
        double read_shockley(const double theta);
        double site_energy_read_shockley(int i);

        double site_energy_gsh(int i);
        
        // HELPER FUNCTIONS

        // returns the maps (spin to euler, spin to gsh)
        SpinMaps read_spin2angle_map(const char *filePath, int &n_lines, int &n_eul, int &n_gsh);

        // euclidean distance between two double*
        double euclideanDistance(const double* array1, const double* array2, const int size);

        // convert euler angles into quaternion
        void euler2quaternion(const double euler_in[3], double quat_out[4]);

        double angle_between(const double quat1[4], const double quat2[4]);
        
    };
}  // namespace SPPARKS_NS

#endif
#endif
