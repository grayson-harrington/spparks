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
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(potts/rs,AppPottsRS)

#else

#ifndef SPK_APP_POTTS_RS_H
#define SPK_APP_POTTS_RS_H

#include <map>
#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsRS : public AppPotts {
 public:

  AppPottsRS(class SPPARKS *, int, char **);
  ~AppPottsRS();

  void init_app();
  void grow_app();
  void input_app(char *, int, char **);

  virtual void site_event_rejection(int, class RandomPark *);

  double site_energy(int);

 protected:
  double *phi1,*phi2,*Phi; //pointer to Q-states for the 3 rotation angles
  int *spin;
  double thetam;           //High-low angle divider
  double Jij;              //Interaction energy
  int Osym;                //Symmetry Operator flag
  double **symquat;        //Symmetry Operator in quaternion space

  //Mobility = Mm * [1 - exp(-B * {theta/theham}^n)]
  double nmob;            //Mobility exponential power, n
  double bmob;            //Mobility exponential scaling, B

  //Get misorientation angle from quaternions
  double quaternions(const double qi[4], const double qj[4]);

  //Multiplication between quaternion vectors
  void quat_mult(const double qi[4], const double qj[4], double q[4]);

  //Define the symmetry operator based on symmetry flag
  void symmat(double ***);

  //Convert symmetry operator into quaternion space
  void mat2quat(const double O[3][3], double q[4]);
  void euler2quat(int i, double q[4]);

  //map to store misorientations
  std::map<std::pair<int,int>, double> misos;
};
 

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

*/
