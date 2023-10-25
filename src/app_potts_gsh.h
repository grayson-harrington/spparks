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

#include <string>
#include <vector>

#include "app_potts.h"

namespace SPPARKS_NS {

class NeuralNetwork {
   public:
    NeuralNetwork();
    NeuralNetwork(char *filePath);
    ~NeuralNetwork();

    // weights and biases for the hidden and output layers
    int inputSize;
    int hiddenSize;
    int outputSize;

    std::vector<std::vector<double>> weightsInputHidden;
    std::vector<std::vector<double>> weightsHiddenOutput;
    std::vector<double> biasHidden;
    std::vector<double> biasOutput;

    // activation functions
    std::string hidden_activation;
    std::string output_activation;

    double sigmoid(double x);
    double relu(double x);

    // Forward pass through the network
    std::vector<double> forward(const std::vector<double> &input);
};

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

    // rKMC functions
    virtual double site_energy(int);
    virtual void site_event_rejection(int, class RandomPark *);

   protected:
    int nspins;
    int *spin;

    bool init_micro = false;

    // spin maps
    int n_euler_angles, n_gsh_coef;
    double **spin2euler, **spin2gsh;
    // returns the maps (spin to euler, spin to gsh)
    SpinMaps read_spin2angle_map(const char *filePath, int &n_lines, int &n_eul, int &n_gsh);

    // nn_energy_function
    NeuralNetwork nn;

    // rbf energy function
    std::vector<std::vector<double>> rbf_positions;
    std::vector<double> rbf_scales;
    void read_rbf_input(const std::string& filename, std::vector<std::vector<double>>& positions, std::vector<double>& scales);
    double rbf_energy_function(std::vector<std::vector<double>> positions, std::vector<double> scales, std::vector<double> x);

    // HELPER FUNCTIONS

    // euclidean distance between two double*
    double euclideanDistance(const double *array1, const double *array2, const int size);
};

}  // namespace SPPARKS_NS

#endif
#endif
