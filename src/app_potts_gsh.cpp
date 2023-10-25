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
#include <sstream>

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

    if (narg != 3) error->all(FLERR, "Illegal app_style command");

    // get the map which is used to convert spin to euler angles and gsh coef
    SpinMaps result = read_spin2angle_map(arg[1], nspins, n_euler_angles, n_gsh_coef);
    spin2euler = result.spin2euler;
    spin2gsh = result.spin2gsh;


    // TODO:  decide nn or rbf

    // // load in neural network energy function
    // nn = NeuralNetwork(arg[2]);

    // load in the rbf energy function data
    read_rbf_input(arg[2], rbf_positions, rbf_scales);

    dt_sweep = 1.0;
    if (nspins <= 0) error->all(FLERR, "Illegal app_style command");
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
    if (narg < 1) {
        error->all(FLERR, "Invalid command for app_style");
    }

    if (strcmp(command, "read_sites") == 0) {
        init_micro = true;
    }
    // add parsing here
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

// gsh euclidean distance site energy
// double AppPottsGSH::site_energy(int i) {

//     double gsh_dist;
//     double eng = 0.0;

//     double* gsh_isite = spin2gsh[spin[i]];
//     double* gsh_jsite;

//     int nei;
//     for (int j = 0; j < numneigh[i]; j++) {
//         nei = neighbor[i][j];
//         gsh_jsite = spin2gsh[spin[nei]];

//         // Can choose between euclidean distance on just 3 values or on all 130

//         // eng += (gsh_isite[1]-gsh_jsite[1])*(gsh_isite[1]-gsh_jsite[1]) + (gsh_isite[5]-gsh_jsite[5])*(gsh_isite[5]-gsh_jsite[5]) + (gsh_isite[105]-gsh_jsite[105])*(gsh_isite[105]-gsh_jsite[105]);
//         eng += euclideanDistance(gsh_isite, gsh_jsite, n_gsh_coef);
//     }

//     return eng;
// }

// gsh energy from nn
// double AppPottsGSH::site_energy(int i) {

//     double nn_energy;
//     double eng = 0.0;
//     std::vector<double> nn_input;

//     double* gsh_isite = spin2gsh[spin[i]];
//     double* gsh_jsite;

//     int nei;
//     for (int j = 0; j < numneigh[i]; j++) {
//         nei = neighbor[i][j];
//         gsh_jsite = spin2gsh[spin[nei]];
        
//         // gsh misorientation depends on 1, 5, and 105 from the single crystral gsh values
//         nn_input = {gsh_isite[1], gsh_isite[5], gsh_isite[105],
//                       gsh_jsite[1], gsh_jsite[5], gsh_jsite[105]};
        
//         bool same = true;
//         for (int i = 0; i < 3; i++) {
//             if (nn_input[i] != nn_input[i+3]) {
//                 same = false;
//             }
//         }

//         if (same) {
//             eng += 0;
//         } else {
//             eng += nn.forward(nn_input)[0];
//         }
//     }

//     return eng;
// }

// gsh energy from rbfs
double AppPottsGSH::site_energy(int i) {

    double nn_energy;
    double eng = 0.0;
    std::vector<double> nn_input;

    double* gsh_isite = spin2gsh[spin[i]];
    double* gsh_jsite;

    int nei;
    for (int j = 0; j < numneigh[i]; j++) {
        nei = neighbor[i][j];
        gsh_jsite = spin2gsh[spin[nei]];

        // get rbf input
        std::vector<double> x;
        bool same = true;
        for (int i = 1; i <= 9; i++) {
            x.push_back(gsh_isite[i]);
        }
        for (int i = 1; i <= 9; i++) {
            if (x[i-1] != gsh_jsite[i]) {
                same = false;
            }
            x.push_back(gsh_jsite[i]);
        }

        if (same) {
            eng += 0;
        } else {
            eng += rbf_energy_function(rbf_positions, rbf_scales, x);
        }
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
    int j = (int)(numneigh[i] * random->uniform());
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

    for (int i = 1; i < size; ++i) {
        double diff = array1[i] - array2[i];
        sum += diff * diff;

        if (i == 9) break; // first 9 gsh values (ignoring the 1)
    }

    return sqrt(sum);
}

void AppPottsGSH::read_rbf_input(const std::string& filename, std::vector<std::vector<double>>& positions, std::vector<double>& scales) {
    positions.clear();
    scales.clear();
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }
    
    std::string line;
    int lineCount = 0;
    int n_positions = 0;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        double value;
        std::vector<double> data;

        while (ss >> value) {
            data.push_back(value);
        }

        if (lineCount == 0) {
            if (data.size() == 1) {
                n_positions = data[0];
            } else {
                std::cerr << "Error: Invalid input format." << std::endl;
                return;
            }
        } else {
            if (lineCount <= n_positions) {
                positions.push_back(data);
            } else {
                scales = data;
            }
        }
        
        lineCount++;
    }
}

double AppPottsGSH::rbf_energy_function(std::vector<std::vector<double>> positions, std::vector<double> scales, std::vector<double> x) {
    double result = 0.0;

    int n_rbfs = positions.size();
    
    for (size_t i = 0; i < positions.size(); i++) {
        double n = 0.0;
        for (size_t j = 0; j < positions[i].size(); j++) {
            n += scales[i] * pow(x[j] - positions[i][j], 2);
        }
        result -= exp(-pow(n, 2));
    }

    double eng = (result + n_rbfs)/n_rbfs;

    return eng;
}

/*
-------------------------------------------------------------------------
-------------------------------------------------------------------------
------------------ Feed Forward Neural Network --------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
*/

// function to create a neural network from file
// this shouldn't be used, but was necessary for compiling
NeuralNetwork::NeuralNetwork() {}

// function to create a neural network from file
NeuralNetwork::NeuralNetwork(char* filePath) {
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;

        iss >> token;

        if (token == "nn_shape") {
            iss >> inputSize >> hiddenSize >> outputSize;

            weightsInputHidden.resize(inputSize, std::vector<double>(hiddenSize, 0.0));
            weightsHiddenOutput.resize(hiddenSize, std::vector<double>(outputSize, 0.0));
            biasHidden.resize(hiddenSize, 0.0);
            biasOutput.resize(outputSize, 0.0);

        } else if (token == "hidden_activation") {
            iss >> hidden_activation;
            if (hidden_activation != "relu" && hidden_activation != "sigmoid" && hidden_activation != "none") {
                std::cerr << "Invalid activation type for input to hidden layer. Only relu, sigmoid, and none are available" << std::endl;
                exit(1);
            }
        } else if (token == "output_activation") {
            iss >> output_activation; // Use std::string directly
            if (output_activation != "relu" && output_activation != "sigmoid" && output_activation != "none") {
                std::cerr << "Invalid activation type for input to hidden layer. Only relu, sigmoid, and none are available" << std::endl;
                exit(1);
            }
        } else if (token == "input_hidden_weights") {
            for (int i = 0; i < inputSize; ++i) {
                std::getline(file, line);
                std::istringstream iss(line);
                for (int j = 0; j < hiddenSize; ++j) {
                    double weight;
                    iss >> weight;
                    weightsInputHidden[i][j] = weight;
                }
            }
        } else if (token == "hidden_bias") {
            std::getline(file, line);
            std::istringstream iss(line);
            for (int i = 0; i < hiddenSize; ++i) {
                double bias;
                iss >> bias;
                biasHidden[i] = bias;
            }
        } else if (token == "hidden_output_weights") {
            for (int i = 0; i < hiddenSize; ++i) {
                std::getline(file, line);
                std::istringstream iss(line);
                for (int j = 0; j < outputSize; ++j) {
                    double weight;
                    iss >> weight;
                    weightsHiddenOutput[i][j] = weight;
                }
            }
        } else if (token == "output_bias") {
            std::getline(file, line);
            std::istringstream iss(line);
            for (int i = 0; i < outputSize; ++i) {
                double bias;
                iss >> bias;
                biasOutput[i] = bias;
            }
        }
    }

    file.close();
}

NeuralNetwork::~NeuralNetwork() {
    // strings and vectors are memory managed already by their 
    // own classes, so no need to deallocate them here.
}

double NeuralNetwork::sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}

double NeuralNetwork::relu(double x) {
    return (x > 0.0) ? x : 0.0;
}

std::vector<double> NeuralNetwork::forward(const std::vector<double>& input) {
    // Forward pass from input to hidden layer
    std::vector<double> hidden(hiddenSize);
    for (int i = 0; i < hiddenSize; ++i) {
        hidden[i] = 0.0;
        for (int j = 0; j < inputSize; ++j) {
            hidden[i] += input[j] * weightsInputHidden[j][i];
        }
        hidden[i] += biasHidden[i];
        if (hidden_activation == "relu") {
            hidden[i] = relu(hidden[i]);
        } else if (hidden_activation == "sigmoid") {
            hidden[i] = sigmoid(hidden[i]);
        } else {
            // nothing to do here. no activation function applied
        }
    }

    // Forward pass from hidden to output layer
    std::vector<double> output(outputSize);
    for (int i = 0; i < outputSize; ++i) {
        output[i] = 0.0;
        for (int j = 0; j < hiddenSize; ++j) {
            output[i] += hidden[j] * weightsHiddenOutput[j][i];
        }
        output[i] += biasOutput[i];
        if (output_activation == "relu") {
            output[i] = relu(output[i]);
        } else if (output_activation == "sigmoid") {
            output[i] = sigmoid(output[i]);
        } else {
            // nothing to do here. no activation function applied
        }
    }

    return output;
}
