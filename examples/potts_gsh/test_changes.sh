#!/bin/bash

cd ../../src
make clean
make mpi
cd ../examples/potts_gsh/
mpirun -np 4 ../../src/spk_mpi < in.potts_gsh

