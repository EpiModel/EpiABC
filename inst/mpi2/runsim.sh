#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time=03:00:00
#SBATCH --mem=55G

source ~/loadR.sh
Rscript sim.R

mpirun -np 1 Rscript sim.R
