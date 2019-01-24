#!/bin/bash

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:30:00
#SBATCH --mem=50G

source ~/loadR.sh
Rscript sim.R

mpirun -np 1 Rscript sim.R
