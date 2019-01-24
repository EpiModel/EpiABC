#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time=03:00:00
#SBATCH --mem=55G

. /suppscr/csde/sjenness/spack/share/spack/setup-env.sh
module load gcc_4.4.7-ompi_1.8.8
module load r-3.5.1-gcc-8.2.0-jjtszuq

mpirun -np 1 Rscript sim.R
