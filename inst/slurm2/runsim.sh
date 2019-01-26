#!/bin/bash

#SBATCH --nodes=1

source ~/loadR.sh
Rscript sim.R
