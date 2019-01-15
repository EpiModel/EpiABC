#!/bin/bash

cd inst/slurm2
scp *.* hyak:/gscratch/csde/sjenness/epiabc

sbatch -p csde -A csde --array=1-32 --export=ALL,wave=0 runsim.sh
sbatch -p ckpt -A csde-ckpt --array=1-32 --export=ALL,wave=0 runsim.sh
