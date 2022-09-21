#!/bin/bash

#SBATCH -n 56 
#SBATCH -N 2
#SBATCH -t 8:00:00
##SBATCH -p htcgpu2
#SBATCH -q sawilso6
##SBATCH -A sawilso6_sec

module purge all

module load matlab/2021a

matlab -nosplash -nodesktop -noFigureWindows < baseline_script.m               
