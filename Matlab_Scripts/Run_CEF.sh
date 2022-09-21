#!/bin/bash

#SBATCH -n 7 
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH -q normal 

module purge all

module load matlab/2021a

matlab -nosplash -nodesktop -noFigureWindows < baseline_script.m
matlab -nosplash -nodesktop -noFigureWindows < get_errors.m 
