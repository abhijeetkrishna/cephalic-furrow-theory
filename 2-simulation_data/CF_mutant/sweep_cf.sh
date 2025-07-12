#!/bin/bash
#SBATCH -J sweep_cf     # the job's name
#SBATCH -t 10:00:00       # max. wall clock time 5s
#SBATCH -n 1              # number of tasks
#SBATCH -o log_files/sweep_out  # output file
#SBATCH -e log_files/sweep_err  # output file
#SBATCH --partition=batch

python3 postprocess_sweep.py