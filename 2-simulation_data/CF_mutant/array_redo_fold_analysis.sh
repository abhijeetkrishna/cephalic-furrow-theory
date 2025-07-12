#!/bin/bash

#SBATCH -J wd
#SBATCH -o log_files/redo_out_%A_%a
#SBATCH -e log_files/redo_err_%A_%a
#SBATCH -t 10:00:00
#SBATCH -a 0-1620%700


#py_map_creator=map_index_cf.py
#initiator=initialize.py
#simulation=gradDescent.jl
map_file=map_index_47446914.csv
postprocess=postprocess_redo.py

#if [ ! -f $map_file ]
#then
#	python3 $py_map_creator $SLURM_ARRAY_JOB_ID
#fi

#python3 $initiator $map_file $SLURM_ARRAY_TASK_ID
#julia $simulation $map_file $SLURM_ARRAY_TASK_ID
python3 $postprocess $map_file $SLURM_ARRAY_TASK_ID