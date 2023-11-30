#!/bin/bash

#SBATCH --partition=medium                   ## Which Partition/Queue to use
##SBATCH --partition=short                   ## Which Partition/Queue to use
##SBATCH --account=short                    ## must match partition
#SBATCH --account=medium                    ## must match partition
#SBATCH --time=48:00:00         # adjust this to match the walltime of your job
##SBATCH --time=00:20:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36     # adjust this if you are using parallel commands
#SBATCH -w node004

matlab -nosplash -nodisplay -r run_get_shelf_collapse_time_HPC_variablem
