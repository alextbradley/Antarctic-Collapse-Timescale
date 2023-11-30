#!/bin/bash 

#SBATCH --partition=medium                   ## Which Partition/Queue to use
#SBATCH --account=medium                    ## must match partition
#SBATCH --time=48:00:00         # adjust this to match the walltime of your job
##SBATCH --time=00:20:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32     # adjust this if you are using parallel commands

module load hpc/matlab
FNAME="runsim_shelf_hpc($SHELFNO,$STEP)"
echo $FNAME
matlab -nosplash -nodisplay -r $FNAME
