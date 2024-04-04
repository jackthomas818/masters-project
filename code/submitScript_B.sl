#!/bin/bash -l

#SBATCH -J simStudy10
#SBATCH --array=1-1728
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --mail-user=jackthom@uvic.ca
#SBATCH --mail-type=ALL
#SBATCH --output=slurm.out
#SBATCH --account=def-lcowen

echo "The current directory is `pwd`"
echo "Start at `date`"

module load r/4.1.2
export R_LIBS=/home/jackthom/R/

Rscript ./runsim_B.R $SLURM_ARRAY_TASK_ID 50 >& ./log_${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}.txt

echo "Done at `date`"
