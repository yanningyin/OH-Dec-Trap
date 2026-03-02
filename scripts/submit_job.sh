#!/bin/bash
#SBATCH --job-name=matlab_job            # Job name
#SBATCH --output=logs/output_%j.out           # Standard output and error log
#SBATCH --error=logs/error_%j.err             # Error log
#SBATCH --nodes=1
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G                        # Memory per node
#SBATCH --partition=wcpu

# export COMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Send some noteworthy information to the output log
echo "Running on node: $(hostname)"
echo "In directory:    $(pwd)"
echo "Starting on:     $(date)"
echo "SLURM_JOB_ID:    ${SLURM_JOB_ID}"

# Load the MATLAB module
# module load matlab/matlab-R2021b
module load matlab

# Change to the directory where the MATLAB script is located
# cd ..

# Run MATLAB script without desktop and no display
matlab -nodisplay -nosplash -r "run('main.m'); exit;"

# Send more noteworthy information to the output log
echo "Finished at:     $(date)"

# End the script with exit code 0
exit 0
