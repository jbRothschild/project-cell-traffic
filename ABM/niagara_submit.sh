#!/bin/bash
#SBATCH --nodes=40
#SBATCH --ntasks-per-node=40
#SBATCH --time=6:00:00
#SBATCH --job-name=ECvEC

#run this code using jbroths:~$ sbatch *script_name.sh*

# DIRECTORY TO RUN - $SLURM_SUBMIT_DIR is directory job was submitted from
cd $SLURM_SUBMIT_DIR

# load modules (must match modules used for compilation)
module load NiaEnv/2019b gnu-parallel

# Turn off implicit threading in Python, R
export OMP_NUM_THREADS=1

# EXECUTION COMMAND; ampersand off 40 jobs and wait

# mkdir -p ${RESULTS_DIR}/${SIM_DIR}

parallel --joblog slurm-$SLURM_JOBID.log --sshdelay 0.1 --wd $PWD "./des42.o 102 {}" ::: {0..3199}
