#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=chipexoqual-%A.out
#SBATCH --error=chipexoqual-%A.out

args=()
args+=("$@")
if [ ! -z "$SLURM_CPUS_PER_TASK" ]
then
  args+=("--threads" "$SLURM_CPUS_PER_TASK")
fi

Rscript "$CHIPEXOQUAL_BASE"/chipexoqual.R "${args[@]}"

