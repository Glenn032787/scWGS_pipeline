#!/bin/bash
#SBATCH -p upgrade
#SBATCH --job-name='run_snakemake'
#SBATCH --mail-type=ALL
#SBATCH --mem=1G
#SBATCH -n 1

eval "$(conda shell.bash hook)"
conda activate snakemake

if [[ $1 == "slurm" ]]; then
    slurm="--profile=slurm"
    mkdir -p slurm_logs
fi

snakemake ${slurm} -k -c 72 --use-singularity --singularity-args "-B /projects,/gsc,/home"  --rerun-incomplete
