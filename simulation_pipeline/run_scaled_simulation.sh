#!/bin/bash
set -euo pipefail
unset SLURM_JOB_ID

snakemake -ks snakemake_sim_scale.py --profile profiles -d ./results/
