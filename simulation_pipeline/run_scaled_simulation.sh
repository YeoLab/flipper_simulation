#!/bin/bash
unset SLURM_JOB_ID

snakemake -ks snakemake_sim_scale.py \
    --configfile configs/test.yaml \
    --profile profiles
