# Simulations for Flipper publication

This repository holds all code used to generate the simulation data for [https://doi.org/10.64898/2026.03.13.711628](https://doi.org/10.64898/2026.03.13.711628)

This includes both the jupyter notebooks used for running individual simulations and comparing simulated and real data distributions (supplementary figure 2), and a snakemake pipeline used for running many simultaneous simulations and quantifying the performance of different differential methods (figure 5). 

# Simulation_notebooks

## run_simulation
This repository includes a notebook for running the simulation as well as 4 additional notebooks for running differential eCLIP analyses and calculating performance metrics.

### 1_sim_notebook.ipynb
This notebook runs the core simulation. For details of the simulation methodology, see the methods section of [https://doi.org/10.64898/2026.03.13.711628](https://doi.org/10.64898/2026.03.13.711628). 

Simulating eCLIP is complex, and requires many seperate parameters. The tables below break these parameters up into 4 main sections. The parameters used to generate the results from the manuscript can be found in the notebook. 

#### Simulation setup parameters
| Parameter    | Description |
| -------- | ------- |
| seed | A random seed used to ensure reproducibility of the simulation. Can be in order to observe variation in simulation runs |
| conditions | The names to use for each condition. Should generally be left as the default ("control", "treatment") |
| output_dir | Where to save the simulation outputs |
| reps  |  The number of reps to simulate (minimum of 2) |
| n_genes | Total number of genes to simulate |
| n_bind | Total number of genes that contain a binding site |

#### Expression simulation parameters. 
| Parameter    | Description |
| -------- | ------- |
| beta_shift | Controls the direction of the treatment effect on expression. For example, a beta shift value of 1 corresponds to a net increase in expression following treatment |
| beta_sd | Controls the deviation of the treatment's effect on expression |
| Mean_ex | Controls the mean gene expression across all genes. |
| SD_ex | Controls the standard deviation of the gene expression across all genes.  |
| dispMeanRel | Controls the dispersion of the mean counts between different samples. |
| ex_size_factors | A list of size factors to be applied to the expression data |

**IMPORTANT NOTE:** Samples are always simulated in a specific pattern based on the number of reps used. For example, if you are running a simulation with reps = 2 and conditions = ("ctrl", "trt"), then you will be simulating 8 total samples with the following format:
| ctrl IN 1  | ctrl IP 1 | ctrl IN 2 | ctrl IP 2 | trt IN 1 | trt IP 1 | trt IN 2 | trt IP 2 |
| -| -| -| - | - | - | - | - |
| | | | | | | | |

This is important to remember when specifying size factors, as the size factors used will be applied based on the same order. 

#### eCLIP background simulation parameters. 
| Parameter    | Description |
| -------- | ------- |
| exon_frac | Fraction of the genome to be treated as exonic |
| exon_peakiness | Controls the "spikyness" of the data in exonic windows. A smaller value will result in a less even window distribution with higher maxima |
| intron_peakiness | Controls the "spikyness" of the data in intronic windows. A smaller value will result in a less even window distribution with higher maxima |
| intron_exon_balance | Controls the proportion of reads that are assigned to intronic versus exonic regions. HIGHER values correspond to more reads being assigned to INTRONS while LOWER values correspond to more EXON reads. Scale is from 0.0 to 1.0  |
| mean_exon_windows | Mean number of consecutive exonic windows (controls simulated exon lengths) |
| mean_intron_windows | Mean number of consecutive intronic windows (controls simulated intron lengths) |
| window_noise_sd | Controls the variance of the splitting behavior across replicates |

#### binding/pulldown simulation parameters. 
| Parameter    | Description |
| -------- | ------- |
| zeta_mean | Controls the "direction" of the treatment effect on binding. For example, a value of 1 will result in a net increase in binding, while a value of -1 will result in a net decrease in binding. Set this value to 0 for a balanced change in binding dependent only on zeta_sd |
| zeta_sd | Controls the deviation of the treatments effect on binding across windows |
| zeta_rep_sd | Controls the deviation of the treatments effect on binding across replicates | 
| Mean_b | Controls the mean of the gamma distribution from which the binding effect will be drawn from. The binding effect is multiplied by the total expression of the gene in order to get the pulldown value for a specific window |
| Shape_b | Controls the shape parameter of the gamma distribution from which the mean binding effect will be drawn  |
| Min_b | Controls the minimum binding effect. This value is important as we do not want binding/pulldown to decrease the number of reads in a window |
| pulldown_a | Controls the dispersion of the mean pulldown between different samples |
|pulldown_b | Controls the dispersion of the mean pulldown between different samples |
| ip_size_factors | A vector of size factors for the IP pulldown. This is seperate from the expression size factors, and is meant for simulating major changes in pulldown efficiency between samples rather than changes in read-depth |
| diff_frac | Controls the fraction of windows affected by the differential treatment |

**IMPORTANT NOTE:** Since IP size factors are only applied to IP samples, the ip size factors are applied in the following order. 
| ctrl IP 1 | ctrl IP 2 | trt IP 1 | trt IP 2 |
| -------- | ------- | ------- | ------- |
| | | | |

Once all parameters for the simulation have been selected, you can simply run all cells of the notebook to complete the simulation. 

### Running differential eCLIP analysis and validation (notebooks 2-5). 

`2_differential_skipper.ipynb`, `3_run_skipper.ipynb`, `4_run_Flipper_on_sim.ipynb`, and `5_sim_validate.ipynb` do not require any additional parameters to run. These notebooks can be executed as is, in order, to perform differential analysis on the simulated data and calculate performance metrics.

## comparison_to_real_data
This notebook generates the comparisons between real and simulated eCLIP data shown in figure 5 and supplemental figures 2 and 3 of [the Flipper manuscript](https://doi.org/10.64898/2026.03.13.711628). The real data used is from the NONO dataset described in the manuscript and originally generated in [https://doi.org/10.1038/s41589-023-01270-0](https://doi.org/10.1038/s41589-023-01270-0). A copy of this data including only the samples/metadata needed to run the comparisons notebook is provided in the `NONO_data` folder. 

## demonstrate_hier_norm_need.ipynb
This notebook generates plots used in supplemental figure 1 of [the Flipper manuscript](https://doi.org/10.64898/2026.03.13.711628).

# Simulation_pipeline

## Install Snakemake 9

If you need Conda, follow the [Miniconda installation guide](https://www.anaconda.com/docs/getting-started/miniconda/install).
Create the same Snakemake environment used by Flipper:

```bash
conda create -n snakemake9 -c conda-forge -c bioconda snakemake=9.12.0
```

Snakemake creates the simulation's Conda environment on the first run.

## Choose a pipeline config

- [manuscript.yaml](simulation_pipeline/configs/manuscript.yaml) runs expression-effect SDs 0 and 2
  for five repetitions each and writes to `simulation_pipeline/results/`.
- [test.yaml](simulation_pipeline/configs/test.yaml) runs one 800-gene simulation with 150 bound
  sites and writes to `simulation_pipeline/test_results/`.

The batch runs each value in `BETA_SD` for `SIM_REPEATS` independent repetitions. Set the core model
parameters in the Simulation parameters section. Relative values of `RUN_DIR` are resolved from the directory
containing the Snakefile. Use an absolute `RUN_DIR` to write the batch to scratch.
Snakemake finds the R scripts and bundled reference data relative to the Snakefile.

## Run on Slurm

Install the Slurm executor in the Snakemake environment:

```bash
conda install -n snakemake9 -c conda-forge -c bioconda snakemake-executor-plugin-slurm=1.4.0
```

Set the Conda prefix, Slurm account, partition, and QoS for your cluster in
[`simulation_pipeline/profiles/config.yaml`](simulation_pipeline/profiles/config.yaml).
Set the `--configfile` line in the launcher to `configs/test.yaml` or `configs/manuscript.yaml`.
From the repository root, run:

```bash
conda activate snakemake9
cd simulation_pipeline
bash run_scaled_simulation.sh
```

## Run locally

From the repository root, activate the same environment and choose a config explicitly:

```bash
conda activate snakemake9
cd simulation_pipeline
snakemake -s snakemake_sim_scale.py --configfile configs/test.yaml --cores 1 --use-conda
```

The benchmark notebooks use the manuscript run names and paths. If you change the output directory
or scenarios, update the notebook paths and scenario settings before analyzing the results.
