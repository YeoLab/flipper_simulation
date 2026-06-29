# Import necessary packages. 
import re
import os

SIMS = ["beta_0_rep1","beta_0_rep2","beta_0_rep3","beta_0_rep4","beta_0_rep5",
        "beta_2_rep1","beta_2_rep2","beta_2_rep3","beta_2_rep4","beta_2_rep5"]

# Create each directory and an additional folder to hold the log data. 
for sim in SIMS:
    os.makedirs(f"{sim}/logs", exist_ok=True)

# Create a small function to extract the appropriate parameters from each 
def parse_sim_name(sim_name: str):
    m = re.match(r"^beta_(-?\d+(?:\.\d+)?)_rep\d+$", sim_name)
    beta_sd = float(m.group(1))
    return beta_sd

# Create rule all (Calls to all other rules. )
rule all:
    input:
        expand("{sim}/init_sim.done", sim=SIMS),
        expand("{sim}/skipper.done", sim=SIMS),
        expand("{sim}/diffskip.done", sim=SIMS),
        expand("{sim}/flipper/Flipper.done", sim=SIMS)

# Rule for running the simulation. 
rule init_sim:
    output:
        done = "{sim}/init_sim.done",
        background = "{sim}/gene_counts/background_count.tsv",
        metrics = "{sim}/true_values/binding_metrics.tsv",
        window_counts = "{sim}/window_counts/window_counts.tsv"
    params:
        beta_sd = lambda wc: parse_sim_name(wc.sim),
        utils = "../../tools/simulation_utils.r"
    resources:
        mem_mb = 16000,
        runtime = 60,
    log:
        "{sim}/logs/init_sim.log"
    conda:
        "envs/environment.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla ../scripts/run_simulation.R \
          {wildcards.sim} \
          {params.beta_sd} \
          {params.utils} \
        &> {log}
        """

# Rule for running normal skipper on all outputs. 
rule run_skipper:
    input:
        window_counts = "{sim}/window_counts/window_counts.tsv",
    output:
        done = "{sim}/skipper.done",
        rep_control_gz = "{sim}/reproducible_enriched_windows/control.reproducible_enriched_windows.tsv.gz",
        rep_treat_gz   = "{sim}/reproducible_enriched_windows/treatment.reproducible_enriched_windows.tsv.gz",
    resources:
        mem_mb = 24000,
        runtime = 480,
    log:
        "{sim}/logs/run_skipper.log"
    conda:
        "envs/environment.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla ../scripts/run_skipper.R \
          {wildcards.sim} \
          {input.window_counts} \
        &> {log}
        """

# Rule for running differential skipper on all outputs. 
rule run_diffskip:
    input:
        window_counts = "{sim}/window_counts/window_counts.tsv",
    output:
        done = "{sim}/diffskip.done",
        rep_ctrtrt = "{sim}/diff_skipper/reproducible_enriched_windows/ctr-trt.reproducible_enriched_windows.tsv",
        rep_trtctr = "{sim}/diff_skipper/reproducible_enriched_windows/trt-ctr.reproducible_enriched_windows.tsv",
    resources:
        mem_mb = 24000,
        runtime = 480,
    log:
        "{sim}/logs/run_diffskip.log"
    conda:
        "envs/environment.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla ../scripts/run_diffskip.R \
          {wildcards.sim} \
          {input.window_counts} \
        &> {log}
        """

# Rule for running Flipper on all outputs. 
rule run_Flipper:
    input:
        window_counts = "{sim}/window_counts/window_counts.tsv",
        rep_ctrl = "{sim}/reproducible_enriched_windows/control.reproducible_enriched_windows.tsv.gz",
        rep_trt  = "{sim}/reproducible_enriched_windows/treatment.reproducible_enriched_windows.tsv.gz",
    output:
        done = "{sim}/flipper/Flipper.done",
        all_windows = "{sim}/flipper/tables/all_windows.tsv",
        diff_windows = "{sim}/flipper/tables/differential_windows.tsv",
        counts_INf = "{sim}/flipper/intermediates/window_counts_INf.tsv",
    params:
        utils = "../../tools/Flipper_utils.r", 
        sample_ctrl = "control",
        sample_trt = "treatment"
    resources:
        mem_mb = 24000,
        runtime = 60
    log:
        "{sim}/logs/run_Flipper.log"
    conda:
        "envs/environment.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla ../scripts/run_Flipper.R \
          {wildcards.sim} \
          {params.utils} \
          {params.sample_ctrl} \
          {params.sample_trt} \
        &> {log}
        """