import json
from pathlib import Path

# Resolve the output directory from the config passed to Snakemake.
pipeline_dir = Path(workflow.basedir)
run_dir = pipeline_dir / Path(config["RUN_DIR"]).expanduser()
workdir: str(run_dir)

# Pass the model settings to R and keep batch settings in Snakemake.
batch_keys = ("RUN_DIR", "BETA_SD", "SIM_REPEATS", "SEED")
simulation_settings = {key: value for key, value in config.items() if key not in batch_keys}

# Name each simulation by its expression-effect SD and repetition.
sims = {
    f"beta_{str(beta_sd).removesuffix('.0')}_rep{rep}": (beta_sd, rep)
    for beta_sd in config["BETA_SD"]
    for rep in range(1, config["SIM_REPEATS"] + 1)
}

rule all:
    input:
        expand("{sim}/init_sim.done", sim=sims),
        expand("{sim}/skipper.done", sim=sims),
        expand("{sim}/diffskip.done", sim=sims),
        expand("{sim}/flipper/Flipper.done", sim=sims)

# Rule for running the simulation. 
rule init_sim:
    input:
        r_script = workflow.source_path("scripts/run_simulation.R"),
        utils = workflow.source_path("../tools/simulation_utils.r"),
        gene_lengths = workflow.source_path("../NONO_data/lens_bp.rds")
    output:
        done = "{sim}/init_sim.done",
        background = "{sim}/gene_counts/background_count.tsv",
        metrics = "{sim}/true_values/binding_metrics.tsv",
        window_counts = "{sim}/window_counts/window_counts.tsv"
    params:
        beta_sd = lambda wc: sims[wc.sim][0],
        seed = lambda wc: config["SEED"] + sims[wc.sim][1] - 1,
        settings = json.dumps(simulation_settings, sort_keys=True)
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
        Rscript --vanilla {input.r_script:q} {wildcards.sim:q} {params.beta_sd} \
          {input.utils:q} {input.gene_lengths:q} {params.seed} {params.settings:q} &> {log:q}
        """

# Rule for running normal skipper on all outputs. 
rule run_skipper:
    input:
        r_script = workflow.source_path("scripts/run_skipper.R"),
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
        Rscript --vanilla {input.r_script:q} \
          {wildcards.sim:q} \
          {input.window_counts:q} \
        &> {log:q}
        """

# Rule for running differential skipper on all outputs. 
rule run_diffskip:
    input:
        r_script = workflow.source_path("scripts/run_diffskip.R"),
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
        Rscript --vanilla {input.r_script:q} \
          {wildcards.sim:q} \
          {input.window_counts:q} \
        &> {log:q}
        """

# Rule for running Flipper on all outputs. 
rule run_Flipper:
    input:
        r_script = workflow.source_path("scripts/run_Flipper.R"),
        utils = workflow.source_path("../tools/Flipper_utils.r"),
        window_counts = "{sim}/window_counts/window_counts.tsv",
        rep_ctrl = "{sim}/reproducible_enriched_windows/control.reproducible_enriched_windows.tsv.gz",
        rep_trt  = "{sim}/reproducible_enriched_windows/treatment.reproducible_enriched_windows.tsv.gz",
    output:
        done = "{sim}/flipper/Flipper.done",
        all_windows = "{sim}/flipper/tables/all_windows.tsv",
        diff_windows = "{sim}/flipper/tables/differential_windows.tsv",
        counts_INf = "{sim}/flipper/intermediates/window_counts_INf.tsv",
    params:
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
        Rscript --vanilla {input.r_script:q} \
          {wildcards.sim:q} \
          {input.utils:q} \
          {params.sample_ctrl} \
          {params.sample_trt} \
        &> {log:q}
        """