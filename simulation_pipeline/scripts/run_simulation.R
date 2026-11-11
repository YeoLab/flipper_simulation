############################ Setup #############################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(glue)
  library(tidyverse)
  library(truncnorm)
  library(readr)
  library(glue)
  library(tictoc)
  library(Matrix)
  library(miceadds)
  library(matrixStats)
  library(parallel)
  library(tibble)
  library(RColorBrewer)
  library(ggrepel)
  library(data.table)
})

# Read the resolved settings and seed this simulation.
args = commandArgs(trailingOnly = TRUE)
output_dir = args[1]
beta_sd = as.numeric(args[2])
source(args[3])
lens_bp = readRDS(args[4])
set.seed(as.integer(args[5]))
invisible(list2env(jsonlite::fromJSON(args[6]), envir = environment()))

# Set up the two-condition count model.
reps = 2
n_samples = reps * 4
dispMeanRel = function(x) dispersion_scale / x + dispersion_floor

# Create all necessary directories. 
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "window_counts"), showWarnings = FALSE)
dir.create(file.path(output_dir, "gene_counts"), showWarnings = FALSE)
dir.create(file.path(output_dir, "input_model_coef"), showWarnings = FALSE)
dir.create(file.path(output_dir, "enriched_windows"), showWarnings = FALSE)
dir.create(file.path(output_dir, "reproducible_enriched_windows"), showWarnings = FALSE)
dir.create(file.path(output_dir, "true_values"), showWarnings = FALSE)
dir.create(file.path(output_dir, "diff_skipper"), showWarnings = FALSE)
dir.create(file.path(output_dir, "diff_skipper/input_model_coef"), showWarnings = FALSE)
dir.create(file.path(output_dir, "diff_skipper/enriched_windows"), showWarnings = FALSE)
dir.create(file.path(output_dir, "diff_skipper/reproducible_enriched_windows"), showWarnings = FALSE)
dir.create(file.path(output_dir, "DIBS/tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "DIBS/intermediates"), showWarnings = FALSE)

############################ Run the simulation #############################

# Simulate the initial reads. 
out = simulate_gene_counts_with_length(n_genes = n_genes, n_samples = n_samples, beta_shift = beta_shift,
                                        beta_sd = beta_sd, Mean_ex = Mean_ex, SD_ex = SD_ex,
                                        dispMeanRel = dispMeanRel, lens_bp = lens_bp,
                                        sizeFactors = ex_size_factors)
gene_counts = out$gene_counts

# Downsample reads to get background. 
background_counts = downsample_counts(subset(gene_counts, select = -c(gene_name)),
                                      mean_prop = background_fraction)
background_counts$gene_name = gene_counts$gene_name
write_tsv(background_counts, paste0(output_dir, "/gene_counts/background_count.tsv"))

# Split the reads into windows. 
window_counts <- split_counts_into_windows(df = background_counts, gene_length_bp = out$gene_length_bp,
                                           exon_frac = exon_frac, exon_peakiness = exon_peakiness,  
                                           intron_exon_balance = intron_exon_balance,
                                           intron_peakiness = intron_peakiness, mean_exon_windows = mean_exon_windows,
                                           mean_intron_windows = mean_intron_windows, window_noise_sd = window_noise_sd)

# Calculate feature counts from the windows. 
feature_counts = window_counts %>%
group_by(feature_id) %>%
summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
ungroup()

# Calculate the IP counts to add to random windows. 
IP_ex = IPDataSpecZ(feature_counts = feature_counts, n_bind = n_bind,
                    n_samples = n_samples, zetaMean = zeta_mean, zetaSD = zeta_sd, interceptMode = Mean_b, 
                    interceptShape = Shape_b, interceptMin = Min_b, pulldown_a = pulldown_a, pulldown_b = pulldown_b, diff_frac = diff_frac,
                    ip_size_factors = ip_size_factors, zeta_rep_sd = zeta_rep_sd)
IPData = IP_ex$IPData

count_cols <- seq(from = 5, by = 2, length.out = reps*2)

features <- IP_ex$to_bind

# Precompute group ids for all rows once (fast).
gid <- match(window_counts$feature_id, features)   # NA if row not in bound set

# Keep only rows belonging to bound genes.
rows <- which(!is.na(gid))
g <- gid[rows]  # group id in 1..length(gene_names)

# Randomize order once, then take the first row per gene = random pick per gene.
o <- sample.int(length(rows))
rows_shuf <- rows[o]
g_shuf <- g[o]

first_pos <- match(seq_along(features), g_shuf)  # position in shuffled list
subchoices <- rows_shuf[first_pos]

# Safety: if any gene has zero windows, you'll get NA here.
if (anyNA(subchoices)) {
  missing <- features[is.na(subchoices)]
  stop(sprintf("Some bound genes have 0 windows in window_counts: %s",
               paste(missing, collapse = ", ")))
}

# Add pulldown values (IPData row i goes to gene_names[i]'s chosen window).
window_counts_enriched = copy(window_counts)

window_counts_enriched[subchoices, count_cols] <-
  window_counts_enriched[subchoices, count_cols] + IPData

# Add gene and window name information. 
IP_ex$trueVals$feature_id =  IP_ex$to_bind
IP_ex$trueVals$window_name =  subchoices

# Observe the metrics
# Save the metrics
write_tsv(IP_ex$trueVals, paste0(output_dir,"/true_values/binding_metrics.tsv"))

# Simulate GC content for the windows
window_counts_enriched$gc <- rtruncnorm(n = nrow(window_counts_enriched), a = 0, b = 1, mean = 0.4, sd = 0.1)

# Add arbitrary names to all windows
window_counts_enriched$window_name = 1:nrow(window_counts_enriched)

rep_labels_control <- as.vector(rbind(paste0("control_IN_", 1:reps), paste0("control_IP_", 1:reps)))

rep_labels_treatment <- as.vector(rbind(paste0("treatment_IN_", 1:reps), paste0("treatment_IP_", 1:reps)))

colnames(window_counts_enriched) <- c("gene_name", "feature_type", "feature_id", rep_labels_control,
                                      rep_labels_treatment, "gc", "window_name")

# Show and save the final dataframe. 
head(window_counts_enriched)
write_tsv(window_counts_enriched, paste0(output_dir,"/window_counts/window_counts.tsv"))


writeLines("OK", file.path(output_dir, "/init_sim.done"))