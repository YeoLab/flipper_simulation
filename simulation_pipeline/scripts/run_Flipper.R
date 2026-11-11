############################### Setup ##############################

# Load necessary packages.
suppressPackageStartupMessages({
    library(glue)
    library(tidyverse)
    library(Matrix)
    library(tictoc)
    library(parallel)
    library(tibble)
    library(data.table)
    library(DESeq2)
    library(data.table)
    library(purrr)
    library(BiocManager)
})

# Load arguments
args <- commandArgs(trailingOnly = TRUE)
output_dir = args[1]
utils_path = args[2]
sample_ctrl = args[3]
sample_trt = args[4]

# Load up functions from utils. 
source(utils_path)

############################### Preprocessing ##############################

# Load up count data. 
count_data = read.csv(paste0(output_dir,"/window_counts/window_counts.tsv"),sep = "\t")

ctrl = sample_ctrl
trt = sample_trt

# Load up the window data. 
window_data = load_and_combine_window(output_dir, sample_ctrl, sample_trt)

# Get the individual count data from each rep for each window. 
window_counts_INw = make_INw(ctrl, trt, window_data, count_data, output_dir) 

# Calculate genewise IN by summing over gene name in the count data. 
feature_wise_IN = get_feature_wise(count_data, ctrl, trt)$feature_wise_IN

# Calculate genewise IN by summing over gene name in the count data. 
gene_wise = get_gene_wise(count_data, ctrl, trt)
gene_wise_IN = gene_wise$gene_wise_IN
gene_wise_IP = gene_wise$gene_wise_IP

# Replace the window wise IN counts with the gene wise IN counts and save. 
window_counts = replace_window_IN_count(window_counts_INw, feature_wise_IN, ctrl, trt, paste0(output_dir, "/flipper"))

############################### Normalization ##############################

# A weak filter is required to avoid potential errors from many zero count rows. 
ctrl_cols  <- get_reps(count_data, ctrl)  
trt_cols  <- get_reps(count_data, trt)  
just_counts <- c(ctrl_cols, trt_cols)
filtered_counts = count_data[rowMeans(count_data[just_counts]) > 1,]

# Load up the reproducible windows for each treatment group. 
IP_ctrl = get_window_counts(filtered_counts, output_dir, sample_ctrl, ctrl)
IP_trt = get_window_counts(filtered_counts, output_dir, sample_trt, trt)
window_ctrl = IP_ctrl$window_counts
window_trt = IP_trt$window_counts
background_ctrl = IP_ctrl$background_counts
background_trt = IP_trt$background_counts

# Get the rep counts for later.
reps_trt <- as.integer(
  sub(".*_(IN|IP)_", "", grep(paste0("^", trt, "_"), names(filtered_counts), value = TRUE))
) |> unique()
reps_ctrl <- as.integer(
  sub(".*_(IN|IP)_", "", grep(paste0("^", ctrl, "_"), names(filtered_counts), value = TRUE))
) |> unique()

# Calculate size factors for the input. 
IN_sf = get_sf_MOR(gene_wise_IN[, setdiff(names(gene_wise_IN), c("gene_name","length", "gc", "window_name"))])

# Calculate norm factors for the treatment and ctrl IPs seperately. 
IP_ctrl_sf <- get_sf_MOR(window_ctrl[, setdiff(names(window_ctrl), c("gene_name", "window_name"))])
IP_trt_sf <- get_sf_MOR(window_trt[, setdiff(names(window_trt), c("gene_name","window_name"))])

# Combine into one vector.
temp_norm_factors <- c(IP_trt_sf, IP_ctrl_sf, IN_sf)

# Build the vector of column names for control IP replicates
ctrl_cols = names(IP_ctrl_sf)

# Same for treatment IP replicates
trt_cols  = names(IP_trt_sf)

# Keep only background windows present in both conditions. 
background_both = inner_join(background_ctrl, background_trt, by = "window_name", suffix = c("_ctrl", "_trt"))

# Sum them as new columns.
background_both[[paste0(ctrl, "_IP_sum")]] = rowSums(background_both[paste0(ctrl_cols)], na.rm = TRUE)
background_both[[paste0(trt, "_IP_sum")]]  = rowSums(background_both[paste0(trt_cols)],  na.rm = TRUE)

# Collapse background windows to genes.
background_gene = background_both %>%
    group_by(gene_name_ctrl) %>%
    summarise(
        across(all_of(c(paste0(ctrl, "_IP_sum"), paste0(trt, "_IP_sum"))), ~ sum(.x, na.rm = TRUE)),
        .groups = "drop"
    ) %>%
    dplyr::rename(gene_name = gene_name_ctrl)

# Size factors for average.
ctrl_sum_col = paste0(ctrl, "_IP_sum")
trt_sum_col  = paste0(trt, "_IP_sum")
IP_backg_sf = get_sf_MOR(dplyr::select(background_gene, gene_name, all_of(c(ctrl_sum_col, trt_sum_col))))

# Combine the rep specific size factors with the average background size factors. 
norm_factors = combine_ip_size_factors(temp_norm_factors, IP_backg_sf, ctrl, trt)

############################### Differential analysis ##############################

# Define the columns you want to keep in the order you would like them. 
ctrl_cols  = get_reps(window_counts, ctrl)
trt_cols  = get_reps(window_counts, trt)   
select_cols = c(ctrl_cols, trt_cols)

# Run the differential analysis. 
result = process_differential(window_counts, norm_factors, ctrl, trt, select_cols, 0.05, 0.1) # Keep log fold change cutoff low to not interfere with results. 
groups=c("IP","Trt") # for PCA analysis

# Save results. 
write.table(result$result_check, file = paste0(output_dir, "/flipper/tables/all_windows.tsv"), sep="\t", row.names = FALSE, quote = FALSE)
write.table(result$sigs, file = paste0(output_dir, "/flipper/tables/differential_windows.tsv"), sep="\t", row.names = FALSE, quote = FALSE)

writeLines("OK", file.path(output_dir, "/flipper/Flipper.done"))
