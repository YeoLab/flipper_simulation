######################### Quality of life #########################
add_x_prefix = function(cond) {
  if (grepl("^[0-9]", cond)) {
    return(paste0("X", cond))
  } else {
    return(cond)
  }
}

outersect = function(x, y) {
  sort(union(setdiff(x, y), setdiff(y, x)))
}

get_reps <- function(df, condition) {
  grep(paste0("^", condition, "_(IN|IP)_"), names(df), value = TRUE)
}

get_reps_IN_only <- function(df, condition) {
  grep(paste0("^", condition, "_IN_"), names(df), value = TRUE)
}

get_reps_IP_only <- function(df, condition) {
  grep(paste0("^", condition, "_IP_"), names(df), value = TRUE)
}

######################### For preprocessing #########################
# automatically extracts out the condition names from the count data. 
extract_condition = function(count_data) {
  colnames_data = colnames(count_data)
  count_cols = colnames_data[8:length(colnames_data)]
  base_names = gsub("_IN_\\d+|_IP_\\d+", "", count_cols)
  cond_name = unique(base_names)
  if (length(cond_name) != 1) {
    stop(glue("Ambiguous or multiple condition names detected: {toString(cond_name)}"))
  }
  return(cond_name)
}

load_and_combine_skipper_count = function(skipper_dir, sample_ctrl, sample_trt, anno_path) {
    
    # Load the complete count data. 
    count_ctrl = read.csv(paste0(skipper_dir, '/secondary_results/counts/genome/tables/',sample_ctrl, '.tsv.gz'),
                              sep = '\t')
    count_trt = read.csv(paste0(skipper_dir, '/secondary_results/counts/genome/tables/',sample_trt, '.tsv.gz'),
                              sep = '\t')

    # Extract out the ctrl and treatment names. 
    ctrl = extract_condition(count_ctrl)
    trt = extract_condition(count_trt)

    # Create a list of the ctrl columns of interest.  
    ctrl_cols  <- get_reps(count_ctrl, ctrl)   
    select_col <- c("chr","start","end","strand","gc","name", ctrl_cols)

    # Subset the count data to the columns of interest. 
    count_data = subset(count_ctrl, select=select_col)

    # Create a list of the treatment columns of interest. 
    trt_cols = get_reps(count_trt, trt)   
    
    # Add the count data 2 to count data 1
    count_data[trt_cols] = count_trt[trt_cols]

    # Add in length and gene id data. 
    annotations = read_tsv(anno_path)
    if (length(annotations$start) != length(count_data$start)) {
        stop(
            paste0(
                "Error: Mismatch in lengths.\n",
                "annotations: ", length(annotations$start), "\n",
                "count_data: ", length(count_data$start)
            )
        )
    }
    count_data$gene_id = annotations$gene_id
    count_data$base_feature_id = annotations$base_feature_id
    count_data$length = count_data$end - count_data$start
    
    return(list(count_data = count_data, ctrl = ctrl, trt = trt))
}

load_and_combine_window = function(skipper_dir, sample_ctrl, sample_trt) {

    # Load up the reproducible windows for each treatment group. 
    window_ctrl = read.csv(glue('{skipper_dir}/reproducible_enriched_windows/{sample_ctrl}.reproducible_enriched_windows.tsv.gz'),
                              sep = '\t')
    window_trt = read.csv(glue('{skipper_dir}/reproducible_enriched_windows/{sample_trt}.reproducible_enriched_windows.tsv.gz'),
                              sep = '\t')

    # Combine the ctrl and treatment windows. 
    window_data = rbind(window_ctrl,window_trt)

    # Remove duplicated windows (no need to double count windows found in both). 
    window_data <- window_data[!duplicated(window_data$window_name), ]
    return(window_data)
}

# Get the individual counts at every window. 
make_INw = function(ctrl, trt, window_data, count_data, output_dir) {
    
    # Get the counts for all of the enriched window. 
    window_counts_INw = left_join(window_data, count_data)

    ctrl_cols  <- get_reps(window_counts_INw, ctrl)   
    trt_cols  <- get_reps(window_counts_INw, trt) 
    
    # establish a consistent column order for all datasets. 
    col_order = c( 'window_name','feature_id',ctrl_cols, trt_cols)
    window_counts_INw = window_counts_INw[, col_order]
    
    return(window_counts_INw)
}

get_feature_wise = function(count_data, ctrl, trt) {
    # Sum up 
    feature_wise_counts <- count_data %>%
      group_by(feature_id) %>%
      summarise(
        across(
          .cols = where(is.numeric) & !matches("^gc$"),
          .fns  = ~ sum(.x, na.rm = TRUE)
        ),
        gc = mean(gc, na.rm = TRUE)  # special case
      ) %>%
      ungroup()
    
    # Remove these broken columns. 
    feature_wise_counts$start = NULL
    feature_wise_counts$end = NULL

    # Subset to just the IN data. 
    trt_IN = get_reps_IN_only(feature_wise_counts, trt)
    ctrl_IN =  get_reps_IN_only(feature_wise_counts, ctrl)
    IN_subset <- c("feature_id",ctrl_IN, trt_IN)

    # Subset to just the IN columns and save. 
    feature_wise_IN = feature_wise_counts[, IN_subset]

    # Subset to just the IP data. 
    trt_IP = get_reps_IP_only(feature_wise_counts, trt)
    ctrl_IP =  get_reps_IP_only(feature_wise_counts, ctrl)
    IP_subset <- c("feature_id", ctrl_IP, trt_IP)

    # Subset to just the IP columns and save. 
    feature_wise_IP = feature_wise_counts[, IP_subset]

    return(list(feature_wise_IN = feature_wise_IN, feature_wise_IP = feature_wise_IP))
}

get_gene_wise = function(count_data, ctrl, trt) {
    # Sum up 
    gene_wise_counts <- count_data %>%
      group_by(gene_name) %>%
      summarise(
        across(
          .cols = where(is.numeric) & !matches("^gc$"),
          .fns  = ~ sum(.x, na.rm = TRUE)
        ),
        gc = mean(gc, na.rm = TRUE)  # special case
      ) %>%
      ungroup()
    
    # Remove these broken columns. 
    gene_wise_counts$start = NULL
    gene_wise_counts$end = NULL

    # Subset to just the IN data. 
    trt_IN = get_reps_IN_only(gene_wise_counts, trt)
    ctrl_IN =  get_reps_IN_only(gene_wise_counts, ctrl)
    IN_subset <- c("gene_name",ctrl_IN, trt_IN)

    # Subset to just the IN columns and save. 
    gene_wise_IN = gene_wise_counts[, IN_subset]

    # Subset to just the IP data. 
    trt_IP = get_reps_IP_only(gene_wise_counts, trt)
    ctrl_IP =  get_reps_IP_only(gene_wise_counts, ctrl)
    IP_subset <- c("gene_name", ctrl_IP, trt_IP)

    # Subset to just the IP columns and save. 
    gene_wise_IP = gene_wise_counts[, IP_subset]

    return(list(gene_wise_IN = gene_wise_IN, gene_wise_IP = gene_wise_IP))
}

# Replace the window wise IN counts (INw) with gene wise IN counts (INf)
replace_window_IN_count = function(window_counts_INw, feature_wise_IN, ctrl, trt, output_dir) {

    # define the columns to keep for the final dataset. 
    trt_cols = get_reps(window_counts_INw, trt)
    ctrl_cols = get_reps(window_counts_INw, ctrl)
    col_order = c('window_name','feature_id', trt_cols, ctrl_cols)

    # Get a list of the INw columns to drop. 
    trt_cols = get_reps_IN_only(window_counts_INw, trt)
    ctrl_cols = get_reps_IN_only(window_counts_INw, ctrl)
    drops = c(trt_cols, ctrl_cols)
    
    # replace the window-wise IN counts with the gene-wise IN counts.
    temp = window_counts_INw[ , !(names(window_counts_INw) %in% drops)]
    window_counts_INf = left_join(temp, feature_wise_IN, by = "feature_id")
    window_counts_INf = window_counts_INf[, col_order]
    write.table(window_counts_INf, glue("{output_dir}/intermediates/window_counts_INf.tsv"), sep="\t", row.names = FALSE)
    return(window_counts_INf)
}

######################### For size factors #########################

get_window_counts = function(count_data, skipper_dir, sample_cond, cond) {
    windows = read.csv(paste0(skipper_dir, '/reproducible_enriched_windows/',sample_cond,'.reproducible_enriched_windows.tsv.gz'),
                              sep = '\t')
    trt_cols = get_reps_IP_only(count_data, cond)   
    IP_cols <- c(trt_cols, "gene_name", "window_name")
    window_counts = dplyr::select(left_join(windows, count_data), all_of(IP_cols))
    background_counts = dplyr::select(anti_join(count_data, windows), all_of(IP_cols))
    return(list(window_counts = window_counts, background_counts = background_counts))
}

# Rescale a size-factor vector so its mean equals 1. 
rescale_to_mean <- function(x, target = 1) {
  x * (target / mean(x, na.rm = TRUE))
}

# Rescale a numeric data frame (ignoring "gene_id") so its total sum equals `target`
rescale_df_total <- function(df, target) {
  s <- mean(as.matrix(df), na.rm = TRUE) 
  scaled_num <- df * (target / s)
  return(scaled_num)
}

combine_ip_size_factors <- function(norm_factors, IP_backg_sf, ctrl, trt, ip_pattern = "_IP_") {
    # Helper for geometric mean.
    GM <- function(x) exp(mean(log(x)))
    
    # Get the sample/replicate names.
    nms <- names(norm_factors)
    
    # Find IP replicate columns for ctrl / trt.
    ctrl_ip_cols <- grep(paste0("^", ctrl, ip_pattern), nms, value = TRUE)
    trt_ip_cols  <- grep(paste0("^", trt,  ip_pattern), nms, value = TRUE)
    
    # Pull summed IP size factors.
    ctrl_sum_name <- paste0(ctrl, ip_pattern, "sum")
    trt_sum_name  <- paste0(trt,  ip_pattern, "sum")
    
    ctrl_sum_sf <- IP_backg_sf[ctrl_sum_name]
    trt_sum_sf  <- IP_backg_sf[trt_sum_name]
    
    # Target ratio from summed libraries.
    R <- as.numeric(ctrl_sum_sf / trt_sum_sf)

    if (is.null(dim(norm_factors))) {
        d <- GM(norm_factors[ctrl_ip_cols])
        o <- GM(norm_factors[trt_ip_cols])

        # Scaling factors.
        a <- sqrt(R * (o / d))
        b <- 1 / a

        scaled <- norm_factors
        scaled[ctrl_ip_cols] <- a * scaled[ctrl_ip_cols]
        scaled[trt_ip_cols]  <- b * scaled[trt_ip_cols]        

    } else {
        d <- as.numeric(GM(colMeans(norm_factors[, ctrl_ip_cols, drop = FALSE])))
        o <- as.numeric(GM(colMeans(norm_factors[, trt_ip_cols,  drop = FALSE])))

        # Scaling factors.
        a <- sqrt(R * (o / d))
        b <- 1 / a

        scaled <- norm_factors
        scaled[, ctrl_ip_cols] <- a * scaled[, ctrl_ip_cols, drop = FALSE]
        scaled[, trt_ip_cols]  <- b * scaled[, trt_ip_cols,  drop = FALSE]        
        }
    return(scaled)
}

get_sf_MOR = function(data) {
    # Pick only count columns. 
    count_cols <- setdiff(names(data),
                          c("gene_name","window_name"))
    
    # Coerce to plain numeric matrix.
    counts_df <- as.data.frame(data[, count_cols, drop = FALSE])
    counts_df[] <- lapply(counts_df, function(x) as.numeric(as.character(x)))
    counts <- as.matrix(counts_df)
    
    # Set gene rownames. 
    rownames(counts) <- as.character(data$gene_name)
    
    # Remove rows with any NA. 
    ok <- complete.cases(counts)
    counts <- counts[ok, , drop = FALSE]
    storage.mode(counts) <- "numeric"  # ensure numeric
    
    # Get size factors. 
    sf <- DESeq2::estimateSizeFactorsForMatrix(counts, type = "poscounts")
    return(sf)
}
                          
######################### For differential analysis #########################
                          
remove_na_rows <- function(sf, counts) {

  # Find complete rows in each dataframe.
  complete_sf   <- complete.cases(sf)
  complete_counts <- complete.cases(counts)
  
  # Keep rows that are complete in both.
  keep <- complete_sf & complete_counts
  
  return(list(sf = sf[keep, , drop = FALSE], counts = counts[keep, , drop = FALSE]))
}

process_differential = function(window_counts, sf, ctrl, trt, select_cols, pval_threshold, log2fc_threshold){

    if (is.null(dim(sf))) {
        wc_sub = tidyr::drop_na(window_counts)
        wc_mat = as.matrix(dplyr::select(wc_sub, all_of(select_cols)))

        # Ensure proper ordering. 
        sf = sf[select_cols]
        wc_mat = wc_mat[, select_cols, drop = FALSE]
    } else {
        # Decide which windows to keep (shared between the two tables).
        keep_windows = intersect(sf$name, window_counts$name)
        
        # Reorder both tables to the same window order, then subset to the same sample columns.
        sf_sub = sf[match(keep_windows, sf$name), colnames(sf), drop = FALSE]
        wc_sub = window_counts[match(keep_windows, window_counts$name), colnames(window_counts), drop = FALSE]
        
        # Convert to matrices.
        sf_mat = as.matrix(dplyr::select(sf_sub, all_of(select_cols)))
        wc_mat = as.matrix(dplyr::select(wc_sub, all_of(select_cols)))
        
        # Remove rownames
        rownames(sf_mat) = NULL
        rownames(wc_mat) = NULL

        # Ensure proper ordering. 
        sf_mat = sf_mat[, select_cols, drop = FALSE]
        wc_mat = wc_mat[, select_cols, drop = FALSE]
    }

    # Define the model and design for deseq. 
    sample_names = colnames(wc_mat)
    model = data.frame(
      Reps = paste0("Rep", sub(".*_(IN|IP)_", "", sample_names)),
      IP   = ifelse(grepl("_IN_", sample_names), "Input",
                    ifelse(grepl("_IP_", sample_names), "IP", NA)),
      Trt  = ifelse(grepl(paste0("^", trt, "_"),  sample_names), "Trt",
                    ifelse(grepl(paste0("^", ctrl, "_"), sample_names), "Ctrl", NA)),
      row.names = sample_names,
      stringsAsFactors = FALSE
    )
    design = ~Reps + IP*Trt

    # Create the deseq object. 
    dds = DESeqDataSetFromMatrix(countData = wc_mat,
                                colData = model,
                                design = design)

    # If sf is a vector, treat it as size factors. If a matrix, as normalization factors.
    if (is.null(dim(sf))) {
        sizeFactors(dds) = sf
    } else {
        normalizationFactors(dds) = sf_mat # Note that the offset alterations reccomended by EDAseq were already performed in calc_norm. 
    }

    # Run deseq.
    dds = DESeq2::estimateDispersions(dds, fitType = "local")
    dds = DESeq2::nbinomWaldTest(dds)

    # Extract deseq results into a cleaned table.
    res = DESeq2::results(dds, name = "IPIP.TrtTrt") 
    dat = data.frame(padj = res$padj, pval = res$pvalue, stat = res$stat, log2FoldChange = res$log2FoldChange)

    # Extract deseq results into a cleaned table. 
    result_check = cbind(wc_sub, dat)
    
    # Sort be the adjusted p-value. 
    result_check = result_check %>% arrange(padj)

    # Extract out the significant results. 
    sigs = dplyr::filter(result_check, padj < pval_threshold, abs(log2FoldChange) > log2fc_threshold)

    # Return the final results. 
    return(list(dds = dds, window_counts = window_counts, wc_mat = wc_mat, result_check = result_check, sigs = sigs))
}