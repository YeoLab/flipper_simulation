suppressPackageStartupMessages({
    library(DESeq2)
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
})

# Load arguments.
args <- commandArgs(trailingOnly = TRUE)
output_dir = args[1]
window_counts_path = args[2]
options(dplyr.summarise.inform = FALSE)
window_counts <- readr::read_tsv(window_counts_path, show_col_types = FALSE)
reps = (length(window_counts) - 5)/4

dif_skip = select(window_counts, c('gene_name','control_IP_1','control_IP_2','treatment_IP_1','treatment_IP_2','gc','window_name'))
               
# Define conditions and replicates
for (i in c("ctr-trt", "trt-ctr")) {
    if (i == "ctr-trt") {
        colnames(dif_skip) = c('gene_name','ctr-trt_IN_1','ctr-trt_IN_2','ctr-trt_IP_1','ctr-trt_IP_2','gc','window_name')
    } else {
        colnames(dif_skip) = c('gene_name','trt-ctr_IP_1','trt-ctr_IP_2','trt-ctr_IN_1','trt-ctr_IN_2','gc','window_name')
    }
    for (j in 1:2) {
        given_input_replicate = paste(i, "_IN_", j, sep = "")
        count_data = dif_skip %>% select(matches(glue("^gc$|{i}")))
        other_input_replicates = setdiff(grep("IN_[0-9]+$", names(count_data), value=TRUE), given_input_replicate)
        count_gc_data = count_data[select(count_data, matches("(IP|IN)_[0-9]*$")) %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(gc, 2))
    
        input_betabinom_fit_data = lapply(other_input_replicates, function(other_input_replicate) {
            input_gc_data = count_gc_data[as.logical(count_gc_data[, given_input_replicate] + count_gc_data[, other_input_replicate] > 0), ]
            processed_input_data = input_gc_data %>%
                summarize(
                    given_input_fraction = mean(.data[[given_input_replicate]] / (.data[[other_input_replicate]] + .data[[given_input_replicate]])) %>% (VGAM::logitlink)
                ) %>%
                inner_join(input_gc_data)

            betabinom_fit = suppressWarnings(
                VGAM::vglm(
                    cbind(
                        processed_input_data[[given_input_replicate]],
                        processed_input_data[[other_input_replicate]]
                    ) ~ processed_input_data$given_input_fraction,
                    VGAM::betabinomial,
                    trace = TRUE
                )
            )

            betabinom_coefs = betabinom_fit %>% (VGAM::Coef)

            global_given_input_fraction = sum(input_gc_data[[given_input_replicate]]) / sum(input_gc_data[[other_input_replicate]] + input_gc_data[[given_input_replicate]])
            distribution_data = processed_input_data %>%
                group_by(
                    input_total = .data[[given_input_replicate]] + .data[[other_input_replicate]],
                    .data[[given_input_replicate]],
                    given_input_fraction,
                    gc_bin
                ) %>%
                count(name = "count") %>%
                group_by(input_total, gc_bin) %>%
                mutate(pdf = count / sum(count)) %>%
                rowwise() %>%
                mutate(
                    binomial = dbinom(
                        x = .data[[given_input_replicate]],
                        size = input_total,
                        prob = global_given_input_fraction
                    ),
                    betabinomial = VGAM::dbetabinom(
                        x = .data[[given_input_replicate]],
                        size = input_total,
                        prob = (betabinom_coefs[1] + given_input_fraction * betabinom_coefs[3]) %>% (VGAM::logitlink)(inverse = TRUE),
                        rho = betabinom_coefs[2] %>% (VGAM::logitlink)(inverse = TRUE)
                    )
                )

            betabinom_coefs %>%
                as_tibble(rownames = "coef") %>%
                transmute(coef = c("mu", "rho", "gc_bias_term"), value) %>%
                pivot_wider(names_from = coef, values_from = value) %>%
                mutate(
                    global_given_replicate_fraction = global_given_input_fraction,
                    given_replicate = given_input_replicate,
                    other_replicate = other_input_replicate,
                    experiment = i
                )
        }) %>% bind_rows()

        input_betabinom_fit_data$rho = input_betabinom_fit_data$rho - 1
        write_tsv(input_betabinom_fit_data, paste0(output_dir, "/diff_skipper/input_model_coef/", i, ".", given_input_replicate, ".tsv"))
    }
}

logisticb2 = function(x) 1 / (1 + 2**-x)
logitb2 = function(x) log2(x / (1 - x))

reps = c("1", "2")

for (cond in c("ctr-trt", "trt-ctr")) {
    if (cond == "ctr-trt") {
        colnames(dif_skip) = c('gene_name','ctr-trt_IN_1','ctr-trt_IN_2','ctr-trt_IP_1','ctr-trt_IP_2','gc','window_name')
    } else {
        colnames(dif_skip) = c('gene_name','trt-ctr_IP_1','trt-ctr_IP_2','trt-ctr_IN_1','trt-ctr_IN_2','gc','window_name')
    }
    for (rep in reps) {
        model_data = read_tsv(glue(paste0(output_dir,"/diff_skipper/input_model_coef/{cond}.{cond}_IN_{rep}.tsv")), show_col_types = FALSE)
        input_replicate_label = glue("{cond}_IN_{rep}")
        clip_replicate_label = glue("{cond}_IP_{rep}")
        output_stem = glue("{cond}.{cond}_IP_{rep}")
        
        # Compile GC content information
        count_gc_data = dif_skip[select(dif_skip, matches("(IP|IN)_[0-9]+$")) %>% rowSums > 0, ]
        count_gc_data = group_by(count_gc_data, gc_bin = cut_number(gc, 2))
        count_gc_data = filter(count_gc_data, .data[[clip_replicate_label]] + .data[[input_replicate_label]] > 0)
    
        prior_mass = 5
    
        processed_count_data = count_gc_data %>%
            rename(input = all_of(input_replicate_label),
                   clip = all_of(clip_replicate_label)) %>%
            group_by(gc_bin) %>%
            mutate(
                total_input_gc_bin = sum(input),
                total_clip_gc_bin  = sum(clip),
        
                gc_bin_in_ip_ratio = total_input_gc_bin / total_clip_gc_bin,
        
                gc_bin_input_frac = total_input_gc_bin / (total_input_gc_bin + total_clip_gc_bin),
                gc_bin_clip_frac  = total_clip_gc_bin  / (total_input_gc_bin + total_clip_gc_bin),
        
                # Compute the original symmetric size factors from the GC-bin totals.
                input_sf = sqrt(total_input_gc_bin / total_clip_gc_bin),
                clip_sf  = sqrt(total_clip_gc_bin  / total_input_gc_bin),
        
                # Add a small prior mass split according to the GC-bin background ratio.
                stabilized_input = input + prior_mass * gc_bin_input_frac,
                stabilized_clip  = clip  + prior_mass * gc_bin_clip_frac
            ) %>%
            ungroup() %>%
            mutate(
                original_total = input + clip,
        
                # Apply the original symmetric normalization to the stabilized counts.
                sym_norm_input = stabilized_input / input_sf,
                sym_norm_clip  = stabilized_clip  / clip_sf,
        
                sym_norm_total = sym_norm_input + sym_norm_clip,
        
                # Rescale both sides by a common factor so the total stays near the original total.
                rescale_factor = if_else(sym_norm_total > 0,
                                         original_total / sym_norm_total,
                                         0),
        
                target_norm_input = sym_norm_input * rescale_factor,
                target_norm_clip  = sym_norm_clip  * rescale_factor,
        
                # Round to integer pseudocounts while preserving the original total exactly.
                norm_input = round(target_norm_input),
                norm_clip  = original_total - norm_input
            ) %>%
            select(-matches("(IP|IN)_[0-9]+$"),
                   -gc_bin_input_frac,
                   -gc_bin_clip_frac,
                   -stabilized_input,
                   -stabilized_clip,
                   -original_total,
                   -sym_norm_input,
                   -sym_norm_clip,
                   -sym_norm_total,
                   -rescale_factor,
                   -target_norm_input,
                   -target_norm_clip)
        
        # Set beta-binomial overdispersion (rho) from model, on probability scale. 
        model_overdispersion = VGAM::logitlink(median(model_data$rho), inverse=TRUE)
        
        # Collapse to unique normalized count combinations, compute enrichment and p-values,
        # and rejoin to per-window rows.
        p_data = processed_count_data %>%
            group_by(norm_clip, norm_input, gc_bin) %>%
            summarize(.groups = "drop") %>%
            mutate(
                enrichment_l2or = log2(norm_clip / norm_input),
                pvalue = pmax(
                    1e-12,
                    1 - VGAM::pbetabinom(
                        q = norm_clip - 1,
                        size = norm_clip + norm_input,
                        prob = 0.5,
                        rho = model_overdispersion
                    )
                )
            ) %>%
            inner_join(processed_count_data, ., by = c("norm_clip", "norm_input", "gc_bin"))
    
        p_clip = with(p_data, sum(clip) / sum(clip + input))
        # Optimize coverage threshold for calling enriched windows by 20% FDR
        threshold_max = p_data %>% 
          mutate(total_counts = input + clip) %>% 
          arrange(desc(total_counts)) %>% 
          head(100) %>% 
          tail(1) %>% pull(total_counts)
          
        threshold_data = tibble(
          threshold = seq(2, threshold_max),
          n_enriched = sapply(seq(2, threshold_max), function(i) {
            sum((p_data %>% filter(input + clip >= i) %>% pull(pvalue) %>% p.adjust(., "fdr")) < 0.2)
          })
        )
        
        optimized_threshold = threshold_data %>% 
          arrange(desc(n_enriched)) %>% 
          head(1) %>% pull(threshold)
        
        q_data = p_data %>% 
          group_by(above_threshold = input + clip >= optimized_threshold) %>% 
          mutate(qvalue = ifelse(above_threshold, p.adjust(pvalue, "fdr"), NA)) %>% 
          ungroup
        
        enriched_window_data = q_data %>% 
          filter(above_threshold, qvalue < 0.2) %>% 
          select(-above_threshold) %>% 
          arrange(qvalue)
        write_tsv(enriched_window_data, glue("{output_dir}/diff_skipper/enriched_windows/{output_stem}.enriched_windows.tsv"))
        
    }
}

data_directory = paste0(output_dir, "/diff_skipper/enriched_windows/")
prefixes = c("ctr-trt", "trt-ctr")

for (prefix in prefixes) {
    blacklist = tibble(chr=character(),start=numeric(),end=numeric(),window_name=character(),score=numeric(),strand=character())
    enriched_window_files = list.files(path = data_directory, pattern = paste0("^", prefix, "\\..*enriched_windows.tsv"), full.names = TRUE)
    
    # specify col_type to handle case of no enriched windows
    enriched_window_data = enriched_window_files %>%
        setNames(sub("\\.enriched_windows\\.tsv", "", basename(.))) %>% 
        map(function(x) read_tsv(x, show_col_types = FALSE) %>% mutate(window_name = as.character(window_name)) %>% anti_join(blacklist, by = join_by(window_name))) %>% 
        Filter(function(x) nrow(x) > 0, .) %>%
        bind_rows(.id = "clip_replicate_label") 
    
    reproducible_enriched_window_data = enriched_window_data %>% group_by(across(-c(clip_replicate_label,input,clip,enrichment_l2or,pvalue,qvalue, norm_input, norm_clip,input_sf, clip_sf, gc, gc_bin, total_input_gc_bin, total_clip_gc_bin, gc_bin_in_ip_ratio))) %>%
    	summarize(
    		input_sum = sum(input),
    		clip_sum = sum(clip),
    		enrichment_n = sum(qvalue < 0.2),
    		enrichment_l2or_min = min(enrichment_l2or),
    		enrichment_l2or_mean = mean(enrichment_l2or),
    		enrichment_l2or_max = max(enrichment_l2or),
    		p_max = max(pvalue),
    		p_min = min(pvalue),
    		q_max = max(qvalue),
    		q_min = min(qvalue)
    	) %>% filter(enrichment_n > 1) %>% arrange(enrichment_l2or_mean %>% desc)
    
    write_tsv(reproducible_enriched_window_data, paste0(output_dir,"/diff_skipper/reproducible_enriched_windows/", prefix, ".reproducible_enriched_windows.tsv"))
}

writeLines("OK", file.path(output_dir, "diffskip.done"))