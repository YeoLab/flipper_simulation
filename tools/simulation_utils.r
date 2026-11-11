######################### Simulating the initial Gene-counts #########################

# Simulate the initial reads using a method identical to the one used by deseq. 
simulate_deseq_mu_and_counts = function(n_genes = 1000, n_samples = 8, betaShift = 0, betaSD = 0,
                                         interceptMean = 10, interceptSD = 3, dispMeanRel = function(x) 2/x + 0.1,
                                         sizeFactors = rep(1, n_samples)) {
    
    # Simulate gene-wise coefficients (counts and treatment effect beta)
    beta = cbind(rnorm(n_genes, interceptMean, interceptSD), rnorm(n_genes, betaShift, betaSD))
    
    # Convert intercepts to base mean expression on count scale.
    base_mean = 2^(beta[, 1])
    
    # Compute gene-wise dispersion using the mean–dispersion relationship.
    dispersion = dispMeanRel(base_mean)

    # Create sample condition labels.
    condition = factor(rep(c("A", "B"), times = c(ceiling(n_samples/2), floor(n_samples/2))))
    
    # Build design matrix for condition effect.
    x = if (n_samples > 1) stats::model.matrix(~ condition) else cbind(rep(1, n_samples), rep(0, n_samples))

    # Compute expected mean counts (mu) for each gene and sample.
    mu = t(2^(x %*% t(beta)) * sizeFactors)

    # Simulate negative binomial counts using gene-wise mu and dispersion.
    countData = matrix(rnbinom(n_samples * n_genes, mu = mu, size = 1 / dispersion), ncol = n_samples)
    
    # Store counts as integers to mimic real RNA-seq data.
    storage.mode(countData) = "integer"

    # Assign sample names.
    colnames(countData) = paste0("sample", seq_len(n_samples))

    # Return simulated counts, true mu, and ground-truth parameters.
    list(countData = data.frame(countData, check.names = FALSE), mu = mu, 
         trueVals = data.frame(trueIntercept = beta[, 1], trueBeta = beta[, 2], trueDisp = dispersion))
}

# Couple gene lengths to gene counts (preserves marginals, adds mild correlation).
couple_lengths_to_counts = function(X, L, strength = 0.3) {
    
    # Get number of genes.
    n = length(X)
    
    # Obtain ordering indices for X and L.
    ix = order(X)
    il = order(L)
    
    # Create ideal ranks.
    ranks = seq_len(n)
    
    # Add Gaussian noise to ranks to control strength of association.
    noise = rnorm(n, mean = 0, sd = (1 - strength) * n / 5)
    
    # Apply noise and clamp ranks to valid range.
    noisy_ranks = ranks + noise
    noisy_ranks = pmax(1, pmin(n, round(noisy_ranks)))
    
    # Sort L and remap values according to noisy ranks.
    L_sorted = L[il]
    L_mapped = L_sorted[order(noisy_ranks)]
    
    # Reassign mapped L values to match ordering of X.
    L_out = integer(n)
    L_out[ix] = L_mapped
    
    # Return output vector with controlled correlation to X.
    L_out
}

attach_abundance_and_length = function(mu_mat, gene_length_bp) {
    
    # Get number of genes.
    n_genes = nrow(mu_mat)

    # Calculate reference length. 
    L_ref = stats::median(gene_length_bp)

    # Compute baseline mean expression per gene across samples.
    mu_gene_mean = rowMeans(mu_mat)

    # Define abundance such that abundance. 
    abundance = mu_gene_mean * (L_ref / gene_length_bp)

    # Return gene length metadata and abundance quantities.
    list(gene_length_bp = as.integer(gene_length_bp), L_ref = as.numeric(L_ref),
         abundance = as.numeric(abundance), mu_gene_mean = as.numeric(mu_gene_mean))
}

simulate_gene_counts_with_length = function(n_genes, n_samples, beta_shift,
                                            beta_sd, Mean_ex, SD_ex, dispMeanRel,
                                            lens_bp, sizeFactors = rep(1, n_samples)) {
    
    # Generate DESeq-style gene counts.
    deseq = simulate_deseq_mu_and_counts(n_genes = n_genes, n_samples = n_samples, betaShift = beta_shift,
                                        betaSD = beta_sd, interceptMean = Mean_ex, interceptSD = SD_ex,
                                        dispMeanRel = dispMeanRel,sizeFactors = sizeFactors)
    
    # Sample gene lengths directly from empirical distribution.
    gene_length_bp = as.integer(sample(lens_bp, size = n_genes, replace = TRUE))
    
    # Use baseline mean per gene as the coupling signal.
    X_signal = rowMeans(deseq$mu)
    gene_length_bp = couple_lengths_to_counts(X = X_signal, L = gene_length_bp, strength = 0.3)
    
    # Attach abundance so mu is preserved exactly when combined with length later.
    latent = attach_abundance_and_length(deseq$mu, gene_length_bp)
    
    # Add gene_name and return everything.
    gene_counts = deseq$countData
    gene_counts$gene_name = paste0("gene_", seq_len(nrow(gene_counts)))

    # Return results. 
    list(gene_counts = gene_counts, mu = deseq$mu, trueVals = deseq$trueVals,
         gene_length_bp = latent$gene_length_bp, abundance = latent$abundance,
         L_ref = latent$L_ref, mu_gene_mean = latent$mu_gene_mean
        )
}

downsample_counts= function(counts, mean_prop = 0.1) {
  
    # Identify numeric count columns (e.g., exclude gene_name).
    is_count_col = vapply(counts, is.numeric, logical(1))
    count_df = counts[, is_count_col, drop = FALSE]

    # Number of genes.
    n_genes = nrow(count_df)

    # Draw one proportion per gene.
    p_i = rep(mean_prop, n_genes)

    # Organize counts into a matrix.
    mat = as.matrix(count_df)

    # Perform the downsampling. 
    down = t(mapply(
        FUN = function(row_counts, p) {
        rbinom(length(row_counts), size = row_counts, prob = p)
        },
        split(mat, row(mat)),
        p_i
    ))

    # Restore dims & names.
    dim(down) = dim(mat)
    rownames(down) = rownames(mat)
    colnames(down) = colnames(mat)
    
    # Convert to dataframe. 
    down_df = as.data.frame(down)

    # Put non-count columns back. 
    if (any(!is_count_col)) {
    out = counts
    out[, is_count_col] = down_df
    return(out)
    }

    return(down_df)
}

######################### Splitting counts into windows #########################

# Partition windows into features. 
partition_into_features = function(total_windows, n_features) {
    # Handle cases of 0 or 1. 
    if (total_windows <= 0 || n_features <= 0) return(integer(0))
    if (n_features == 1) return(as.integer(total_windows))

    # Give each feature at least one window, then distribute the remainder.
    remainder = total_windows - n_features
    probs = rgamma(n_features, shape = 1, rate = 1)
    probs = probs / sum(probs)

    extra = as.vector(rmultinom(1, size = remainder, prob = probs))
    return(as.integer(1 + extra))
}
# This functions purpose is to randomly partition the windows into features after the size of the 
# features has been decided. 

# Build contiguous exon/intron assignments for one gene.
make_feature_layout = function(n_win, active_windows, mean_exon_windows, mean_intron_windows) {
    
    # calcualte inactive windows as a function of active windows. 
    inactive_windows = n_win - active_windows

    # Handle the all-exon case.
    if (inactive_windows <= 0) {
        return(data.frame(feature_type = rep("exon", n_win),
                          feature_id = rep("exon_1", n_win),
                          stringsAsFactors = FALSE))
    }

    # Choose the number of exon and intron features as a function of the mean number of windows per feature. 
    n_exon_features = max(1, as.integer(round(active_windows / mean_exon_windows)))
    n_intron_features = max(1, as.integer(round(inactive_windows / mean_intron_windows)))

    # Make the features alternate between introns and exons.
    while (abs(n_exon_features - n_intron_features) > 1) {
        if (n_exon_features > n_intron_features) {
            n_exon_features = n_exon_features - 1
        } else {
            n_intron_features = n_intron_features - 1
        }
    }

    # Partition the windows into each feature (will end up with exon lengths around the mean values). 
    exon_lengths = partition_into_features(active_windows, n_exon_features)
    intron_lengths = partition_into_features(inactive_windows, n_intron_features)

    # Initialize variables for while loop. 
    feature_type = character(0)
    feature_id = character(0)
    exon_counter = 0
    intron_counter = 0
    e = 1
    j = 1

    # Determine whether to start with an exon or intron (must line up correctly for alternating introns/exons). 
    start_with_exon = n_exon_features >= n_intron_features
    current_type = if (start_with_exon) "exon" else "intron"

    # Create a while loop that assigns feature types and names to each window.  
    while (e <= length(exon_lengths) || j <= length(intron_lengths)) {
        if (current_type == "exon") {
            if (e <= length(exon_lengths)) {
                exon_counter = exon_counter + 1
                feature_type = c(feature_type, rep("exon", exon_lengths[e]))
                feature_id = c(feature_id, rep(paste0("exon_", exon_counter), exon_lengths[e]))
                e = e + 1
            }
            current_type = "intron"
        } else {
            if (j <= length(intron_lengths)) {
                intron_counter = intron_counter + 1
                feature_type = c(feature_type, rep("intron", intron_lengths[j]))
                feature_id = c(feature_id, rep(paste0("intron_", intron_counter), intron_lengths[j]))
                j = j + 1
            }
            current_type = "exon"
        }
    }

    return(data.frame(feature_type = feature_type, feature_id = feature_id, stringsAsFactors = FALSE))
}

split_counts_into_windows = function(df, gene_length_bp = NULL, window_bp = 100, length_col = "gene_length_bp",
                                      exon_frac = 0.10, exon_peakiness = 0.2, intron_exon_balance = 0.02,
                                      intron_peakiness = 5, mean_exon_windows = 2, mean_intron_windows = 8,
                                      window_noise_sd = 0.6) {

    # Get the number of windows per gene based on gene length.
    n_windows_vec = pmax(1, ceiling(gene_length_bp / window_bp))
    gene_names = df$gene_name

    # Extract out just the count columns. 
    is_count_col = vapply(df, is.numeric, logical(1))
    is_count_col[names(is_count_col) == length_col] = FALSE
    count_df = df[, is_count_col, drop = FALSE]

    # Get the sample names.
    sample_cols = colnames(count_df)

    # Convert to a matrix.
    count_mat = as.matrix(count_df)

    # Get the total number of genes and samples.
    n_genes = nrow(df)
    n_samp = ncol(count_mat)

    # Create a list to store results.
    expanded_list = vector("list", length = n_genes)

    # Iterate over all genes. 
    for (i in seq_len(n_genes)) {
        n_win = n_windows_vec[i]

        # Get the fraction of windows that are active while ensuring at least one active window.
        exon_frac_i = max(1 / n_win, min(exon_frac, 1.0))

        # Determine how many windows are active.
        k = as.integer(max(1, round(n_win * exon_frac_i)))

        # Create contiguous exon/intron structure.
        feature_layout = make_feature_layout(n_win = n_win, active_windows = k,
                                             mean_exon_windows = mean_exon_windows,
                                             mean_intron_windows = mean_intron_windows)

        # Combine the feature layouts with the gene names for proper feature IDs. 
        feature_layout$feature_id = paste0(gene_names[i], "_", feature_layout$feature_id)

        # Feature types are used to split the gene into broad exonic and intronic regions. 
        active_idx = which(feature_layout$feature_type == "exon")
        inactive_idx = which(feature_layout$feature_type == "intron")

        # Allocate all reads to exons when the gene has no intron windows.
        intron_prob = if (length(inactive_idx) > 0) intron_exon_balance else 0

        # Generate a gene-level baseline distribution over active windows.
        alpha_gene = rep(exon_peakiness, length(active_idx))
        p0_active = rgamma(length(active_idx), shape = alpha_gene, rate = 1)
        p0_active = p0_active / sum(p0_active)

        # Generate a baseline distribution over inactive windows if present.
        if (length(inactive_idx) > 0) {
            alpha_inactive = rep(intron_peakiness, length(inactive_idx))
            p0_inactive = rgamma(length(inactive_idx), shape = alpha_inactive, rate = 1)
            p0_inactive = p0_inactive / sum(p0_inactive)
        } else {
            p0_inactive = numeric(0)
        }

        # Store per-sample splits.
        splits = vector("list", length = n_samp)

        # Iterate over samples.
        for (s in seq_len(n_samp)) {
            X = as.integer(count_mat[i, s])

            # Handle zero or missing counts.
            if (is.na(X) || X <= 0) {
                full_counts = integer(n_win)
            } else {
                
                # Partition counts into inactive/background and active/enriched components.
                X_inactive = rbinom(1, size = X, prob = intron_prob)
                X_active = X - X_inactive

                full_counts = integer(n_win)

                # Add sample-specific noise to the active windows.
                p_active_s = p0_active * exp(rnorm(length(p0_active), mean = 0, sd = window_noise_sd))
                p_active_s = p_active_s / sum(p_active_s)

                # Add sample-specific noise to the inactive windows if intronic windows exist.
                if (length(p0_inactive) > 0) {
                    p_inactive_s = p0_inactive * exp(rnorm(length(p0_inactive), mean = 0, sd = window_noise_sd))
                    p_inactive_s = p_inactive_s / sum(p_inactive_s)
                } else {
                    p_inactive_s = numeric(0)
                }

                # Allocate counts across intronic windows.
                if (X_inactive > 0) {
                    inactive_counts = as.vector(rmultinom(1, size = X_inactive, prob = p_inactive_s))
                    full_counts[inactive_idx] = inactive_counts
                }

                # Allocate counts across exonic windows. 
                if (X_active > 0) {
                    active_counts = as.vector(rmultinom(1, size = X_active, prob = p_active_s))
                    full_counts[active_idx] = active_counts
                }
            }

            # Store sample-specific window counts.
            splits[[s]] = matrix(full_counts, ncol = 1)
        }

        # Combine sample splits into a matrix and assign column names.
        split_mat = do.call(cbind, splits)
        colnames(split_mat) = sample_cols

        # Store expanded gene-level data frame.
        expanded_list[[i]] = data.frame(gene_name = gene_names[i], feature_type = feature_layout$feature_type, feature_id = feature_layout$feature_id,
                                        split_mat, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    }

    # Combine all genes into the final result.
    result = do.call(rbind, expanded_list)
    rownames(result) = NULL

    return(result)
}

######################### For the simulation of binding and pulldown #########################

IPDataSpecZ = function(feature_counts, n_bind, n_samples,
                        zetaMean = 1.2, zetaSD = 0.01,
                        interceptMode = 0.2, interceptShape = 5,
                        interceptMin = 0.1, pulldown_a = 2, 
                        pulldown_b = 0.02, diff_frac = 0.5,
                        ip_size_factors = rep(1, n_samples/2), zeta_rep_sd = 0.05) {

    # Find which cols are numeric.
    num_cols = which(sapply(feature_counts, is.numeric))
    
    # Then pick the 2nd, 4th, 6th, 8th *numeric* columns.
    cols_to_check = num_cols[c(2, 4, 6, 8)]
    
    # Ensures that in at least 1 treatment group there are over 30 reads. 
    present_features = which(rowMeans(feature_counts[, cols_to_check] >= 10) == 1)
    
    # Sample potential bound genes only from expressed genes.
    to_bind = unname(sample(present_features, n_bind))
    
    # Enforce mean(intercepts) ≈ interceptMean. 
    rate = interceptShape / (interceptMode - interceptMin)
    intercepts = interceptMin + rgamma(n_bind, interceptShape, rate)
    
    # Identify log fold changes based on zetaMean. 
    log_fc = rnorm(n_bind, mean = zetaMean, sd = zetaSD)
    treat_effects = exp(log_fc)   # multiplicative, >0, centered at 1 if mean=0
    
    # Make it so that only a fraction of binding sites are affected by the treatment.
    n_diff = round(n_bind * diff_frac)
    is_diff = rep(FALSE, n_bind)
    is_diff[sample(seq_len(n_bind), n_diff)] = TRUE
    treat_effects[!is_diff] = 1
    
    # Combine into zeta matrix.
    zeta = cbind(intercept = intercepts, treatment = treat_effects)
    
    # Condition vector.
    colData = DataFrame(condition = factor(rep(c("A", "B"),
                                             times = c(ceiling(n_samples / 4), floor(n_samples / 4)))))
    x = stats::model.matrix.default(~ colData$condition)
    
    intercepts = zeta[, 1]
    multipliers = zeta[, 2]
    treat = x[, 2] 

    mu_base = matrix(intercepts, nrow = length(intercepts), ncol = length(treat), byrow = FALSE)
    
    mult_mat = matrix(1, nrow = n_bind, ncol = length(treat))
    mult_mat[, treat == 1] = multipliers

    # Add small replicate-specific variability to treatment effects.
    # This makes treatment IP replicates wobble around the same true site effect.
    if (zeta_rep_sd > 0) {
        trt_cols = which(treat == 1)
        zeta_rep_noise = matrix(
            exp(rnorm(n_bind * length(trt_cols), mean = 0, sd = zeta_rep_sd)),
            nrow = n_bind,
            ncol = length(trt_cols)
        )
        mult_mat[, trt_cols] = mult_mat[, trt_cols] * zeta_rep_noise
    }
    
    mu = mu_base * mult_mat
    mu[mu < 0] = 0
    
    # Pull out IP columns from gene_counts. 
    ip_cols = seq(2, n_samples, by = 2)
    feature_ip = data.matrix(feature_counts[, ip_cols])
    
    # Mean binding hits (n_bind x n_ip_samples)
    mean_binding_hits = feature_ip[to_bind, , drop = FALSE] * mu

    # Calculate pulldown noise as a function of mean binding hits and 
    pulldown_sd_mat = sqrt(pulldown_a / pmax(mean_binding_hits, 1) + pulldown_b)
    
    pulldown_noise = matrix(rnorm(n_bind * ncol(mean_binding_hits), mean = 0, sd = as.vector(pulldown_sd_mat)),
                            nrow = n_bind, ncol = ncol(mean_binding_hits))
    
    pulldown_noise = exp(pulldown_noise)
    mean_binding_hits = mean_binding_hits * pulldown_noise

    # Add baseline binding offset while preserving treatment effects.
    offset_mat = matrix(10, nrow = n_bind, ncol = n_samples/2)
    
    # Apply treatment scaling to treatment samples only.
    offset_mat[, treat == 1] = 10 * mult_mat[, treat == 1]
    
    mean_binding_hits = mean_binding_hits + offset_mat
    
    # Apply IP pulldown size factors (column-wise).
    mean_binding_hits = sweep(mean_binding_hits, 2, as.numeric(ip_size_factors), `*`)
    
    # Generate IP data with scaled means.
    IPData = matrix(as.integer(round(as.vector(mean_binding_hits))), ncol = n_samples/2, byrow = FALSE)
    
    trueVals = data.frame(trueIntercept = zeta[, 1], trueZeta = zeta[, 2])
    
    return(list(IPData = IPData, trueVals = trueVals, to_bind = feature_counts$feature_id[c(to_bind)]))
}