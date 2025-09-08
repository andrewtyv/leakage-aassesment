# Kolmogorov–Smirnov Leakage Assessment (KSLA)  
# “KSLA” ≈ TVLA but using the two-sample KS test instead of t-test

ksla_from_inputs <- function(
    name,
    traces_path,
    inputs_file,
    v         = 0.5,
    threshold = 0.2,    # example KS-statistic threshold
    min_fixed = 10
) {
  cat("Running KSLA for:", name, "\n\n")
  
  # 1) Read inputs matrix N×M (columns V1..V7). We'll group by V1 == v
  inputs_mat <- as.matrix(
    read.table(inputs_file, header=FALSE, sep="", col.names=paste0("V",1:7))
  )
  N <- nrow(inputs_mat)
  
  # 2) Collect & sort trace filenames by numeric index in trace_<idx>.txt
  trace_files <- list.files(traces_path, "^trace_\\d+\\.txt$", full.names=TRUE)
  extract_idx <- function(f) as.integer(sub("^trace_(\\d+)\\.txt$","\\1",basename(f)))
  idx <- sapply(trace_files, extract_idx)
  ord <- order(idx)
  trace_files <- trace_files[ord]
  idx         <- idx[ord]
  
  # 3) Read all traces into matrix N×S (rows=traces, cols=samples)
  trace_list <- lapply(trace_files, scan, quiet=TRUE)
  traces     <- do.call(rbind, trace_list)
  
  # 4) Split into fixed vs random by first input value (V1)
  first_vals  <- inputs_mat[idx+1,1]     # align inputs (file indices start at 0)
  fixed_idx   <- which(first_vals == v)
  random_idx  <- which(first_vals != v)
  if (length(fixed_idx) < min_fixed || length(random_idx) < 1)
    stop("Not enough traces in one of the groups.")
  
  fixed_mat   <- traces[fixed_idx, , drop=FALSE]
  random_mat  <- traces[random_idx,, drop=FALSE]
  
  # 5) Compute two-sample KS statistic per sample column (time index)
  S <- ncol(traces)
  ks_values <- numeric(S)
  for (t in seq_len(S)) {
    ks_values[t] <- ks.test(fixed_mat[,t], random_mat[,t])$statistic
  }
  plot_dir <- "/Users/andrew/Desktop/protectedvsunprotected" 
  
  # 6) Plot KS curve over time with a horizontal threshold
  pdf(file.path(plot_dir, paste0("ksla_", name, ".pdf")), width=10, height=5)
  plot(ks_values, type="l", lwd=1.5,
       main=paste("KSLA –", name),
       xlab="Sample index", ylab="KS statistic")
  abline(h=threshold, col="red", lty=2)
  grid()
  dev.off()
  
  # 7) Save KS values to CSV (sample index + KS statistic)
  write.csv(data.frame(sample=seq_len(S), ks=ks_values),
            file=paste0("ksla_values_", name, ".csv"),
            row.names=FALSE)
  
  # Report samples exceeding the KS threshold (potential leakage)
  leaks <- which(ks_values > threshold)
  cat(name, "- detected", length(leaks),
      "samples with KS >", threshold, "\n\n")
  
  invisible(list(
    ks_values    = ks_values,
    leakage_pts  = leaks,
    fixed_count  = length(fixed_idx),
    random_count = length(random_idx)
  ))
}

# --- Example runs on unprotected/protected datasets ---
res_unprot <- ksla_from_inputs(
  name        = "unprotected_new_nn",
  traces_path = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected",
  inputs_file = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected/inputs.txt",
  v           = 0.5,
  threshold   = 0.2
)

res_prot <- ksla_from_inputs(
  name        = "protected_new_nn",
  traces_path = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/protected0.5",
  inputs_file = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/protected0.5/inputs.txt", 
  v           = 0.5,
  threshold   = 0.2
)
