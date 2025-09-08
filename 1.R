# TVLA
# -----------------------------------------------------------------------------
# tvla_from_inputs()
# Runs a windowed, non-specific t-test (TVLA) over a selected sample window
# [qs:qe] for two trace groups:
#   - "fixed": traces where the first input value equals v
#   - "random": traces where the first input value differs from v
#
# It:
#   1) Loads inputs and aligns them with trace files by numeric index.
#   2) Splits traces into fixed vs random groups based on inputs[,1] == v.
#   3) Computes a Welch's t-value per sample in the window (column-wise).
#   4) Saves a PDF plot of t-values with ±threshold overlays and a CSV of values.
#   5) Builds a "power curve" by recomputing max |t| as traces-per-group increases.
#      (Helps estimate how many traces are needed to exceed the threshold.)
#
# Arguments:
#   name        : Label used in output filenames and plot titles.
#   traces_path : Directory with trace_*.txt files; each contains one trace vector.
#   inputs_file : Path to 'inputs.txt' with 7 columns (V1..V7), one row per trace.
#   out_dir     : Output directory for PDFs/CSVs.
#   v           : The constant value used to define the "fixed" group (default 0.5).
#   threshold   : TVLA pass/fail line; |t| above this indicates leakage (default 4.5).
#   min_fixed   : Minimum traces per group to start the power-curve analysis.
#   qs, qe      : 1-based inclusive sample window; pick 1..24430 to cover full trace.
#
# Returns (invisible):
#   A list with t-values, leakage points, counts, window, and output file paths,
#   plus power-curve (m values and max |t| per m).
# -----------------------------------------------------------------------------


tvla_from_inputs <- function(
    name,
    traces_path,
    inputs_file,
    out_dir     = "/Users/andrew/Desktop/protectedvsunprotected/",
    v           = 0.5, # const value of neuron 
    threshold   = 4.5, 
    min_fixed   = 10, 
    qs, # qs and qe is for windowed TVLA, basicaly u can go from 1 to 24430 to work with all file 
    qe
) {
  cat("Running windowed TVLA for:", name, "\n\n")
  
  
  #ensure output directory exists 
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    cat("Created output directory:", out_dir, "\n\n")
  }
  
  #load inputs as matrix with columns V1..V7 
  inputs_mat <- as.matrix(
    read.table(inputs_file, header = FALSE, sep = " ",
               col.names = paste0("V", 1:7))
  )
  N <- nrow(inputs_mat) # nomber of traces expected from inputs 
  
  #collect and sort all trace files trace_*.txt by indx
  trace_files <- list.files(traces_path, "^trace_\\d+\\.txt$",
                            full.names = TRUE)
  extract_index <- function(f)
    as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(f)))
  idx <- sapply(trace_files, extract_index)
  ord <- order(idx)
  trace_files <- trace_files[ord]
  idx <- idx[ord]
  
  #read all traces to matrix rows = traces, col = samples
  traces <- do.call(rbind, lapply(trace_files, scan, quiet = TRUE))
  
  # Align inputs to traces by index; inputs file is 0-based, so +1 
  inputs <- inputs_mat[idx + 1, , drop = FALSE]
  
  #split into fixed vs random 
  fixed_idx  <- which(inputs[,1] == v)
  random_idx <- which(inputs[,1] != v)
 
  #extract the sample window [qs:qe] for both groups 
  fixed_win  <- traces[fixed_idx, qs:qe, drop = FALSE]
  random_win <- traces[random_idx, qs:qe, drop = FALSE]
  
  #compute Welch's t-value per sample
  n_win    <- ncol(fixed_win)
  t_values <- sapply(seq_len(n_win), function(i) {
    t.test(fixed_win[,i], random_win[,i], var.equal = FALSE)$statistic
  })
  
  #save TVLA plot
  pdf_file <- file.path(out_dir, paste0("tvla2-neuron_", name, ".pdf"))
  pdf(pdf_file, width = 10, height = 5)
  plot(t_values, type="l", lwd=1.5,
       main = paste(" TVLA -", name),
       xlab = sprintf("Sample index (%d-%d)", qs, qe),
       ylab = "t-value")
  abline(h = c(-threshold, threshold), col = "red", lty = 2)
  grid()
  dev.off()
  cat("Saved plot to:", pdf_file, "\n")
  
  #save raw t-values to CSV
  csv_file <- file.path(out_dir, paste0("tvalues2-neuron_", name, ".csv"))
  write.csv(
    data.frame(
      sample  = seq_len(n_win) + qs - 1,
      t_value = t_values
    ),
    file      = csv_file,
    row.names = FALSE
  )
  cat("Saved data to:", csv_file, "\n")
  
  #Identify leakage points where |t| exceeds threshold
  leaks <- which(abs(t_values) > threshold)
  cat(name, "-leakage at", length(leaks),
      "samples in [", qs, "-", qe, "]\n\n")
  
  #power curve 
  n_fixed  <- nrow(fixed_win)
  n_random <- nrow(random_win)
  m_max    <- min(n_fixed, n_random)
  
  if (m_max < min_fixed) stop("Not enough traces per group for the power curve.")
  
  # Choose up to 100 steps from min_fixed..m_max (unique rounded values)  
  steps <- min(100L, m_max - min_fixed + 1L)
  m_seq <- unique(round(seq(min_fixed, m_max, length.out = steps)))
  
  # For each m, compute max |t| over the window using the first m traces per group
  max_tvals <- numeric(length(m_seq))
  
  for (j in seq_along(m_seq)) {
    m <- m_seq[j]
    fw <- fixed_win[1:m, , drop = FALSE]
    rw <- random_win[1:m, , drop = FALSE]
    t_vals_m <- sapply(seq_len(n_win), function(i) {
      t.test(fw[, i], rw[, i], var.equal = FALSE)$statistic
    })
    max_tvals[j] <- max(abs(t_vals_m), na.rm = TRUE)   # <<--- max |t|
  }
  
  # Save power-curve as a simple barplot (max |t| vs traces-per-group)
  
  power_pdf <- file.path(out_dir, paste0("tvla_power_curve_", name, ".pdf"))
  pdf(power_pdf, width = 9, height = 5)
  barplot(height = max_tvals, names.arg = m_seq,
          xlab = "Number of traces per group (fixed & random)",
          ylab = "Max |t| value in window",
          main = paste("TVLA max-|t| histogram –", name),
          col = "skyblue")
  grid(nx = NA, ny = NULL)
  dev.off()
  cat("Saved power-curve histogram to:", power_pdf, "\n")
  
  #save CSV
    power_csv <- file.path(out_dir, paste0("tvla_power_curve_", name, ".csv"))
  write.csv(
    data.frame(
      traces_per_group = m_seq,
      max_abs_t        = max_tvals
    ),
    file = power_csv, row.names = FALSE
  )
  cat("Saved power-curve data to:", power_csv, "\n\n")
  
  invisible(list(
    t_values            = t_values,
    leakage_points      = leaks,
    fixed_count         = length(fixed_idx),
    random_count        = length(random_idx),
    window              = c(qs, qe),
    pdf_tvalues         = pdf_file,
    csv_tvalues         = csv_file,
    power_curve_m       = m_seq,            
    pdf_power_curve     = power_pdf,        
    csv_power_curve     = power_csv         
  ))
}
result_unprot <- tvla_from_inputs(
  name        = "unprotected_new_nn",
  traces_path = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected",
  inputs_file = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected/inputs.txt",
  qs = 1,
  qe = 24430
)

result_prot <- tvla_from_inputs(
  name        = "protected_new_nn",
  traces_path = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/protected0.5",
  inputs_file = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/protected0.5/inputs.txt",
  qs = 1,
  qe = 24430
)