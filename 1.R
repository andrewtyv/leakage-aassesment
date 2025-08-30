tvla_from_inputs <- function(
    name,
    traces_path,
    inputs_file,
    out_dir     = "/Users/andrew/Desktop/protectedvsunprotected/",
    v           = 0.5,
    threshold   = 4.5,
    min_fixed   = 10,
    qs,
    qe
) {
  cat("Running windowed TVLA for:", name, "\n\n")
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    cat("Created output directory:", out_dir, "\n\n")
  }
  
  inputs_mat <- as.matrix(
    read.table(inputs_file, header = FALSE, sep = " ",
               col.names = paste0("V", 1:7))
  )
  N <- nrow(inputs_mat)
  
  trace_files <- list.files(traces_path, "^trace_\\d+\\.txt$",
                            full.names = TRUE)
  extract_index <- function(f)
    as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(f)))
  idx <- sapply(trace_files, extract_index)
  ord <- order(idx)
  trace_files <- trace_files[ord]
  idx <- idx[ord]
  
  traces <- do.call(rbind, lapply(trace_files, scan, quiet = TRUE))
  
  inputs <- inputs_mat[idx + 1, , drop = FALSE]
  fixed_idx  <- which(inputs[,1] == v)
  random_idx <- which(inputs[,1] != v)
 
  fixed_win  <- traces[fixed_idx, qs:qe, drop = FALSE]
  random_win <- traces[random_idx, qs:qe, drop = FALSE]
  
  n_win    <- ncol(fixed_win)
  t_values <- sapply(seq_len(n_win), function(i) {
    t.test(fixed_win[,i], random_win[,i], var.equal = FALSE)$statistic
  })
  
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
  
  leaks <- which(abs(t_values) > threshold)
  cat(name, "-leakage at", length(leaks),
      "samples in [", qs, "-", qe, "]\n\n")
  
  
  n_fixed  <- nrow(fixed_win)
  n_random <- nrow(random_win)
  m_max    <- min(n_fixed, n_random)
  
  if (m_max < min_fixed) stop("Not enough traces per group for the power curve.")
  
  steps <- min(100L, m_max - min_fixed + 1L)
  m_seq <- unique(round(seq(min_fixed, m_max, length.out = steps)))
  
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
  
  # --- побудова гістограми ---
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
  
  # --- збереження CSV ---
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
    power_curve_m       = m_seq,            # NEW
    pdf_power_curve     = power_pdf,        # NEW
    csv_power_curve     = power_csv         # NEW
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