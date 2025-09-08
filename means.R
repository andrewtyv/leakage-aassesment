set.seed(7)

# ---------- helpers ----------
read_traces_matrix <- function(traces_path) {
  files <- list.files(traces_path, "^trace_\\d+\\.txt$", full.names = TRUE)
  stopifnot(length(files) > 0)
  idx <- as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(files)))
  ord <- order(idx); files <- files[ord]
  X <- do.call(rbind, lapply(files, scan, quiet = TRUE))
  list(X = X, files = files, idx = idx[ord])
}

center_mean0 <- function(X, qs=1, qe=ncol(X)) {
  Xw <- X[, qs:qe, drop = FALSE]
  mu <- colMeans(Xw)
  sweep(Xw, 2, mu, "-")
}

# --- MBB,  ---
mbb_generate_from_rows_cols <- function(X0, row_pool, col_pool, n_out, B=64) {
  # X0: N×S zero-centered
  # row_pool: indx of rows
  # col_pool: ondx of cols
  N <- nrow(X0); S <- ncol(X0); K <- ceiling(S / B)
  syn <- matrix(0, nrow = n_out, ncol = S)
  
    #valid starts
  max_start <- S - B + 1
  col_pool  <- col_pool[col_pool >= 1 & col_pool <= S]
  valid_starts <- col_pool[col_pool <= max_start]
  if (length(valid_starts) == 0) valid_starts <- 1:max_start
  
  for (r in seq_len(n_out)) {
    acc <- numeric(0)
    for (k in 1:K) {
      i  <- sample(row_pool, 1)
      s0 <- sample(valid_starts, 1)
      acc <- c(acc, X0[i, s0:(s0+B-1)])
    }
    syn[r, ] <- acc[1:S]
  }
  syn
}

# optional: non-linear shape tweak (keeps mean close to 0)
shape_tweak <- function(X, delta = 0.15) {
  # odd, mean-preserving-ish transform: x -> x + delta*sign(x)*(abs(x)-median(|x|))
  med_abs <- apply(abs(X), 2, median)
  X + sweep(sign(X) * (abs(X) - rep(med_abs, each = nrow(X))), 2, delta, `*`)
}

# ---------- TVLA (Welch) & KSLA per sample ----------
tvla_welch <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))
  sapply(seq_len(ncol(A)), function(i) t.test(A[,i], B[,i], var.equal = FALSE)$statistic)
}
ksla_stat <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))
  sapply(seq_len(ncol(A)), function(i) suppressWarnings(ks.test(A[,i], B[,i])$statistic))
}

# ---------- MAIN: proper TVLA with fixed vs random inside ONE dataset (row-pooled MBB) ----------
run_synth_tvla_from_real <- function(
    traces_path,
    inputs_file,          
    v           = 0.5,   
    n_fixed     = 7000,
    n_random    = 1000,
    qs          = 1,
    qe          = NULL,
    B           = 64,
    win_fixed   = NULL,   # 
    win_random  = NULL,   # 
    random_transform = FALSE, #
    out_pdf     = "tvla_ksla_synth_fixed_vs_random.pdf",
    thr_t       = 4.5,
    thr_ks_q    = 0.99     # 
) {
  # 1) read & window
  rd <- read_traces_matrix(traces_path)
  X  <- rd$X
  idx_files <- rd$idx
  if (is.null(qe)) qe <- ncol(X)
  X0 <- center_mean0(X, qs, qe)      # mean = 0, keep shape/correlation
  S  <- ncol(X0)
  
  # 1b) align inputs to traces by numeric index (inputs 0-based → +1)
  inputs_mat <- as.matrix(read.table(inputs_file, header = FALSE, sep = " ",
                                     col.names = paste0("V", 1:7)))
  inputs <- inputs_mat[idx_files + 1, , drop = FALSE]
  
  fixed_rows  <- which(inputs[,1] == v)
  random_rows <- which(inputs[,1] != v)
  if (length(fixed_rows) < 2 || length(random_rows) < 2)
    stop("Not enough real traces in fixed or random group to build row-pooled MBB.")
  
  # 2) define COLUMN pools
  if (is.null(win_fixed))  win_fixed  <- 1:S
  if (is.null(win_random)) win_random <- 1:S
  
  # 3) synthetic fixed & random via row-conditioned MBB
  fixed_mat  <- mbb_generate_from_rows_cols(X0, fixed_rows,  win_fixed,  n_out = n_fixed,  B = B)
  random_mat <- mbb_generate_from_rows_cols(X0, random_rows, win_random, n_out = n_random, B = B)
  
  if (random_transform) {
    random_mat <- shape_tweak(random_mat, delta = 0.15)  # small shape change, mean≈0
    # re-center tiny drift if any:
    random_mat <- sweep(random_mat, 2, colMeans(random_mat), "-")
  }
  
  # 4) (labels kept for completeness; not used in tests below)
  inputs_V1 <- c(rep(0.5, n_fixed), runif(n_random, -1, 1))
  
  # 5) TVLA & KSLA
  t_vals  <- tvla_welch(fixed_mat, random_mat)
  ks_vals <- ksla_stat(fixed_mat, random_mat)
  
  ks_thr <- quantile(ks_vals, thr_ks_q)
  
  # 6) plots
  pdf(out_pdf, width = 10, height = 6)
  par(mfrow = c(2,1))
  plot(t_vals, type="l", col="steelblue", lwd=1.2,
       main = sprintf("TVLA (Welch) — fixed vs random   Nf=%d, Nr=%d", n_fixed, n_random),
       xlab = "Sample index", ylab = "t-value")
  abline(h=c(-thr_t, thr_t), col="red", lty=2)
  
  plot(ks_vals, type="l", col="darkgreen", lwd=1.2,
       main = sprintf("KSLA — fixed vs random   KS threshold = q%.2f", thr_ks_q),
       xlab = "Sample index", ylab = "KS statistic")
  abline(h=ks_thr, col="red", lty=2)
  dev.off()
  
  cat("TVLA: |t| >", thr_t, "at", sum(abs(t_vals) > thr_t), "samples\n")
  cat("KSLA: KS >", round(ks_thr, 4), "at", sum(ks_vals > ks_thr), "samples\n")
  
  invisible(list(
    t_values  = t_vals,
    ks_values = ks_vals,
    ks_thr    = ks_thr,
    fixed_n   = n_fixed,
    random_n  = n_random,
    fixed_rows_count  = length(fixed_rows),
    random_rows_count = length(random_rows),
    inputs_V1 = inputs_V1
  ))
}

# ---------- EXAMPLES ----------
path_unprotected <- "/Users/andrew/Desktop/thesis/only-traces/capture_traces/unprotected"
inputs_unprotected <- "/Users/andrew/Desktop/thesis/only-traces/capture_traces/unprotected/inputs.txt"

qs <- 1; qe <- 24430
poi    <- 534:24430
nonpoi <- 1:533

# no  shape tweak
resA <- run_synth_tvla_from_real(
  traces_path   = path_unprotected,
  inputs_file   = inputs_unprotected,
  v             = 0.5,
  n_fixed       = 1715,
  n_random      = 285,
  qs = qs, qe = qe,
  B = 64,
  win_fixed  = nonpoi,         # COLUMN pool for fixed
  win_random = 1:24430,        # COLUMN pool for random
  random_transform = FALSE,
  out_pdf = "tvla_ksla_fixed_random_from_real_POI_mixed.pdf"
)

#shape tweak 
resB <- run_synth_tvla_from_real(
  traces_path   = path_unprotected,
  inputs_file   = inputs_unprotected,
  v             = 0.5,
  n_fixed       = 1715,
  n_random      = 285,
  qs = qs, qe = qe,
  B = 64,
  win_fixed  = nonpoi,
  win_random = 1:24430,
  random_transform = TRUE,
  out_pdf = "tvla_ksla_fixed_random_from_real_POI_mixed_tweak.pdf"
)
