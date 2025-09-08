set.seed(7)

# ---------- helpers ----------
read_traces_matrix <- function(traces_path) {
  files <- list.files(traces_path, "^trace_\\d+\\.txt$", full.names = TRUE)
  stopifnot(length(files) > 0)
  idx <- as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(files)))
  ord <- order(idx); files <- files[ord]
  X <- do.call(rbind, lapply(files, scan, quiet = TRUE))
  list(X = X, files = files)
}

center_mean0 <- function(X, qs=1, qe=ncol(X)) {
  Xw <- X[, qs:qe, drop = FALSE]
  mu <- colMeans(Xw)
  sweep(Xw, 2, mu, "-")
}

mbb_generate_from_pool <- function(X0, pool_cols, n_out, B=64) {
  # X0: N×S zero-centered; pool_cols: allowed start indices for blocks
  N <- nrow(X0); S <- ncol(X0); K <- ceiling(S / B)
  syn <- matrix(0, nrow = n_out, ncol = S)
  valid_starts <- pool_cols[pool_cols <= (max(pool_cols) - B + 1)]
  if (length(valid_starts) == 0) valid_starts <- 1:(S - B + 1)
  for (r in 1:n_out) {
    acc <- numeric(0)
    for (k in 1:K) {
      i  <- sample.int(N, 1)
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

# ---------- MAIN: proper TVLA with fixed vs random inside ONE dataset ----------
run_synth_tvla_from_real <- function(
    traces_path,
    n_fixed    = 7000,
    n_random   = 1000,
    qs         = 1,
    qe         = NULL,
    B          = 64,
    win_fixed  = NULL,   # integer vector of columns (e.g., non-POI)
    win_random = NULL,   # integer vector of columns (e.g., POI ∪ non-POI)
    random_transform = FALSE, # apply shape tweak to random after MBB
    out_pdf    = "tvla_ksla_synth_fixed_vs_random.pdf",
    thr_t      = 4.5,
    thr_ks_q   = 0.99     # KS threshold as quantile of KS statistics
) {
  # 1) read & window
  rd <- read_traces_matrix(traces_path)
  X  <- rd$X
  if (is.null(qe)) qe <- ncol(X)
  X0 <- center_mean0(X, qs, qe)      # mean≈0, keep shape/correlation
  S  <- ncol(X0)
  
  # 2) define pools
  if (is.null(win_fixed))  win_fixed  <- 1:S
  if (is.null(win_random)) win_random <- 1:S
  
  # 3) synthetic fixed & random via MBB
  fixed_mat  <- mbb_generate_from_pool(X0, win_fixed,  n_out = n_fixed,  B = B)
  random_mat <- mbb_generate_from_pool(X0, win_random, n_out = n_random, B = B)
  
  if (random_transform) {
    random_mat <- shape_tweak(random_mat, delta = 0.15)  # small shape change, mean≈0
    # re-center small drift if any:
    mcol <- colMeans(random_mat)
    random_mat <- sweep(random_mat, 2, mcol, "-")
  }
  
  # 4) build "inputs": V1=0.5 for fixed, random U for random (classic TVLA labelling)
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
    t_values = t_vals,
    ks_values = ks_vals,
    ks_thr = ks_thr,
    fixed_n = n_fixed,
    random_n = n_random,
    inputs_V1 = inputs_V1
  ))
}
# 1) шляхи до РЕАЛЬНИХ трас (один сет, наприклад unprotected)
path_unprotected <- "/Users/andrew/Desktop/thesis/only-traces/capture_traces/unprotected"

# 2) задай вікно (якщо треба) і POI/non-POI для басейнів блоків
qs <- 1; qe <- 24430  # або твоє реальне вікно
poi    <- 534:24430
nonpoi <- 1:533

# 3) варіант A: fixed з non-POI, random з POI∪non-POI (більш «колюча» форма у random)
resA <- run_synth_tvla_from_real(
  traces_path   = path_unprotected,
  n_fixed       = 8600,
  n_random      = 1400,
  qs = qs, qe = qe,
  B = 64,
  win_fixed  = nonpoi,
  win_random = 1:24430,
  random_transform = FALSE,                    # поки без дод. нелінійності
  out_pdf = "tvla_ksla_fixed_random_from_real_POI_mixed.pdf"
)

# 4) варіант B: як вище + легка нелінійна «shape tweak» у random,
#    щоб mean залишився ≈0, але хвости/форма трохи змінились
resB <- run_synth_tvla_from_real(
  traces_path   = path_unprotected,
  n_fixed       = 8600,
  n_random      = 1400,
  qs = qs, qe = qe,
  B = 64,
  win_fixed  = nonpoi,
  win_random = 1:24430,
  random_transform = TRUE,                     # вмикаємо м’яку зміну форми
  out_pdf = "tvla_ksla_fixed_random_from_real_POI_mixed_tweak.pdf"
)
