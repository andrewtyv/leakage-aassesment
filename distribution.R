if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")
if (!requireNamespace("e1071", quietly = TRUE)) install.packages("e1071")
library(nortest)
library(e1071)

analyze_window_to_pdf <- function(from_sample = 1000, to_sample = 2000, v = 0.5,
                                  traces_path, inputs_file, plot_dir) {
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  inputs_mat <- as.matrix(read.table(inputs_file, header = FALSE, sep = " ",
                                     col.names = paste0("V", 1:7)))
  
  trace_files <- list.files(traces_path, "^trace_\\d+\\.txt$", full.names = TRUE)
  extract_index <- function(f) as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(f)))
  idx <- sapply(trace_files, extract_index)
  ord <- order(idx)
  trace_files <- trace_files[ord]
  idx         <- idx[ord]
  
  read_partial_trace <- function(file, from = from_sample, to = to_sample) {
    trace <- scan(file, quiet = TRUE)
    if (length(trace) < to) return(rep(NA, to - from + 1))
    trace[from:to]
  }
  trace_list <- lapply(trace_files, read_partial_trace)
  traces     <- do.call(rbind, trace_list)
  
  first_vals <- inputs_mat[idx + 1, 1]
  fixed_idx  <- which(first_vals == v)
  random_idx <- which(first_vals != v)
  fixed_mat  <- traces[fixed_idx, , drop = FALSE]
  random_mat <- traces[random_idx, , drop = FALSE]
  fixed_all  <- as.vector(fixed_mat)
  random_all <- as.vector(random_mat)
  
  pdf(file.path(plot_dir, "qq_fixed.pdf"))
  qqnorm(fixed_all, main = "QQ-plot fixed")
  qqline(fixed_all, col = "red", lwd = 2)
  dev.off()
  
  pdf(file.path(plot_dir, "qq_random.pdf"))
  qqnorm(random_all, main = "QQ-plot random")
  qqline(random_all, col = "red", lwd = 2)
  dev.off()
  
  pdf(file.path(plot_dir, "hist_fixed.pdf"))
  hist(fixed_all, breaks = 50, freq = FALSE, main = "Hist fixed", xlab = "Amplitude")
  if (sd(fixed_all) > 0) curve(dnorm(x, mean(fixed_all), sd(fixed_all)),
                               add = TRUE, col = "blue", lwd = 2)
  dev.off()
  
  pdf(file.path(plot_dir, "hist_random.pdf"))
  hist(random_all, breaks = 50, freq = FALSE, main = "Hist random", xlab = "Amplitude")
  if (sd(random_all) > 0) curve(dnorm(x, mean(random_all), sd(random_all)),
                                add = TRUE, col = "blue", lwd = 2)
  dev.off()
  
  # --- Тести на нормальність для вікна ---
  cat("\n=== Normality tests for window", from_sample, "to", to_sample, "===\n")
  sw_sample_size <- 5000
  if (length(fixed_all) > sw_sample_size) {
    set.seed(1)
    sw_fixed  <- shapiro.test(sample(fixed_all, sw_sample_size))
    sw_random <- shapiro.test(sample(random_all, sw_sample_size))
  } else {
    sw_fixed  <- shapiro.test(fixed_all)
    sw_random <- shapiro.test(random_all)
  }
  ks_fixed  <- ks.test(fixed_all,  "pnorm", mean(fixed_all),  sd(fixed_all))
  ks_random <- ks.test(random_all, "pnorm", mean(random_all), sd(random_all))
  ad_fixed  <- nortest::ad.test(fixed_all)
  ad_random <- nortest::ad.test(random_all)
  
  cat("Fixed group:\n")
  cat("  Shapiro-Wilk p =", format.pval(sw_fixed$p.value), "\n")
  cat("  KS vs Normal  p =", format.pval(ks_fixed$p.value), "\n")
  cat("  AD test       p =", format.pval(ad_fixed$p.value), "\n\n")
  
  cat("Random group:\n")
  cat("  Shapiro-Wilk p =", format.pval(sw_random$p.value), "\n")
  cat("  KS vs Normal  p =", format.pval(ks_random$p.value), "\n")
  cat("  AD test       p =", format.pval(ad_random$p.value), "\n")
  
  skew_fixed  <- skewness(fixed_all)
  kurt_fixed  <- kurtosis(fixed_all)
  skew_random <- skewness(random_all)
  kurt_random <- kurtosis(random_all)
  
  cat("Skewness & Kurtosis:\n")
  cat("  Fixed:  skew =", skew_fixed, ", kurt =", kurt_fixed, "\n")
  cat("  Random: skew =", skew_random, ", kurt =", kurt_random, "\n\n")
}


# --- Друга функція: алгебраїчна перевірка нормальності на всіх даних ---
test_normality_full_data <- function(v = 0.5,
                                     traces_path,
                                     inputs_file) {
  cat("\n=== Starting full-data normality analysis ===\n")
  
  # 1. Зчитування inputs
  inputs_mat <- as.matrix(read.table(inputs_file, header = FALSE, sep = " ",
                                     col.names = paste0("V", 1:7)))
  
  # 2. Збір і сортування trace-файлів
  trace_files <- list.files(traces_path, "^trace_\\d+\\.txt$", full.names = TRUE)
  extract_index <- function(f) as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(f)))
  idx <- sapply(trace_files, extract_index)
  ord <- order(idx)
  trace_files <- trace_files[ord]
  idx         <- idx[ord]
  
  # 3. Зчитування всіх трас повністю
  cat("Reading all traces (this may take a while)...\n")
  trace_list <- lapply(trace_files, scan, quiet = TRUE)
  traces     <- do.call(rbind, trace_list)
  
  # 4. Matching: inputs ↔ traces
  first_vals <- inputs_mat[idx + 1, 1]
  fixed_idx  <- which(first_vals == v)
  random_idx <- which(first_vals != v)
  fixed_mat  <- traces[fixed_idx, , drop = FALSE]
  random_mat <- traces[random_idx, , drop = FALSE]
  fixed_all  <- as.vector(fixed_mat)
  random_all <- as.vector(random_mat)
  
  cat("Fixed samples:", length(fixed_all), " | Random samples:", length(random_all), "\n")
  
  # 5. Безпека: попередження про великі обсяги
  if (length(fixed_all) > 1e7 || length(random_all) > 1e7) {
    cat("⚠️ Warning: you're running full tests on >10 million points per group.\n",
        "This may be slow or memory-intensive.\n")
  }
  
  # 6. Алгебраїчні тести (усі — на full data)
  cat("\n=== Running tests (full data, no sampling) ===\n")
  
  # Kolmogorov–Smirnov
  ks_fixed  <- ks.test(fixed_all,  "pnorm", mean(fixed_all),  sd(fixed_all))
  ks_random <- ks.test(random_all, "pnorm", mean(random_all), sd(random_all))
  
  # Skewness / Kurtosis
  library(e1071)
  skew_fixed  <- skewness(fixed_all)
  kurt_fixed  <- kurtosis(fixed_all)
  skew_random <- skewness(random_all)
  kurt_random <- kurtosis(random_all)
  
  # Anderson–Darling — повна версія
  library(nortest)
  cat("Running Anderson–Darling test on full fixed data...\n")
  ad_fixed  <- nortest::ad.test(fixed_all)
  cat("Running Anderson–Darling test on full random data...\n")
  ad_random <- nortest::ad.test(random_all)
  
  # 7. Вивід
  cat("\n=== Full-data normality test results ===\n")
  cat("Kolmogorov–Smirnov test:\n")
  cat("  Fixed:  p =", format.pval(ks_fixed$p.value), "\n")
  cat("  Random: p =", format.pval(ks_random$p.value), "\n\n")
  
  cat("Skewness & Kurtosis:\n")
  cat("  Fixed:  skew =", skew_fixed, ", kurt =", kurt_fixed, "\n")
  cat("  Random: skew =", skew_random, ", kurt =", kurt_random, "\n\n")
  
  cat("Anderson–Darling test:\n")
  cat("  Fixed:  p =", format.pval(ad_fixed$p.value), "\n")
  cat("  Random: p =", format.pval(ad_random$p.value), "\n")
}

analyze_window_to_pdf(
  from_sample = 1000,
  to_sample   = 2000,
  v           = 0.5,
  traces_path = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected",
  inputs_file = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected/inputs.txt",
  plot_dir    = "/Users/andrew/Desktop/protectedvsunprotected/plots_unprot"
)


analyze_window_to_pdf(
  from_sample = 1000,
  to_sample   = 2000,
  v           = 0.5,
  traces_path = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/protected0.5",
  inputs_file = "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/protected0.5/inputs.txt",
  plot_dir    = "/Users/andrew/Desktop/protectedvsunprotected/plots_prot"
)
