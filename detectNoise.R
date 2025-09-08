library(ggplot2)

traces_path <- "/Users/andrew/Desktop/thesis/only-traces/capture_traces/unprotected"   # folder with trace_*.txt
inputs_file <- "/Users/andrew/Desktop/thesis/only-traces/capture_traces/unprotected/inputs.txt"      # inputs file (N×7)
v            <- 0.5                           # value that defines the "fixed" group (V1 == v)
tvla_thresh  <- 4.5                           # TVLA threshold for |t| (leakage if exceeded)
diff_frac    <- 0.5                           # fraction of the peak |diff_wave| to define a wide window

# Read inputs as a numeric matrix with columns V1..V7
inputs_mat <- as.matrix(
  read.table(inputs_file, header = FALSE, sep = "", col.names = paste0("V", 1:7))
)
N <- nrow(inputs_mat)
cat("Loaded inputs:", N, "rows ×", ncol(inputs_mat), "columns\n")

# Discover and order trace files by numeric index trace_<idx>.txt
trace_files <- list.files(traces_path, "^trace_\\d+\\.txt$", full.names = TRUE)
extract_idx <- function(f) as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(f)))
idx <- sapply(trace_files, extract_idx)
ord <- order(idx)
trace_files <- trace_files[ord]
idx         <- idx[ord]

# Read all traces into a matrix: each row = one trace (scan -> rbind)
trace_list <- lapply(trace_files, scan, quiet = TRUE)
traces     <- do.call(rbind, trace_list)  
S <- ncol(traces)                            # number of samples per trace
cat("Loaded traces:", nrow(traces), "rows ×", S, "samples per trace\n\n")

# Build "fixed" vs "random" groups based on first input value (V1)
first_vals <- inputs_mat[idx + 1, 1]         # align inputs with traces (inputs are 0-based indexed in file naming)
fixed_idx  <- which(first_vals == v)
random_idx <- which(first_vals != v)
cat("Fixed group size :", length(fixed_idx), "\n")
cat("Random group size:", length(random_idx), "\n\n")

# Slice the trace matrix into two groups
fixed_mat  <- traces[fixed_idx, , drop = FALSE]
random_mat <- traces[random_idx, , drop = FALSE]

# Compute average waveforms and their difference
mean_fixed  <- colMeans(fixed_mat)
mean_random <- colMeans(random_mat)
diff_wave   <- mean_fixed - mean_random      # amplitude difference per sample

# Compute Welch t-statistic per sample (column-wise TVLA)
t_values <- sapply(seq_len(S), function(i) {
  t.test(fixed_mat[, i], random_mat[, i], var.equal = FALSE)$statistic
})

# -------- Windowing by "difference wave" (for visualization/heuristics) --------
# Take a wide window around where |diff_wave| is relatively large.
peak_idx <- which.max(abs(diff_wave))        # index of the largest absolute difference
thr_diff <- diff_frac * max(abs(diff_wave))  # threshold = diff_frac * peak |diff|
sel_diff <- which(abs(diff_wave) >= thr_diff)
qs_diff  <- min(sel_diff)                    # start of "diff-based" window
qe_diff  <- max(sel_diff)                    # end of "diff-based" window
cat("diff_wave peak at sample:", peak_idx, "\n")
cat("diff_wave window:", qs_diff, "–", qe_diff, "\n\n")

# -------- Windowing by TVLA (where |t| exceeds the threshold) --------
leaks    <- which(abs(t_values) > tvla_thresh)
qs_tval  <- min(leaks)                       # start of "TVLA-based" window
qe_tval  <- max(leaks)                       # end of "TVLA-based" window
cat("TVLA window (|t|>", tvla_thresh, "):", qs_tval, "–", qe_tval, "\n\n")

# -------- Plots --------
# 1) diff_wave with the diff-based window shown by dashed blue lines
df1 <- data.frame(sample = seq_len(S), diff = diff_wave)
p1 <- ggplot(df1, aes(sample, diff)) +
  geom_line() +
  geom_vline(xintercept = c(qs_diff, qe_diff), linetype="dashed", color="blue") +
  ggtitle("diff_wave (mean_fixed - mean_random)") +
  xlab("Sample index") + ylab("Amplitude difference") +
  theme_minimal()

# 2) t-values with TVLA threshold (red dashed) and TVLA-based window (blue dashed)
df2 <- data.frame(sample = seq_len(S), tval = t_values)
p2 <- ggplot(df2, aes(sample, tval)) +
  geom_line() +
  geom_hline(yintercept = c(-tvla_thresh, tvla_thresh), linetype="dashed", color="red") +
  geom_vline(xintercept = c(qs_tval, qe_tval), linetype="dashed", color="blue") +
  ggtitle("TVLA t-values") +
  xlab("Sample index") + ylab("t-value") +
  theme_minimal()

print(p1)
print(p2)

# Echo the discovered diff-based window in the console (so it prints as values)
qs_diff
qe_diff 
