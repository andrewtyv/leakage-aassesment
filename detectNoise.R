library(ggplot2)

traces_path <- "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected"   # папка з trace_*.txt
inputs_file <- "/Users/andrew/Desktop/protectedvsunprotected/only-traces/capture_traces/unprotected/inputs.txt"      # файл з входами (N×7)
v            <- 0.5                           
tvla_thresh  <- 4.5                             
diff_frac    <- 0.5                          
inputs_mat <- as.matrix(
  read.table(inputs_file, header = FALSE, sep = "", col.names = paste0("V", 1:7))
)
N <- nrow(inputs_mat)
cat("Loaded inputs:", N, "rows ×", ncol(inputs_mat), "columns\n")

trace_files <- list.files(traces_path, "^trace_\\d+\\.txt$", full.names = TRUE)
extract_idx <- function(f) as.integer(sub("^trace_(\\d+)\\.txt$", "\\1", basename(f)))
idx <- sapply(trace_files, extract_idx)
ord <- order(idx)
trace_files <- trace_files[ord]
idx         <- idx[ord]

trace_list <- lapply(trace_files, scan, quiet = TRUE)
traces     <- do.call(rbind, trace_list)  
S <- ncol(traces)
cat("Loaded traces:", nrow(traces), "rows ×", S, "samples per trace\n\n")


first_vals <- inputs_mat[idx + 1, 1]
fixed_idx  <- which(first_vals == v)
random_idx <- which(first_vals != v)
cat("Fixed group size :", length(fixed_idx), "\n")
cat("Random group size:", length(random_idx), "\n\n")

fixed_mat  <- traces[fixed_idx, , drop = FALSE]
random_mat <- traces[random_idx, , drop = FALSE]

mean_fixed  <- colMeans(fixed_mat)
mean_random <- colMeans(random_mat)
diff_wave   <- mean_fixed - mean_random

t_values <- sapply(seq_len(S), function(i) {
  t.test(fixed_mat[, i], random_mat[, i], var.equal = FALSE)$statistic
})

peak_idx <- which.max(abs(diff_wave))
thr_diff <- diff_frac * max(abs(diff_wave))
sel_diff <- which(abs(diff_wave) >= thr_diff)
qs_diff  <- min(sel_diff)
qe_diff  <- max(sel_diff)
cat("diff_wave peak at sample:", peak_idx, "\n")
cat("diff_wave window:", qs_diff, "–", qe_diff, "\n\n")

leaks    <- which(abs(t_values) > tvla_thresh)
qs_tval  <- min(leaks)
qe_tval  <- max(leaks)
cat("TVLA window (|t|>", tvla_thresh, "):", qs_tval, "–", qe_tval, "\n\n")


df1 <- data.frame(sample = seq_len(S), diff = diff_wave)
p1 <- ggplot(df1, aes(sample, diff)) +
  geom_line() +
  geom_vline(xintercept = c(qs_diff, qe_diff), linetype="dashed", color="blue") +
  ggtitle("diff_wave (mean_fixed - mean_random)") +
  xlab("Sample index") + ylab("Amplitude difference") +
  theme_minimal()

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

qs_diff
qe_diff