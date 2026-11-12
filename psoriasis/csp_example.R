if (!file.exists("csp_functions.R")) {
  stop("csp_functions.R not found. Set your working directory to the folder containing these files.")
}

source("csp_functions.R")

# ============================================================
# Illustrative example dataset (arm-level data)
# ============================================================
#
# Outcome:
#   r = events, n = randomized
#
# Treatment coding:
#   A = Treatment A (Control)
#   B = Treatment B
#   C = Treatment C
#   D = Treatment D
#   E = Treatment E
#
# Studies:
#   S1, S2, S3, S4, S5
# ============================================================

example_data <- data.frame(
  study = c(
    "S1","S1","S1","S1",
    "S2","S2","S2","S2",
    "S3","S3",
    "S4","S4",
    "S5","S5"
  ),
  id = c(
    1,1,1,1,
    2,2,2,2,
    3,3,
    4,4,
    5,5
  ),
  t = c(
    "A","B","C","D",
    "A","C","D","E",
    "A","E",
    "A","B",
    "A","B"
  ),
  r = c(
    1110,482,353,396,
    440,104,148,301,
    75,59,
    91,85,
    33,78
  ),
  n = c(
    4321,2104,1596,1561,
    4088,947,1399,2743,
    521,541,
    148,151,
    58,136
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# Run example
# ============================================================
example_fit <- fit_csp_nma(example_data)

# Optional random-effects analysis
# example_fit_random <- fit_csp_nma(example_data, model = "random")

# 1. Projection matrix (P)
# Section 2.1.3: Projection representation of the NMA estimator
cat("\n================ P MATRIX ================\n")
print(example_fit$P)

# 2. Generalized Cochran Q decomposition
# Q_net = Q_het + Q_inc, where Q_het is within-design heterogeneity
# and Q_inc is between-design inconsistency.
cat("\n================ GENERALIZED COCHRAN Q ================\n")
q_stats <- q_decomposition(example_fit)
print(q_stats, row.names = FALSE)

# 3. Contrast decomposition with estimates and weights
# Section 2.2: Canonical study-based decomposition
cat("\n================ DECOMPOSITION: A:E ================\n")
decomp_AE <- contrast_decomposition_table(example_fit, "A:E")
print(decomp_AE, row.names = FALSE)

# 4. Forest plot for decomposition
# Section 2.3.1: Forest plot for a target comparison
p_forest <- plot_csp_forest(
  example_fit,
  target = "A:E",
  show_indirect_paths = TRUE,
  title = "Illustrative example: forest plot for A:E"
)
print(p_forest)

# 5. Tension plot
# Section 2.3.3: Tension plot: direct versus indirect evidence
p_tension <- plot_csp_tension(
  example_fit,
  baseline = "A",
  title = "Illustrative example: direct / indirect / network comparisons"
)
print(p_tension)

# 6. Three-dimensional visualization of canonical decomposition weights
# Corresponds to Figure 2 in the accompanying manuscript
p_3d <- plot_csp_3d(
  example_fit,
  title = "Illustrative example: direct-study and overall indirect weights"
)
p_3d
