# remove_outliers.R — outlier detection (3-sigma / MAD, auto-selected by normality)
# Dependencies: none (base R only)


remove_outliers <- function(data, column, alpha = 0.05,
                            sigma_threshold = 3, mad_threshold = 3,
                            set_na = TRUE, verbose = TRUE) {

  data_copy <- data.frame(data)

  if (!column %in% names(data_copy))
    stop("Column '", column, "' not found in data frame")
  if (!is.numeric(data_copy[[column]]))
    stop("Column '", column, "' is not numeric")

  valid_rows <- which(!is.na(data_copy[[column]]))
  clean_vector <- data_copy[[column]][valid_rows]
  original_length <- length(clean_vector)

  if (original_length < 3)
    stop("Too few data points (less than 3) for normality test")

  # ── Shapiro-Wilk ────────────────────────────────────────────
  if (length(unique(clean_vector)) > 1) {
    shapiro_test <- shapiro.test(clean_vector)
    is_normal <- shapiro_test$p.value > alpha
  } else {
    if (verbose) message("All values identical, skipping: ", column)
    return(data_copy)
  }

  fmt_p <- function(p) {
    if (p < 0.0001) format(p, scientific = TRUE, digits = 3) else round(p, 4)
  }

  if (verbose) {
    cat("=== Outlier Detection Report ===\n")
    cat("Column:", column, "\n")
    cat("Data points:", original_length, "\n")
    cat("Shapiro-Wilk p-value:", fmt_p(shapiro_test$p.value), "\n")
    cat("Normal distribution:", ifelse(is_normal, "Yes", "No"), "\n")
    cat("Action:", ifelse(set_na, "Set outliers to NA",
                          "Remove rows with outliers"), "\n")
  }

  outlier_rows        <- numeric(0)
  outlier_multipliers <- numeric(0)

  if (is_normal) {
    if (verbose) cat("Using 3-sigma method for outlier detection\n")

    mean_val     <- mean(clean_vector)
    sd_val       <- sd(clean_vector)
    lower_bound  <- mean_val - sigma_threshold * sd_val
    upper_bound  <- mean_val + sigma_threshold * sd_val

    is_outlier  <- clean_vector < lower_bound | clean_vector > upper_bound
    outliers    <- clean_vector[is_outlier]
    outlier_multipliers <- abs(outliers - mean_val) / sd_val
    outlier_rows <- valid_rows[is_outlier]

    if (verbose) {
      cat("\n--- Normal Distribution Statistics ---\n")
      cat("Mean:", round(mean_val, 4), "\n")
      cat("SD:",   round(sd_val, 4), "\n")
      cat("Min:",  round(min(clean_vector), 4), "\n")
      cat("Q1:",   round(quantile(clean_vector, 0.25), 4), "\n")
      cat("Median:", round(median(clean_vector), 4), "\n")
      cat("Q3:",   round(quantile(clean_vector, 0.75), 4), "\n")
      cat("Max:",  round(max(clean_vector), 4), "\n")
    }

  } else {
    if (verbose) cat("Using MAD method for outlier detection\n")

    median_val <- median(clean_vector)
    mad_val    <- mad(clean_vector)

    if (mad_val == 0) {
      if (verbose) message("MAD equals 0, skipping: ", column)
      return(data_copy)
    }

    mad_scores  <- abs(clean_vector - median_val) / mad_val
    is_outlier  <- mad_scores > mad_threshold
    outliers    <- clean_vector[is_outlier]
    outlier_multipliers <- mad_scores[is_outlier]
    outlier_rows <- valid_rows[is_outlier]

    if (verbose) {
      cat("\n--- Non-Normal Distribution Statistics ---\n")
      cat("Median:", round(median_val, 4), "\n")
      cat("MAD:",    round(mad_val, 4), "\n")
      cat("Min:",    round(min(clean_vector), 4), "\n")
      cat("Q1:",     round(quantile(clean_vector, 0.25), 4), "\n")
      cat("Q3:",     round(quantile(clean_vector, 0.75), 4), "\n")
      cat("IQR:",    round(IQR(clean_vector), 4), "\n")
      cat("Max:",    round(max(clean_vector), 4), "\n")
    }
  }

  n_outliers <- length(outliers)
  outlier_pct <- round(n_outliers / original_length * 100, 2)

  if (verbose) {
    cat("\n--- Outlier Detection Results ---\n")
    cat("Outliers:", n_outliers, "(", outlier_pct, "%)\n")
  }

  if (n_outliers > 0) {
    if (verbose) {
      cat("Outlier values:", paste(round(outliers, 4), collapse = ", "), "\n")
      cat("Outlier rows:", paste(outlier_rows, collapse = ", "), "\n\n")

      center <- if (is_normal) mean_val else median_val
      unit   <- if (is_normal) "sigma" else "MAD"
      for (i in seq_len(n_outliers)) {
        direction <- ifelse(outliers[i] > center, "above", "below")
        cat(sprintf("  Row %d: Value %.4f is %.2f %s %s the %s\n",
                    outlier_rows[i], outliers[i],
                    outlier_multipliers[i], unit, direction,
                    if (is_normal) "mean" else "median"))
      }
    }

    if (set_na) {
      data_copy[outlier_rows, column] <- NA
    } else {
      data_copy <- data_copy[-outlier_rows, , drop = FALSE]
    }
  }

  if (verbose) cat("====================\n\n")
  return(data_copy)
}


remove_outliers_all <- function(data, alpha = 0.05,
                                sigma_threshold = 3, mad_threshold = 3,
                                set_na = TRUE, col = NULL, verbose = TRUE) {
  if (is.null(col)) col <- colnames(data)
  for (c in col) {
    if (is.numeric(data[[c]])) {
      data <- remove_outliers(
        data = data, column = c,
        alpha = alpha,
        sigma_threshold = sigma_threshold,
        mad_threshold = mad_threshold,
        set_na = set_na,
        verbose = verbose
      )
    }
  }
  return(data)
}


# ══════════════════════════════════════════════════════════════════════════════
# 参考文献
# ══════════════════════════════════════════════════════════════════════════════
#
# ── Shapiro-Wilk 正态性检验 ─────────────────────────────────────────────────
# Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for
#   normality (complete samples). Biometrika, 52(3/4), 591–611.
#   doi:10.1093/biomet/52.3-4.591
#   → 用于自动选择离群值检测方法: 正态 → 3-sigma; 非正态 → MAD。
#     注意: 离群值本身会导致 Shapiro-Wilk 拒绝正态性，因此非正态数据
#     自动走 MAD 路径是合理的设计选择。
#
# ── 3-sigma 法则 (Three-Sigma Rule) ─────────────────────────────────────────
# Pukelsheim, F. (1994). The three sigma rule. The American Statistician,
#   48(2), 88–91. doi:10.1080/00031305.1994.10476030
#   → 对任意分布，至少 88.9% 的数据落在 μ ± 3σ 内；对正态分布约 99.7%。
#     默认阈值 σ = 3 基于正态假设。对于小样本或非正态分布，MAD 方法更稳健。
#
# ── MAD (Median Absolute Deviation) ─────────────────────────────────────────
# Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting
#   outliers: Do not use standard deviation around the mean, use absolute
#   deviation around the median. Journal of Experimental Social Psychology,
#   49(4), 764–766. doi:10.1016/j.jesp.2013.03.013
#   → 推荐 MAD 替代 3-sigma 用于偏态分布。MAD 对离群值本身的污染比 SD
#     更稳健。方法: 计算 |Xᵢ - median|，取其 median 作为 MAD。
#     R 的 mad() 默认 constant = 1.4826 (正态一致性常数)。
#
# Hampel, F. R. (1974). The influence curve and its role in robust estimation.
#   Journal of the American Statistical Association, 69(346), 383–393.
#   doi:10.1080/01621459.1974.10482962
#   → 引入 MAD 作为稳健尺度估计量。
#
# ── IQR 法则 (备选，本代码未直接使用) ────────────────────────────────────────
# Tukey, J. W. (1977). Exploratory Data Analysis. Addison-Wesley.
#   ISBN: 978-0201076165
#   → 离群值定义为 < Q1 - 1.5×IQR 或 > Q3 + 1.5×IQR。MAD 通常比 IQR 法则
#     对小样本更敏感。
#
# ── 缺失值处理策略 ─────────────────────────────────────────────────────────
# Schafer, J. L., & Graham, J. W. (2002). Missing data: our view of the state
#   of the art. Psychological Methods, 7(2), 147–177.
#   doi:10.1037/1082-989X.7.2.147
#   → set_na = TRUE: 将离群值设为 NA 而非删除整行，保留样本量。
#     set_na = FALSE: 删除离群值所在行。前者适合后续多重插补，后者适合
#     离群值比例很低时的简单处理。
