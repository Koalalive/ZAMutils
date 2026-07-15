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
