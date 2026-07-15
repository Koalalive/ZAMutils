# group_comparison.R — 组间比较（参数/非参数自适应）
# Dependencies: tidyverse, rstatix, PMCMRplus

library(tidyverse)
library(rstatix)
library(PMCMRplus)


# ── 控制台输出版 ──────────────────────────────────────────────────────────────

group_comparison <- function(data, colname) {
  complete_data <- data %>%
    dplyr::select(group, !!sym(colname)) %>%
    na.omit()

  n_complete <- nrow(complete_data)
  group_counts <- complete_data %>%
    count(group) %>%
    pull(n)

  cat(paste("Complete cases for [", colname, "]:", n_complete, "\n"))
  cat(paste("Group sample sizes:", paste(group_counts, collapse = ", "), "\n\n"))

  cat(paste("Descriptive Statistics for [", colname, "]:\n", sep = ""))
  desc_stats <- complete_data %>%
    group_by(group) %>%
    summarise(
      n      = n(),
      mean   = mean(!!sym(colname), na.rm = TRUE),
      sd     = sd(!!sym(colname), na.rm = TRUE),
      Q1     = quantile(!!sym(colname), 0.25, na.rm = TRUE),
      median = median(!!sym(colname), na.rm = TRUE),
      Q3     = quantile(!!sym(colname), 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  print(desc_stats)
  cat("\n")

  form <- reformulate("group", response = colname)
  parametric_test <- TRUE

  # ── Shapiro-Wilk ──────────────────────────────────────────
  cat(paste("Shapiro Wilk test for [", colname, "]:\n", sep = ""))
  p_shapiro <- complete_data %>%
    group_by(group) %>%
    shapiro_test(!!sym(colname)) %>%
    print() %>%
    dplyr::select(any_of(c("p", "p.value")))

  if (all(p_shapiro$p >= 0.05)) {
    cat("Normality test passed.\n")
  } else {
    cat("Normality test failed.\n")
    parametric_test <- FALSE
  }

  # ── Levene ────────────────────────────────────────────────
  cat(paste("Levene test for [", colname, "]:\n", sep = ""))
  p_leven <- complete_data %>%
    levene_test(form) %>%
    print() %>%
    dplyr::select(any_of(c("p", "p.value")))

  if (p_leven$p >= 0.05) {
    cat("Variance homogeneity test passed.\n")
  } else {
    cat("Variance homogeneity test failed.\n")
    parametric_test <- FALSE
  }

  k <- n_distinct(complete_data$group)
  df_between <- k - 1
  df_within  <- n_complete - k

  if (parametric_test) {
    cat("Parametric test is selected.\n")
    cat(paste("One-way ANOVA for [", colname, "]:\n", sep = ""))

    anova_result <- complete_data %>%
      anova_test(form) %>%
      print()

    ss_between  <- anova_result$SSn[1]
    ss_total    <- anova_result$SSn[1] + anova_result$SSd[1]
    eta_squared <- ss_between / ss_total

    cat(paste("Eta squared (η²) =", round(eta_squared, 4), "\n"))
    cat(paste("Degrees of freedom: between =", df_between, ", within =", df_within, "\n"))

    if (anova_result$p[1] < 0.05 && k > 2) {
      cat(paste("Tukey and Bonferroni test for [", colname, "]:\n", sep = ""))

      complete_data %>%
        tukey_hsd(form) %>%
        print()

      bonferroni_result <- complete_data %>%
        t_test(form, p.adjust.method = "bonferroni") %>%
        print()

      if (nrow(bonferroni_result) > 0) {
        cat("\nEffect sizes (r) for pairwise comparisons:\n")
        for (i in seq_len(nrow(bonferroni_result))) {
          g1 <- complete_data %>%
            filter(group == bonferroni_result$group1[i]) %>%
            pull(!!sym(colname))
          g2 <- complete_data %>%
            filter(group == bonferroni_result$group2[i]) %>%
            pull(!!sym(colname))

          mean_diff <- mean(g1) - mean(g2)
          pooled_sd <- sqrt(
            ((length(g1) - 1) * var(g1) + (length(g2) - 1) * var(g2)) /
            (length(g1) + length(g2) - 2)
          )
          cohens_d <- mean_diff / pooled_sd
          r_effect  <- cohens_d / sqrt(cohens_d^2 + 4)

          cat(paste(
            bonferroni_result$group1[i], "vs", bonferroni_result$group2[i],
            ": r =", round(r_effect, 4), "\n"
          ))
        }
      }
    }

  } else {
    cat("Nonparametric test is selected.\n")
    cat(paste("Kruskal-Wallis test for [", colname, "]:\n", sep = ""))

    kw_result <- complete_data %>%
      kruskal_test(form) %>%
      print()

    H_stat  <- kw_result$statistic
    n_total <- n_complete
    epsilon_sq <- max(0, (H_stat - (k - 1)) / (n_total - k))

    cat(paste("Epsilon squared (ε²) =", round(epsilon_sq, 4), "\n"))
    cat(paste("Degrees of freedom =", df_between, "\n"))

    if (kw_result$p < 0.05 && k > 2) {
      cat(paste("Dunn test for [", colname, "]:\n", sep = ""))

      dunn_result <- complete_data %>%
        dunn_test(form) %>%
        print()

      if (nrow(dunn_result) > 0) {
        cat("\nEffect sizes (r) for Dunn test comparisons:\n")
        for (i in seq_len(nrow(dunn_result))) {
          g1 <- complete_data %>%
            filter(group == dunn_result$group1[i]) %>%
            pull(!!sym(colname))
          g2 <- complete_data %>%
            filter(group == dunn_result$group2[i]) %>%
            pull(!!sym(colname))

          z_value  <- dunn_result$statistic[i]
          r_effect <- z_value / sqrt(length(g1) + length(g2))

          cat(paste(
            dunn_result$group1[i], "vs", dunn_result$group2[i],
            ": Z =", round(z_value, 4),
            ", r =", round(r_effect, 4), "\n"
          ))
        }
      }
    }
  }
}


# ── 对象返回版 ────────────────────────────────────────────────────────────────

group_comparison_obj <- function(data, colname) {
  results <- list(
    variable          = colname,
    complete_cases    = NULL,
    group_counts      = NULL,
    descriptive_stats = NULL,
    normality_test    = NULL,
    homogeneity_test  = NULL,
    parametric        = NULL,
    main_test         = NULL,
    effect_size       = NULL,
    posthoc_tests     = NULL,
    pairwise_effects  = NULL
  )

  complete_data <- data %>%
    dplyr::select(group, !!sym(colname)) %>%
    na.omit()

  n_complete <- nrow(complete_data)

  # deframe 保证组名和样本量对齐
  results$complete_cases <- n_complete
  results$group_counts <- complete_data %>%
    count(group) %>%
    tibble::deframe()

  results$descriptive_stats <- complete_data %>%
    group_by(group) %>%
    summarise(
      n      = n(),
      mean   = mean(!!sym(colname), na.rm = TRUE),
      sd     = sd(!!sym(colname), na.rm = TRUE),
      Q1     = quantile(!!sym(colname), 0.25, na.rm = TRUE),
      median = median(!!sym(colname), na.rm = TRUE),
      Q3     = quantile(!!sym(colname), 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  form <- reformulate("group", response = colname)
  parametric_test <- TRUE

  # ── Shapiro-Wilk ──────────────────────────────────────────
  p_shapiro <- complete_data %>%
    group_by(group) %>%
    rstatix::shapiro_test(!!sym(colname))

  results$normality_test <- list(
    test    = "Shapiro-Wilk",
    results = p_shapiro,
    passed  = all(p_shapiro$p >= 0.05)
  )

  if (!all(p_shapiro$p >= 0.05)) parametric_test <- FALSE

  # ── Levene ────────────────────────────────────────────────
  p_leven <- complete_data %>%
    rstatix::levene_test(form)

  results$homogeneity_test <- list(
    test    = "Levene",
    results = p_leven,
    passed  = p_leven$p >= 0.05
  )

  if (p_leven$p < 0.05) parametric_test <- FALSE

  results$parametric <- parametric_test

  k <- n_distinct(complete_data$group)
  df_between <- k - 1
  df_within  <- n_complete - k

  if (parametric_test) {
    # ── One-way ANOVA ───────────────────────────────────────
    anova_result <- complete_data %>%
      rstatix::anova_test(form)

    ss_between  <- anova_result$SSn[1]
    ss_total    <- anova_result$SSn[1] + anova_result$SSd[1]
    eta_squared <- ss_between / ss_total

    results$main_test <- list(
      test       = "One-way ANOVA",
      results    = anova_result,
      df_between = df_between,
      df_within  = df_within
    )

    results$effect_size <- list(
      eta_squared    = eta_squared,
      interpretation = ifelse(eta_squared < 0.01, "negligible",
                       ifelse(eta_squared < 0.06, "small",
                       ifelse(eta_squared < 0.14, "medium", "large")))
    )

    if (anova_result$p[1] < 0.05 && k > 2) {
      tukey_result <- complete_data %>%
        rstatix::tukey_hsd(form)

      bonferroni_result <- complete_data %>%
        rstatix::t_test(form, p.adjust.method = "bonferroni")

      results$posthoc_tests <- list(
        tukey      = tukey_result,
        bonferroni = bonferroni_result
      )

      if (nrow(bonferroni_result) > 0) {
        pairwise_effects <- vector("list", nrow(bonferroni_result))

        for (i in seq_len(nrow(bonferroni_result))) {
          g1 <- complete_data %>%
            filter(group == bonferroni_result$group1[i]) %>%
            pull(!!sym(colname))
          g2 <- complete_data %>%
            filter(group == bonferroni_result$group2[i]) %>%
            pull(!!sym(colname))

          mean_diff <- mean(g1) - mean(g2)
          pooled_sd <- sqrt(
            ((length(g1) - 1) * var(g1) + (length(g2) - 1) * var(g2)) /
            (length(g1) + length(g2) - 2)
          )
          cohens_d <- mean_diff / pooled_sd
          r_effect  <- cohens_d / sqrt(cohens_d^2 + 4)

          pairwise_effects[[i]] <- data.frame(
            group1     = bonferroni_result$group1[i],
            group2     = bonferroni_result$group2[i],
            cohens_d   = cohens_d,
            r_effect   = r_effect,
            interpretation = ifelse(abs(cohens_d) < 0.2, "negligible",
                             ifelse(abs(cohens_d) < 0.5, "small",
                             ifelse(abs(cohens_d) < 0.8, "medium", "large"))),
            stringsAsFactors = FALSE
          )
        }

        results$pairwise_effects <- bind_rows(pairwise_effects)
      }
    }

  } else {
    # ── Kruskal-Wallis ──────────────────────────────────────
    kw_result <- complete_data %>%
      rstatix::kruskal_test(form)

    H_stat  <- kw_result$statistic
    n_total <- n_complete
    epsilon_sq <- max(0, (H_stat - (k - 1)) / (n_total - k))

    results$main_test <- list(
      test    = "Kruskal-Wallis",
      results = kw_result,
      df      = df_between
    )

    results$effect_size <- list(
      epsilon_squared = epsilon_sq,
      interpretation  = ifelse(epsilon_sq < 0.01, "negligible",
                        ifelse(epsilon_sq < 0.06, "small",
                        ifelse(epsilon_sq < 0.14, "medium", "large")))
    )

    if (kw_result$p < 0.05 && k > 2) {
      dunn_result <- complete_data %>%
        rstatix::dunn_test(form)

      results$posthoc_tests <- list(dunn = dunn_result)

      if (nrow(dunn_result) > 0) {
        pairwise_effects <- vector("list", nrow(dunn_result))

        for (i in seq_len(nrow(dunn_result))) {
          g1 <- complete_data %>%
            filter(group == dunn_result$group1[i]) %>%
            pull(!!sym(colname))
          g2 <- complete_data %>%
            filter(group == dunn_result$group2[i]) %>%
            pull(!!sym(colname))

          z_value  <- dunn_result$statistic[i]
          r_effect <- z_value / sqrt(length(g1) + length(g2))

          pairwise_effects[[i]] <- data.frame(
            group1     = dunn_result$group1[i],
            group2     = dunn_result$group2[i],
            z_value    = z_value,
            r_effect   = r_effect,
            interpretation = ifelse(abs(r_effect) < 0.1, "negligible",
                             ifelse(abs(r_effect) < 0.3, "small",
                             ifelse(abs(r_effect) < 0.5, "medium", "large"))),
            stringsAsFactors = FALSE
          )
        }

        results$pairwise_effects <- bind_rows(pairwise_effects)
      }
    }
  }

  # ── 打印方法 ─────────────────────────────────────────────
  results$print_summary <- function() {
    cat(paste("=== Group Comparison Results for [", results$variable, "] ===\n\n"))
    cat(paste("Complete cases:", results$complete_cases, "\n"))
    cat(paste("Group sample sizes:", paste(results$group_counts, collapse = ", "), "\n\n"))
    cat("Descriptive Statistics:\n")
    print(results$descriptive_stats)
    cat("\n")
    cat("Normality Test (Shapiro-Wilk):\n")
    print(results$normality_test$results)
    cat(paste("Passed:", results$normality_test$passed, "\n\n"))
    cat("Homogeneity of Variance Test (Levene):\n")
    print(results$homogeneity_test$results)
    cat(paste("Passed:", results$homogeneity_test$passed, "\n\n"))
    cat(paste("Test selected:", ifelse(results$parametric, "Parametric", "Nonparametric"), "\n"))
    cat(paste("Main test:", results$main_test$test, "\n"))
    print(results$main_test$results)

    if (results$parametric) {
      cat(paste("\nEffect size (η²):", round(results$effect_size$eta_squared, 4),
                " [", results$effect_size$interpretation, "]\n"))
    } else {
      cat(paste("\nEffect size (ε²):", round(results$effect_size$epsilon_squared, 4),
                " [", results$effect_size$interpretation, "]\n"))
    }

    if (!is.null(results$posthoc_tests)) {
      cat("\nPost-hoc Tests:\n")
      if (results$parametric) {
        cat("Tukey HSD:\n")
        print(results$posthoc_tests$tukey)
        cat("\nBonferroni:\n")
        print(results$posthoc_tests$bonferroni)
      } else {
        cat("Dunn test:\n")
        print(results$posthoc_tests$dunn)
      }
      if (!is.null(results$pairwise_effects)) {
        cat("\nPairwise Effect Sizes:\n")
        print(results$pairwise_effects)
      }
    }
  }

  class(results) <- c("group_comparison_results", "list")
  return(results)
}
