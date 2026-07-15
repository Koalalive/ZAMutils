# annova_fdr.R — ANCOVA with FDR correction and effect sizes
# Dependencies: tidyverse, broom, effectsize, car

library(tidyverse)
library(broom)
library(effectsize)
library(car)


annova_FDR <- function(data, group_var, covariates = NULL, outcome_vars = NULL) {

  results_list <- list()

  for (outcome in outcome_vars) {

    # ── 构建公式 ──────────────────────────────────────────────
    if (is.null(covariates) || length(covariates) == 0) {
      formula_str <- paste(outcome, "~", group_var)
    } else {
      formula_str <- paste(outcome, "~", group_var, "+",
                           paste(covariates, collapse = " + "))
    }
    formula_obj <- as.formula(formula_str)

    model <- lm(formula_obj, data = data)
    model_summary <- summary(model)

    # ── ANOVA 表 (Type II, car::Anova) ────────────────────────
    anova_table <- car::Anova(model)

    F_row        <- anova_table[group_var, ]
    F_value      <- F_row[["F value"]]
    df_numerator <- F_row[["Df"]]
    F_p_value    <- F_row[["Pr(>F)"]]
    ss_group     <- F_row[["Sum Sq"]]

    df_residual <- anova_table["Residuals", "Df"]
    ss_residual <- anova_table["Residuals", "Sum Sq"]
    ss_total    <- sum(anova_table[, "Sum Sq"])

    # ── 效应量 ────────────────────────────────────────────────
    eta_squared    <- ss_group / ss_total
    partial_eta_sq <- ss_group / (ss_group + ss_residual)
    cohens_f2      <- partial_eta_sq / (1 - partial_eta_sq)

    # ── 各组样本量 ────────────────────────────────────────────
    group_n <- table(data[[group_var]])
    residual_sd <- sigma(model)

    # ── 提取所有分组系数 (vs 参考水平) ────────────────────────
    coef_df    <- as.data.frame(model_summary$coefficients)
    coef_names <- rownames(coef_df)

    group_coef_names <- grep(paste0("^", group_var), coef_names, value = TRUE)

    if (length(group_coef_names) == 0) {
      warning("No group coefficients found for outcome: ", outcome)
      next
    }

    group_levels <- levels(as.factor(data[[group_var]]))
    ref_level <- group_levels[1]

    for (coef_name in group_coef_names) {
      gc <- coef_df[coef_name, ]

      # 从系数名推断比较的组 (e.g. group_varDrugB → DrugB)
      comp_level <- sub(paste0("^", group_var), "", coef_name)

      adj_mean_diff <- gc[["Estimate"]]
      adj_se        <- gc[["Std. Error"]]
      t_val         <- gc[["t value"]]
      p_val         <- gc[["Pr(>|t|)"]]
      t_df          <- model_summary$df[2]

      n_ref  <- unname(group_n[ref_level])
      n_comp <- unname(group_n[comp_level])

      # ── Cohen's d 变体 ──────────────────────────────────

      # ① d_ancova: ANCOVA 修正 d（控制协变量后的标准化均值差）
      d_ancova <- NA_real_
      if (!is.na(adj_mean_diff) && !is.na(residual_sd) && residual_sd > 0) {
        d_ancova <- adj_mean_diff / residual_sd
      }

      # ② d_raw: 经典 d（原始均值差 / 原始 pooled SD，不控制协变量）
      d_raw <- NA_real_
      raw_ref  <- na.omit(data[[outcome]][data[[group_var]] == ref_level])
      raw_comp <- na.omit(data[[outcome]][data[[group_var]] == comp_level])
      if (length(raw_ref) > 1 && length(raw_comp) > 1) {
        n1 <- length(raw_comp)
        n2 <- length(raw_ref)
        s1 <- var(raw_comp)
        s2 <- var(raw_ref)
        pooled_sd_raw <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
        if (pooled_sd_raw > 0) {
          d_raw <- (mean(raw_comp) - mean(raw_ref)) / pooled_sd_raw
        }
      }

      # ③ d_t: t 值反推（ANCOVA 下 t 已控制协变量，此值仅供参考）
      d_t <- NA_real_
      if (!is.na(t_val) && !is.na(n_ref) && !is.na(n_comp)) {
        d_t <- t_val * sqrt(1 / n_ref + 1 / n_comp)
      }

      # ── 构建输出 ────────────────────────────────────────
      results_list[[paste0(outcome, "|", coef_name)]] <- data.frame(
        Outcome        = outcome,
        Comparison     = paste(comp_level, "-", ref_level),
        Ref_N          = n_ref,
        Comp_N         = n_comp,

        Adj_Mean_Diff  = adj_mean_diff,
        SE             = adj_se,

        t_value        = t_val,
        t_df           = t_df,
        t_p_value      = p_val,

        F_value        = F_value,
        F_df_num       = df_numerator,
        F_df_den       = df_residual,
        F_p_value      = F_p_value,

        d_ancova       = d_ancova,
        d_raw          = d_raw,
        d_t            = d_t,
        Cohens_f2      = cohens_f2,
        Eta_Squared    = eta_squared,
        Partial_Eta_Sq = partial_eta_sq,

        Residual_df    = df_residual,
        Residual_SD    = residual_sd,
        R_squared      = model_summary$r.squared,
        Adj_R_squared  = model_summary$adj.r.squared,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results_list) == 0) {
    warning("No results to return.")
    return(NULL)
  }

  results_df <- bind_rows(results_list)

  # ── FDR 校正 ────────────────────────────────────────────────
  results_df$FDR_t_p   <- p.adjust(results_df$t_p_value, method = "fdr")
  results_df$Sig_FDR_t <- results_df$FDR_t_p < 0.05

  results_df$FDR_F_p   <- p.adjust(results_df$F_p_value, method = "fdr")
  results_df$Sig_FDR_F <- results_df$FDR_F_p < 0.05

  # ── 排列 ────────────────────────────────────────────────────
  results_df <- results_df %>%
    select(
      Outcome, Comparison, Ref_N, Comp_N,
      Adj_Mean_Diff, SE,
      t_value, t_df, t_p_value, FDR_t_p, Sig_FDR_t,
      F_value, F_df_num, F_df_den, F_p_value, FDR_F_p, Sig_FDR_F,
      d_ancova, d_raw, d_t, Cohens_f2, Eta_Squared, Partial_Eta_Sq,
      Residual_df, Residual_SD, R_squared, Adj_R_squared
    )

  return(results_df)
}
