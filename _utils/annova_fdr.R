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


# ══════════════════════════════════════════════════════════════════════════════
# 参考文献
# ══════════════════════════════════════════════════════════════════════════════
#
# ── ANCOVA (协方差分析) ─────────────────────────────────────────────────────
# Fisher, R. A. (1932). Statistical Methods for Research Workers (4th ed.).
#   Oliver & Boyd.
#   → ANCOVA 将连续协变量纳入线性模型，在控制协变量影响后比较调整均值。
#     模型: Y = β₀ + β₁·Group + β₂·Cov₁ + ... + ε
#     调整均值差 = 各组在协变量均值处的预测值之差。
#
# ── Type II Sum of Squares ──────────────────────────────────────────────────
# Fox, J., & Weisberg, S. (2019). An R Companion to Applied Regression
#   (3rd ed.). Sage. ISBN: 978-1544336473
#   → car::Anova() 默认 Type II SS: 每个效应的 SS 在控制所有同阶或更低阶
#     效应（但不控制更高阶交互）后计算。Type II 不依赖变量输入顺序，
#     适合主效应检验。与 Type I (sequential) 不同，各行 SS 不可加。
#
# ── η² 与偏 η² ─────────────────────────────────────────────────────────────
# Cohen, J. (1973). Eta-squared and partial eta-squared in fixed factor ANOVA
#   designs. Educational and Psychological Measurement, 33(1), 107–112.
#   doi:10.1177/001316447303300111
#   η²         = SS_group / SS_total
#   partial η² = SS_group / (SS_group + SS_residual)
#   → 注意: Type II SS 下 SS_total 不是各行的算术和，因此 η² 可能不等于
#     partial η²。ANCOVA 中通常报告 partial η²。
#
# Richardson, J. T. E. (2011). Eta squared and partial eta squared as measures
#   of effect size in educational research. Educational Research Review,
#   6(2), 135–147. doi:10.1016/j.edurev.2010.12.001
#   → 综述了 η² / partial η² 的区别和报告建议。
#
# ── Cohen's f² ──────────────────────────────────────────────────────────────
# Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences
#   (2nd ed.). Lawrence Erlbaum Associates. ISBN: 978-0805802832
#   公式: f² = R² / (1 - R²) = partial_η² / (1 - partial_η²)
#   阈值: f² < 0.02 small, < 0.15 medium, ≥ 0.15 large
#   → f² 用于 G*Power 等软件的样本量规划。在 ANCOVA 中表示分组变量
#     相对于残差的独特解释力。
#
# ── Cohen's d（三种变体）────────────────────────────────────────────────────
# ① d_raw — 经典 pooled SD 版:
#   Cohen, J. (1988).同上.
#   公式: d = (M₁ - M₂) / s_pooled
#         s_pooled = √(((n₁-1)s₁² + (n₂-1)s₂²) / (n₁+n₂-2))
#   → 不控制协变量，反映原始组间差异。与 ANCOVA 同时报告时，d_raw > d_ancova
#     说明协变量解释了部分组间差异。
#
# ② d_ancova — ANCOVA 调整版:
#   Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
#     (2009). Introduction to Meta-Analysis. Wiley. ISBN: 978-0470057247
#     (第 4 章讨论了协变量调整后的标准化均值差)
#   公式: d = (adjusted_M₁ - adjusted_M₂) / √MSE
#   → √MSE = sigma(model)，即残差标准差。这是 ANCOVA 下推荐的效应量，
#     等价于控制了协变量后的标准化均值差。
#
# ③ d_t — t 值反推:
#   Rosenthal, R. (1994). Parametric measures of effect size. In H. Cooper &
#     L. V. Hedges (Eds.), The Handbook of Research Synthesis (pp. 231–244).
#     Russell Sage Foundation.
#   公式: d = t × √(1/n₁ + 1/n₂)
#   → 从 t 检验统计量反推 d。注意: ANCOVA 中的 t 已包含协变量信息，
#     反推得到的 d_t 通常偏大，仅供与文献对比时参考。推荐使用 d_ancova。
#
# ── FDR (False Discovery Rate) ──────────────────────────────────────────────
# Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
#   a practical and powerful approach to multiple testing. Journal of the
#   Royal Statistical Society: Series B (Methodological), 57(1), 289–300.
#   doi:10.1111/j.2517-6161.1995.tb02031.x
#   → 对 t 检验 p 值和 F 检验 p 值分别进行 FDR 校正。
#   FDR 比 Bonferroni 更不保守（检验效能更高），适合大规模多重比较。
#
# ── Hedges' g（偏差校正 d）──────────────────────────────────────────────────
# Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#   size and related estimators. Journal of Educational Statistics, 6(2),
#   107–128. doi:10.3102/10769986006002107
#   公式: g = d × J(N-2)，其中 J(m) ≈ 1 - 3/(4m-1)
#   → 小样本下 d 略微偏高，g 对其进行了偏差校正。n₁+n₂ ≥ 50 时 d ≈ g。
