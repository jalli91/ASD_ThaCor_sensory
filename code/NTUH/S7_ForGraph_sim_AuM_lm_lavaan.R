# =============================================================================
# Simple Mediation: Thalamo-Auditory FC → AASP Auditory → SRS fso_com
# Model 4 (no moderator): X (FC) → M (auditory_c) → Y (fso_com)
#                         Covariates: age_c, meanFD_c (within-sex subgroup)
# Subgroup: females (sex_F = 1) only
#
# SE strategy (PROCESS-equivalent, v4):
#
#   a-path  : lm(M ~ X + covs)           → OLS SE, t-distribution CI and p
#   b-path  : lm(Y ~ M + X + covs)       → OLS SE, t-distribution CI and p
#   c′-path : lm(Y ~ M + X + covs)       → OLS SE, t-distribution CI and p
#   total c : lm(Y ~ X + covs)           → OLS SE, t-distribution CI and p
#   indirect: lavaan sem(se="bootstrap")  → bootstrap SD + percentile CI
#
#   This exactly replicates the PROCESS macro (Hayes, 2022) approach:
#   OLS regression for individual paths, bootstrap only for indirect.
#   Point estimates are identical across lm and lavaan; SE now matches lm.
#
# v4 fixes (relative to v3):
#   [Fix 1] BH-FDR now applied to b-path as well (8 tests, one per FC),
#           matching the treatment of a, cp, and total.
#   [Fix 2] lavaan sem() uses fixed.x = TRUE (covariates treated as fixed
#           exogenous variables, consistent with PROCESS macro and lm()).
#           Reduces bootstrap convergence failures in small subgroups.
#   [Fix 3] Decomposition sanity check added: warns if |c - (c'+ab)| > 1e-6,
#           which would indicate inconsistent complete-cases across models.
#   [Fix 5] Reverted: Z-score standardization is performed on the full sample
#           (before subsetting) so that male and female z-scores remain on the
#           same scale — intentional design for cross-sex comparability.
# =============================================================================

# ---- Section 1: Setup -------------------------------------------------------

pkgs <- c("readxl", "lavaan", "dplyr", "tidyr", "ggplot2", "writexl")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

set.seed(42)

data_file <- "/Users/james/Desktop/CYL/ASD_ThaCor/Score_sheet/drop_nosig_raw/Combine_All_cleaned_0.01-0.08_drop_nosig_raw_AuTha.xlsx"
out_dir   <- "/Users/james/Desktop/CYL/ASD_ThaCor/Score_sheet/drop_nosig_raw/ForGraph_sim_AuM_lm_lavaan"

# ---- Section 2: Data loading and preparation --------------------------------

raw <- read_excel(data_file, sheet = "工作表1")

rename_map <- c(
  sex_M1F0     = "gender(M1F0)",
  fc_thR_aSTGr = "Default_atlas.Thalamus r*Default_atlas.aSTG r (Superior Temporal Gyrus, anterior division Right)",
  fc_thR_aSTGl = "Default_atlas.Thalamus r*Default_atlas.aSTG l (Superior Temporal Gyrus, anterior division Left)",
  fc_thR_pSTGr = "Default_atlas.Thalamus r*Default_atlas.pSTG r (Superior Temporal Gyrus, posterior division Right)",
  fc_thR_pSTGl = "Default_atlas.Thalamus r*Default_atlas.pSTG l (Superior Temporal Gyrus, posterior division Left)",
  fc_thL_aSTGr = "Default_atlas.Thalamus l*Default_atlas.aSTG r (Superior Temporal Gyrus, anterior division Right)",
  fc_thL_aSTGl = "Default_atlas.Thalamus l*Default_atlas.aSTG l (Superior Temporal Gyrus, anterior division Left)",
  fc_thL_pSTGr = "Default_atlas.Thalamus l*Default_atlas.pSTG r (Superior Temporal Gyrus, posterior division Right)",
  fc_thL_pSTGl = "Default_atlas.Thalamus l*Default_atlas.pSTG l (Superior Temporal Gyrus, posterior division Left)"
)

df <- raw |>
  rename(any_of(rename_map)) |>
  mutate(
    sex_F = 1L - as.integer(sex_M1F0),
    group = case_when(
      substr(ID, 2, 2) == "2" ~ "ASD",
      substr(ID, 2, 2) == "1" ~ "NAC",
      TRUE ~ NA_character_
    )
  )

fc_cols <- c("fc_thR_aSTGr","fc_thR_aSTGl","fc_thR_pSTGr","fc_thR_pSTGl",
             "fc_thL_aSTGr","fc_thL_aSTGl","fc_thL_pSTGr","fc_thL_pSTGl")

# ---- Section 3: Z-score standardization -------------------------------------
# Standardization is performed on the full sample (before sex subsetting),
# so that z-scores are on a common scale across sexes. This is intentional:
# male and female z-scores remain comparable on the same metric.

df <- df |>
  mutate(
    auditory_c = as.numeric(scale(auditory)),
    age_c      = as.numeric(scale(age)),
    meanFD_c   = as.numeric(scale(meanFD)),
    fso_com_z  = as.numeric(scale(fso_com))
  )

for (fc in fc_cols) {
  df[[paste0(fc, "_c")]] <- as.numeric(scale(df[[fc]]))
}

fc_c_cols <- paste0(fc_cols, "_c")

# ---- Section 4: PROCESS-equivalent mediation function ----------------------
#
# Path estimation strategy (mirrors PROCESS macro exactly):
#
#   lm_a    : M ~ X + covariates
#             → a-path: coef[X], SE, t, p, 95% t-CI
#
#   lm_y    : Y ~ M + X + covariates
#             → b-path : coef[M], SE, t, p, 95% t-CI
#             → c′-path: coef[X], SE, t, p, 95% t-CI
#
#   lm_c    : Y ~ X + covariates           (no mediator)
#             → total c: coef[X], SE, t, p, 95% t-CI
#
#   lavaan  : same two-equation model, se = "bootstrap", R = 5000
#             → indirect (a×b): bootstrap SD + percentile CI only
#               point estimate also taken from here (= a_lm × b_lm)
#
# All p-values for a, b, c′, total are Wald t-tests from lm().
# FDR-BH correction applied across 8 FCs for a, b, c′, total.
# indirect significance → bootstrap CI excludes zero.

run_simple_med <- function(data, x_var, y_var,
                           m_var      = "auditory_c",
                           covariates = c("age_c", "meanFD_c"),
                           R          = 5000) {

  all_vars <- c(y_var, m_var, x_var, covariates)
  d        <- data[complete.cases(data[, all_vars, drop = FALSE]), ]
  n_used   <- nrow(d)

  if (n_used < length(all_vars) + 2) {
    message("  Skipping ", x_var, ": n too small (", n_used, ")")
    return(NULL)
  }

  cov_rhs <- if (length(covariates) > 0)
    paste("+", paste(covariates, collapse = " + "))
  else ""

  # ------------------------------------------------------------------
  # OLS regressions — individual path SE (identical to PROCESS macro)
  # ------------------------------------------------------------------

  # a-path: M ~ X + covariates
  fml_a  <- as.formula(paste(m_var,  "~", x_var,        cov_rhs))
  fit_a  <- lm(fml_a, data = d)
  sm_a   <- summary(fit_a)$coefficients
  ci_a   <- confint(fit_a, level = 0.95)

  a_est  <- sm_a[x_var, "Estimate"]
  a_se   <- sm_a[x_var, "Std. Error"]
  a_t    <- sm_a[x_var, "t value"]
  a_p    <- sm_a[x_var, "Pr(>|t|)"]
  a_cilo <- ci_a[x_var, 1]
  a_cihi <- ci_a[x_var, 2]

  # Y model: Y ~ M + X + covariates  → b and c′
  fml_y  <- as.formula(paste(y_var, "~", m_var, "+", x_var, cov_rhs))
  fit_y  <- lm(fml_y, data = d)
  sm_y   <- summary(fit_y)$coefficients
  ci_y   <- confint(fit_y, level = 0.95)

  b_est  <- sm_y[m_var,  "Estimate"];   b_se  <- sm_y[m_var,  "Std. Error"]
  b_t    <- sm_y[m_var,  "t value"];    b_p   <- sm_y[m_var,  "Pr(>|t|)"]
  b_cilo <- ci_y[m_var,  1];            b_cihi <- ci_y[m_var,  2]

  cp_est <- sm_y[x_var,  "Estimate"];   cp_se  <- sm_y[x_var,  "Std. Error"]
  cp_t   <- sm_y[x_var,  "t value"];    cp_p   <- sm_y[x_var,  "Pr(>|t|)"]
  cp_cilo <- ci_y[x_var, 1];            cp_cihi <- ci_y[x_var, 2]

  # Total effect: Y ~ X + covariates  (no mediator)
  fml_c  <- as.formula(paste(y_var, "~", x_var, cov_rhs))
  fit_c  <- lm(fml_c, data = d)
  sm_c   <- summary(fit_c)$coefficients
  ci_c   <- confint(fit_c, level = 0.95)

  c_est  <- sm_c[x_var, "Estimate"];    c_se  <- sm_c[x_var, "Std. Error"]
  c_t    <- sm_c[x_var, "t value"];     c_p   <- sm_c[x_var, "Pr(>|t|)"]
  c_cilo <- ci_c[x_var, 1];             c_cihi <- ci_c[x_var, 2]

  # Indirect point estimate from OLS products (= lavaan point est)
  ab_est <- a_est * b_est

  # ------------------------------------------------------------------
  # lavaan bootstrap — indirect percentile CI only
  # ------------------------------------------------------------------
  cov_str <- if (length(covariates) > 0)
    paste(" +", paste(covariates, collapse = " + "))
  else ""

  model_syn <- paste0(
    m_var, " ~ a*", x_var, cov_str, "\n",
    y_var, " ~ b*", m_var, " + cp*", x_var, cov_str, "\n",
    "indirect := a * b\n"
  )

  fit_boot <- tryCatch(
    lavaan::sem(model_syn, data = d,
                se = "bootstrap", bootstrap = R,
                fixed.x = TRUE,   # [Fix 2] treat covariates as fixed exogenous
                                  # (consistent with PROCESS macro / lm(); reduces
                                  #  bootstrap convergence failures in small n)
                missing = "listwise"),
    error = function(e) {
      message("  [boot] lavaan error for ", x_var, ": ", e$message)
      NULL
    }
  )

  # Bootstrap CI for indirect
  ab_boot_se <- NA_real_
  ab_cilo    <- NA_real_
  ab_cihi    <- NA_real_
  ab_sig     <- NA
  n_boot_fail <- NA_integer_

  if (!is.null(fit_boot)) {
    n_boot_fail <- tryCatch({
      b <- lavTech(fit_boot, "boot")
      if (is.matrix(b) && nrow(b) > 0) R - nrow(b) else NA_integer_
    }, error = function(e) NA_integer_)

    pe_boot <- tryCatch(
      lavaan::parameterestimates(fit_boot,
                                 boot.ci.type = "perc",
                                 level        = 0.95,
                                 ci           = TRUE,
                                 standardized = FALSE),
      error = function(e) NULL
    )

    if (!is.null(pe_boot)) {
      row_ab <- pe_boot[pe_boot$op == ":=" & pe_boot$lhs == "indirect", ]
      if (nrow(row_ab) > 0) {
        ab_boot_se <- row_ab$se[1]
        ab_cilo    <- row_ab$ci.lower[1]
        ab_cihi    <- row_ab$ci.upper[1]
        ab_sig     <- (ab_cilo > 0) | (ab_cihi < 0)
      }
    }
  }

  # ------------------------------------------------------------------
  # Assemble output: one row per path
  # ------------------------------------------------------------------
  make_path <- function(label, est, se, t_val, p_val, cilo, cihi,
                        sig_ci = NA) {
    data.frame(
      label      = label,
      est        = est,
      se         = se,
      t_or_z     = t_val,
      pval       = p_val,
      ci_low     = cilo,
      ci_high    = cihi,
      sig_ci     = sig_ci,   # only meaningful for indirect (bootstrap CI)
      x_var      = x_var,
      y_var      = y_var,
      n_used     = n_used,
      n_boot_fail = n_boot_fail,
      se_source  = ifelse(label == "indirect", "bootstrap", "OLS"),
      ci_source  = ifelse(label == "indirect", "percentile_bootstrap",
                          "t_distribution"),
      stringsAsFactors = FALSE
    )
  }

  # ------------------------------------------------------------------
  # [Fix 3] Decomposition sanity check: c must equal c' + a*b exactly
  # under OLS (complete-cases must be identical across all three models).
  # A difference > 1e-6 signals that complete.cases diverged between models.
  # ------------------------------------------------------------------
  decomp_diff <- abs(c_est - (cp_est + ab_est))
  if (!is.na(decomp_diff) && decomp_diff > 1e-6)
    warning(sprintf(
      "[Decomp check FAILED] %s: |c - (c'+ab)| = %.2e  (c=%.6f, c'=%.6f, ab=%.6f)",
      x_var, decomp_diff, c_est, cp_est, ab_est
    ))

  rbind(
    make_path("a",        a_est,  a_se,       a_t,  a_p,  a_cilo,  a_cihi),
    make_path("b",        b_est,  b_se,       b_t,  b_p,  b_cilo,  b_cihi),
    make_path("cp",       cp_est, cp_se,      cp_t, cp_p, cp_cilo, cp_cihi),
    make_path("total",    c_est,  c_se,       c_t,  c_p,  c_cilo,  c_cihi),
    make_path("indirect", ab_est, ab_boot_se, NA,   NA,   ab_cilo, ab_cihi,
              sig_ci = ab_sig)
  )
}

# ---- Section 5: Run female subgroup -----------------------------------------

cat("\n========== Subgroup: Female (sex_F = 1) ==========\n")
df_female <- df[df$sex_F == 1, ]
# z-scored columns (_c, fso_com_z) already computed on the full sample (Section 3).

res_female_list <- lapply(fc_c_cols, function(xc) {
  cat("  Running FC:", xc, "...\n")
  run_simple_med(data       = df_female,
                 x_var      = xc,
                 y_var      = "fso_com_z",
                 covariates = c("age_c", "meanFD_c"),
                 R          = 5000)
})

results_female          <- do.call(rbind, Filter(Negate(is.null), res_female_list))
results_female$subgroup <- "Female"

# ---- Section 5b: BH-FDR correction -----------------------------------------
# [Fix 1] b-path is now included in FDR correction.
#   Each FC produces its own lm_y (Y ~ M + X + covs), so b varies across
#   the 8 FCs and constitutes 8 independent tests — the same as a, cp, total.
# - a, b, cp, total : 8 independent tests → BH across 8 FCs
# - indirect        : significance by bootstrap CI, no p-value → no FDR

results_female$p_fdr <- NA_real_
for (lbl in c("a", "b", "cp", "total")) {   # [Fix 1] "b" added
  idx <- which(results_female$label == lbl)
  if (length(idx) > 0)
    results_female$p_fdr[idx] <- p.adjust(results_female$pval[idx],
                                           method = "BH")
}

# ---- Section 6: Console summary ---------------------------------------------

print_indirect_summary <- function(tbl, grp_label) {
  cat(sprintf("\n--- Indirect effects: %s ---\n", grp_label))
  indir <- tbl[tbl$label == "indirect", ]
  indir <- indir[order(-abs(indir$est)), ]
  cat(sprintf("  %-24s  %7s  %9s  %9s  sig\n",
              "FC", "a×b", "CI_lo", "CI_hi"))
  for (i in seq_len(nrow(indir))) {
    r <- indir[i, ]
    cat(sprintf("  %-24s  %7.4f  %9.4f  %9.4f  %s\n",
                r$x_var, r$est, r$ci_low, r$ci_high,
                ifelse(isTRUE(r$sig_ci), "*", "")))
  }
}

print_indirect_summary(results_female, "Females (sex_F = 1)")

# ---- Section 7: Export raw results ------------------------------------------

label_order <- c("indirect", "total", "a", "b", "cp")

sort_results <- function(tbl)
  tbl[order(tbl$x_var, match(tbl$label, label_order)), ]

write_xlsx(sort_results(results_female),
           file.path(out_dir, "sim_med_female_v4.xlsx"))
cat("\nRaw results saved: sim_med_female_v4.xlsx\n")

# ---- Section 8: Forest plots ------------------------------------------------

fc_label_map <- c(
  fc_thR_aSTGr_c = "Tha R-aSTG R", fc_thR_aSTGl_c = "Tha R-aSTG L",
  fc_thR_pSTGr_c = "Tha R-pSTG R", fc_thR_pSTGl_c = "Tha R-pSTG L",
  fc_thL_aSTGr_c = "Tha L-aSTG R", fc_thL_aSTGl_c = "Tha L-aSTG L",
  fc_thL_pSTGr_c = "Tha L-pSTG R", fc_thL_pSTGl_c = "Tha L-pSTG L"
)

fdr_stars <- function(p) {
  dplyr::case_when(
    is.na(p) | is.nan(p) | is.infinite(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

raw_t_stars <- function(p) {
  dplyr::case_when(
    is.na(p) | is.nan(p) | is.infinite(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# 8a: Indirect effect (bootstrap percentile CI) -------------------------------
plot_df_indir <- results_female[results_female$label == "indirect", ] |>
  mutate(
    fc_label  = factor(fc_label_map[x_var], levels = rev(fc_label_map)),
    pt_fill   = ifelse(isTRUE(sig_ci), "#E69F00", "white"),
    degen_ci  = !is.na(ci_low) & !is.na(ci_high) &
                  abs(ci_high - ci_low) < 1e-8
  )

forest_indir <- ggplot(plot_df_indir, aes(x = est, y = fc_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 colour = "#E69F00", height = 0.25, linewidth = 0.8,
                 na.rm = TRUE) +
  geom_point(aes(fill = pt_fill),
             colour = "#E69F00", shape = 22L, size = 3.5, stroke = 0.9) +
  scale_fill_identity() +
  labs(
    title    = "Simple mediation — indirect effect (FC → AASP Auditory → SRS): Female",
    subtitle = "Outcome: fso_com; 95% percentile bootstrap CI (R=5000); filled = CI excludes zero",
    x        = "Indirect effect (a×b), 95% percentile bootstrap CI",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "sim_med_forest_female_indirect_v4.png"),
       plot = forest_indir, width = 10, height = 6, dpi = 300)
cat("Indirect effect forest plot saved.\n")

# 8b: Direct effect c′ (OLS CI from lm) --------------------------------------
plot_df_cp <- results_female[results_female$label == "cp", ] |>
  mutate(
    fc_label  = factor(fc_label_map[x_var], levels = rev(fc_label_map)),
    sig_point = !is.na(p_fdr) & p_fdr < 0.05,
    pt_fill   = ifelse(sig_point, "#E69F00", "white"),
    stars     = fdr_stars(p_fdr)
  )

forest_cp <- ggplot(plot_df_cp, aes(x = est, y = fc_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 colour = "#E69F00", height = 0.25, linewidth = 0.8) +
  geom_point(aes(fill = pt_fill),
             colour = "#E69F00", shape = 22L, size = 3.5, stroke = 0.9) +
  geom_text(data = ~ subset(.x, stars != ""),
            aes(x = est, label = stars),
            colour = "#E69F00", size = 5, hjust = 0.5, vjust = -0.7) +
  scale_fill_identity() +
  labs(
    title    = "Simple mediation — direct effect c′ (FC → SRS | Auditory): Female",
    subtitle = "Outcome: fso_com; 95% OLS t-CI; BH-FDR: * p<.05  ** p<.01  *** p<.001",
    x        = "Direct effect (c′), 95% t-distribution CI",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "sim_med_forest_female_direct_v4.png"),
       plot = forest_cp, width = 10, height = 6, dpi = 300)
cat("Direct effect (c') forest plot saved.\n")

# 8c: Total effect (OLS CI from lm) ------------------------------------------
plot_df_total <- results_female[results_female$label == "total", ] |>
  mutate(
    fc_label  = factor(fc_label_map[x_var], levels = rev(fc_label_map)),
    sig_point = !is.na(p_fdr) & p_fdr < 0.05,
    pt_fill   = ifelse(sig_point, "#E69F00", "white"),
    stars     = fdr_stars(p_fdr)
  )

forest_total <- ggplot(plot_df_total, aes(x = est, y = fc_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 colour = "#E69F00", height = 0.25, linewidth = 0.8) +
  geom_point(aes(fill = pt_fill),
             colour = "#E69F00", shape = 22L, size = 3.5, stroke = 0.9) +
  geom_text(data = ~ subset(.x, stars != ""),
            aes(x = est, label = stars),
            colour = "#E69F00", size = 5, hjust = 0.5, vjust = -0.7) +
  scale_fill_identity() +
  labs(
    title    = "Simple mediation — total effect c (FC → SRS): Female",
    subtitle = "Outcome: fso_com; 95% OLS t-CI; BH-FDR: * p<.05  ** p<.01  *** p<.001",
    x        = "Total effect (c), 95% t-distribution CI",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "sim_med_forest_female_total_v4.png"),
       plot = forest_total, width = 10, height = 6, dpi = 300)
cat("Total effect forest plot saved.\n")

# ---- Section 9: Summary Excel -----------------------------------------------
# Format:
#   a, b, cp, total : "est ± SE*"  SE = OLS SE, * = BH-FDR from t-test p
#   indirect        : "est [CI_lo, CI_hi]*"  CI = percentile bootstrap

fem_res          <- results_female
fem_res$FC_label <- unname(fc_label_map[fem_res$x_var])

reg_fem      <- fem_res[fem_res$label %in% c("a", "b", "cp", "total"), ]
# [Fix 1] All four paths (a, b, cp, total) now use BH-FDR corrected stars.
reg_fem$cell <- sprintf("%.3f ± %.3f%s", reg_fem$est, reg_fem$se,
                        fdr_stars(reg_fem$p_fdr))

indir_fem       <- fem_res[fem_res$label == "indirect", ]
indir_fem$cell  <- sprintf("%.3f [%.3f, %.3f]%s",
                            indir_fem$est,
                            indir_fem$ci_low,
                            indir_fem$ci_high,
                            ifelse(!is.na(indir_fem$sig_ci) & indir_fem$sig_ci, "*", ""))
indir_fem$label <- "ab"

param_order_fem <- c("a", "b", "cp", "total", "ab")
fc_label_order  <- unname(fc_label_map)

female_summary <- do.call(rbind, list(
  reg_fem  [, c("FC_label", "label", "cell")],
  indir_fem[, c("FC_label", "label", "cell")]
)) |>
  mutate(
    FC_label = factor(FC_label, levels = fc_label_order),
    label    = factor(label,    levels = param_order_fem)
  ) |>
  arrange(FC_label, label) |>
  tidyr::pivot_wider(names_from = label, values_from = cell) |>
  rename(FC = FC_label) |>
  mutate(FC = as.character(FC))

legend_fem <- data.frame(
  Item = c(
    "Format (a, b, cp, total)",
    "Format (ab / indirect)",
    "SE source (a, b, cp, total)",
    "SE source (ab)",
    "CI source (a, b, cp, total)",
    "CI source (ab)",
    "Stars (a, b, cp, total)",
    "Stars (ab)",
    "[v4 Fix 1] b-path FDR",
    "[v4 Fix 2] lavaan fixed.x",
    "[v4 Fix 3] Decomp check",
    "[v4 Fix 5] Standardization"
  ),
  Description = c(
    "estimate ± SE; stars from BH-FDR corrected t-test p-value",
    "estimate [95% CI_lo, CI_hi]; * if bootstrap CI excludes zero",
    "OLS SE from lm() — identical to PROCESS macro output",
    "Bootstrap SD from lavaan (se='bootstrap', R=5000) — reported for reference only",
    "95% t-distribution CI from confint(lm())",
    "95% percentile bootstrap CI from lavaan parameterestimates(boot.ci.type='perc')",
    "* p_fdr<.05  ** p_fdr<.01  *** p_fdr<.001 — BH-FDR across 8 FCs within each label",
    "* = 95% percentile bootstrap CI excludes zero",
    "b-path now included in BH-FDR (8 tests, one per FC); previously uncorrected",
    "fixed.x = TRUE: covariates treated as fixed exogenous (matches PROCESS; reduces bootstrap failures)",
    "c = c' + a*b sanity check; warnings printed to console if |diff| > 1e-6",
    "Z-scores computed on the full sample before sex subsetting (whole-sample SD/mean) — intentional for cross-sex comparability"
  ),
  stringsAsFactors = FALSE
)

writexl::write_xlsx(
  list(
    "Female_Summary" = as.data.frame(female_summary),
    "Legend"         = legend_fem
  ),
  file.path(out_dir, "sim_AuM_female_summary_v4.xlsx")
)
cat("Female summary Excel saved: sim_AuM_female_summary_v4.xlsx\n")

cat("\n====== Done. All outputs written to:\n", out_dir, "\n")
cat("v4 key changes (relative to v3):\n")
cat("  [Fix 1] BH-FDR now applied to b-path (8 tests across FCs).\n")
cat("  [Fix 2] lavaan sem() uses fixed.x = TRUE (matches PROCESS macro).\n")
cat("  [Fix 3] Decomposition sanity check: warns if |c - (c'+ab)| > 1e-6.\n")
cat("  [Fix 5] Z-score standardization: whole-sample (reverted to original design).\n")
sessionInfo()
