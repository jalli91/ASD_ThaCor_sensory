# =============================================================================
# Linear Regression: Thalamo-Auditory FC → AASP Auditory score
#                    Thalamo-Auditory FC → SRS fso_com
# Subgroups: full | males | females | ASD | NAC | ASD×sex | NAC×sex
# Covariates: age_c + meanFD_c + sex_F (mixed-sex groups)
#             age_c + meanFD_c only    (within-sex groups)
# Forest plots styled to match simple_mediation_AuM.R
# =============================================================================

# ---- Section 1: Setup -------------------------------------------------------

pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "writexl")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

set.seed(42)

data_file <- "/Users/james/Desktop/CYL/ASD_ThaCor/Score_sheet/drop_nosig_raw/Combine_All_cleaned_0.01-0.08_drop_nosig_raw_AuTha.xlsx"
out_dir   <- "/Users/james/Desktop/CYL/ASD_ThaCor/Score_sheet/drop_nosig_raw/ForGraph_lm_FC"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
outcomes  <- c("auditory_c", "fso_com_z")

# ---- Section 3: Linear regression function ----------------------------------
#
# Runs lm(y ~ FC + covariates) and returns the FC coefficient row.
# FC (X) is raw.
# Within-sex subgroups omit sex_F (constant within the stratum).

run_lm <- function(data, x_var, y_var,
                   covariates = c("age_c", "meanFD_c", "sex_F")) {
  all_vars <- c(y_var, x_var, covariates)
  d <- data[complete.cases(data[, all_vars, drop = FALSE]), ]
  n_used <- nrow(d)

  if (n_used < length(all_vars) + 2) {
    message("  Skipping ", x_var, " -> ", y_var, ": n too small (", n_used, ")")
    return(NULL)
  }

  rhs <- paste(c(x_var, covariates), collapse = " + ")
  fit <- tryCatch(
    lm(as.formula(paste0(y_var, " ~ ", rhs)), data = d),
    error = function(e) { message("  lm error: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)

  s  <- summary(fit)
  ci <- confint(fit, level = 0.95)

  if (!(x_var %in% rownames(s$coefficients))) return(NULL)

  data.frame(
    x_var   = x_var,
    y_var   = y_var,
    est     = s$coefficients[x_var, "Estimate"],
    se      = s$coefficients[x_var, "Std. Error"],
    t_val   = s$coefficients[x_var, "t value"],
    pval    = s$coefficients[x_var, "Pr(>|t|)"],
    ci_low  = ci[x_var, 1],
    ci_high = ci[x_var, 2],
    n_used  = n_used,
    stringsAsFactors = FALSE
  )
}

run_subgroup <- function(data, covariates = c("age_c", "meanFD_c", "sex_F")) {
  rows <- list()
  for (xc in fc_c_cols)
    for (yv in outcomes)
      rows[[length(rows) + 1]] <- run_lm(data, xc, yv, covariates)
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ---- Section 4: Run all subgroups -------------------------------------------

cat("\n===== Full =====\n")
results_full          <- run_subgroup(df)
results_full$subgroup <- "Full"

cat("\n===== Male (sex_F = 0) =====\n")
df_male <- df[df$sex_F == 0, ]
results_male          <- run_subgroup(df_male, c("age_c", "meanFD_c"))
results_male$subgroup <- "Male"

cat("\n===== Female (sex_F = 1) =====\n")
df_female <- df[df$sex_F == 1, ]
results_female          <- run_subgroup(df_female, c("age_c", "meanFD_c"))
results_female$subgroup <- "Female"

cat("\n===== ASD =====\n")
df_asd <- df[!is.na(df$group) & df$group == "ASD", ]
cat("  n =", nrow(df_asd), "\n")
results_asd          <- run_subgroup(df_asd)
results_asd$subgroup <- "ASD"

cat("\n===== NAC =====\n")
df_nac <- df[!is.na(df$group) & df$group == "NAC", ]
cat("  n =", nrow(df_nac), "\n")
results_nac          <- run_subgroup(df_nac)
results_nac$subgroup <- "NAC"

cat("\n===== ASD Male =====\n")
df_asd_male <- df[!is.na(df$group) & df$group == "ASD" & df$sex_F == 0, ]
cat("  n =", nrow(df_asd_male), "\n")
results_asd_male          <- run_subgroup(df_asd_male, c("age_c", "meanFD_c"))
results_asd_male$subgroup <- "ASD_Male"

cat("\n===== ASD Female =====\n")
df_asd_female <- df[!is.na(df$group) & df$group == "ASD" & df$sex_F == 1, ]
cat("  n =", nrow(df_asd_female), "\n")
results_asd_female          <- run_subgroup(df_asd_female, c("age_c", "meanFD_c"))
results_asd_female$subgroup <- "ASD_Female"

cat("\n===== NAC Male =====\n")
df_nac_male <- df[!is.na(df$group) & df$group == "NAC" & df$sex_F == 0, ]
cat("  n =", nrow(df_nac_male), "\n")
results_nac_male          <- run_subgroup(df_nac_male, c("age_c", "meanFD_c"))
results_nac_male$subgroup <- "NAC_Male"

cat("\n===== NAC Female =====\n")
df_nac_female <- df[!is.na(df$group) & df$group == "NAC" & df$sex_F == 1, ]
cat("  n =", nrow(df_nac_female), "\n")
results_nac_female          <- run_subgroup(df_nac_female, c("age_c", "meanFD_c"))
results_nac_female$subgroup <- "NAC_Female"

# ---- Section 5: FDR correction (BH) ----------------------------------------
# Applied separately per outcome across the 8 FC tests within each subgroup.

apply_fdr <- function(tbl) {
  tbl$pFDR   <- NA_real_
  tbl$sig_fdr <- NA_character_
  for (yv in outcomes) {
    idx <- which(tbl$y_var == yv)
    if (length(idx) > 0) {
      tbl$pFDR[idx]   <- p.adjust(tbl$pval[idx], method = "BH")
      tbl$sig_fdr[idx] <- ifelse(tbl$pFDR[idx] < 0.05, "*", "")
    }
  }
  tbl
}

results_full       <- apply_fdr(results_full)
results_male       <- apply_fdr(results_male)
results_female     <- apply_fdr(results_female)
results_asd        <- apply_fdr(results_asd)
results_nac        <- apply_fdr(results_nac)
results_asd_male   <- apply_fdr(results_asd_male)
results_asd_female <- apply_fdr(results_asd_female)
results_nac_male   <- apply_fdr(results_nac_male)
results_nac_female <- apply_fdr(results_nac_female)

sig_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ""))))
}

# ---- Section 6: Console summary --------------------------------------------

fc_label_map <- c(
  fc_thR_aSTGr_c = "Thalamus R-aSTG R FC", fc_thR_aSTGl_c = "Thalamus R-aSTG L FC",
  fc_thR_pSTGr_c = "Thalamus R-pSTG R FC", fc_thR_pSTGl_c = "Thalamus R-pSTG L FC",
  fc_thL_aSTGr_c = "Thalamus L-aSTG R FC", fc_thL_aSTGl_c = "Thalamus L-aSTG L FC",
  fc_thL_pSTGr_c = "Thalamus L-pSTG R FC", fc_thL_pSTGl_c = "Thalamus L-pSTG L FC"
)

print_lm_summary <- function(tbl, grp_label) {
  cat(sprintf("\n--- %s ---\n", grp_label))
  for (yv in outcomes) {
    rows <- tbl[tbl$y_var == yv, ]
    rows <- rows[order(rows$pFDR), ]
    cat(sprintf("  Outcome: %s\n", yv))
    cat(sprintf("  %-18s  %7s  %7s  %9s  %9s  %7s  %7s  sig\n",
                "FC", "est", "SE", "CI_lo", "CI_hi", "p", "pFDR"))
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      cat(sprintf("  %-18s  %7.4f  %7.4f  %9.4f  %9.4f  %7.4f  %7s   %s\n",
                  fc_label_map[r$x_var], r$est, r$se,
                  r$ci_low, r$ci_high, r$pval,
                  ifelse(is.na(r$pFDR), "  NA  ", sprintf("%7.4f", r$pFDR)),
                  ifelse(is.na(r$sig_fdr), "", r$sig_fdr)))
    }
  }
}

grp_list <- list(
  list(results_full,       "Full sample"),
  list(results_male,       "Males (sex_F = 0)"),
  list(results_female,     "Females (sex_F = 1)"),
  list(results_asd,        "ASD"),
  list(results_nac,        "NAC"),
  list(results_asd_male,   "ASD Male"),
  list(results_asd_female, "ASD Female"),
  list(results_nac_male,   "NAC Male"),
  list(results_nac_female, "NAC Female")
)

for (g in grp_list) print_lm_summary(g[[1]], g[[2]])

# ---- Section 7: Export results ---------------------------------------------

all_results <- do.call(rbind, lapply(grp_list, `[[`, 1))
all_results <- all_results[order(all_results$subgroup, all_results$y_var,
                                 all_results$x_var), ]

write_xlsx(
  list(
    all        = all_results,
    full       = results_full,
    male       = results_male,
    female     = results_female,
    asd        = results_asd,
    nac        = results_nac,
    asd_male   = results_asd_male,
    asd_female = results_asd_female,
    nac_male   = results_nac_male,
    nac_female = results_nac_female
  ),
  file.path(out_dir, "lm_regression_results.xlsx")
)
cat("\nResults saved: lm_regression_results.xlsx\n")

# ---- Section 8: Forest plots -----------------------------------------------

outcome_labels <- c(
  auditory_c = "AASP Auditory score",
  fso_com_z  = "SRS fso_com"
)

# Wong (2011) colorblind-safe palette — same as simple_mediation_AuM.R
pal_sex  <- c("Full" = "#0072B2", "Male" = "#009E73", "Female" = "#D55E00")
shp_sex  <- c("Full" = 21L, "Male" = 24L, "Female" = 22L)
pal_diag <- c("ASD" = "#CC79A7", "NAC" = "#56B4E9")
shp_diag <- c("ASD" = 21L, "NAC" = 22L)
pal_asd_sex <- c("ASD_Male" = "#009E73", "ASD_Female" = "#D55E00")
shp_asd_sex <- c("ASD_Male" = 24L,       "ASD_Female" = 22L)
pal_nac_sex <- c("NAC_Male" = "#56B4E9", "NAC_Female" = "#E69F00")
shp_nac_sex <- c("NAC_Male" = 24L,       "NAC_Female" = 22L)

make_forest <- function(results_list, subgroup_levels, outcome,
                        pal_col, pal_shp, title_str) {
  plot_df <- do.call(rbind, lapply(results_list, function(tbl) {
    tbl[tbl$y_var == outcome, ]
  })) |>
    mutate(
      fc_label  = factor(fc_label_map[x_var], levels = rev(fc_label_map)),
      subgroup  = factor(subgroup, levels = subgroup_levels),
      sig_point = !is.na(pFDR) & pFDR < 0.05,
      pt_fill   = ifelse(sig_point, pal_col[as.character(subgroup)], "white"),
      degen_ci  = abs(ci_high - ci_low) < 1e-8,
      stars     = sig_stars(pFDR)
    )

  ggplot(plot_df, aes(x = est, y = fc_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(
      aes(xmin = ci_low, xmax = ci_high, colour = subgroup),
      height   = 0.25, linewidth = 0.65,
      position = position_dodge(width = 0.6)
    ) +
    geom_point(
      aes(colour = subgroup, fill = pt_fill, shape = subgroup),
      size     = 3.2, stroke = 0.9,
      position = position_dodge(width = 0.6)
    ) +
    geom_text(
      data = ~ subset(.x, stars != ""),
      aes(x = est, label = stars, colour = subgroup),
      size = 5, hjust = 0.5, vjust = -0.7,
      position = position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    geom_text(
      data = ~ subset(.x, degen_ci),
      aes(x = est, label = paste0("n=", n_used, "\nCI n/a"), colour = subgroup),
      size = 2.5, hjust = -0.05,
      position = position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    scale_colour_manual(name = "Subgroup", values = pal_col) +
    scale_fill_identity() +
    scale_shape_manual(name = "Subgroup", values = pal_shp) +
    guides(
      colour = guide_legend(override.aes = list(
        shape = unname(pal_shp), fill = unname(pal_col), size = 3.5
      )),
      shape = "none", fill = "none"
    ) +
    labs(
      title    = title_str,
      subtitle = paste0("Outcome: ", outcome_labels[outcome],
                        "; 95% CI (t-dist); BH-FDR: * p<.05  ** p<.01  *** p<.001"),
      x        = paste0("β (standardised), 95% CI"),
      y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

make_forest_4grp <- function(results_list, outcome, title_str) {
  plot_df <- do.call(rbind, lapply(results_list, function(tbl) {
    tbl[tbl$y_var == outcome, ]
  })) |>
    mutate(
      fc_label  = factor(fc_label_map[x_var], levels = rev(fc_label_map)),
      diagnosis = factor(ifelse(grepl("^ASD", subgroup), "ASD", "NAC"),
                         levels = c("ASD", "NAC")),
      sex       = factor(ifelse(grepl("Male$", subgroup), "Male", "Female"),
                         levels = c("Male", "Female")),
      sig_point = !is.na(pFDR) & pFDR < 0.05,
      pt_fill   = ifelse(sig_point,
                         c("Male" = "#009E73", "Female" = "#D55E00")[
                           as.character(ifelse(grepl("Male$", subgroup), "Male", "Female"))],
                         "white"),
      degen_ci  = abs(ci_high - ci_low) < 1e-8,
      stars     = sig_stars(pFDR)
    )

  pal_s <- c("Male" = "#009E73", "Female" = "#D55E00")
  shp_s <- c("Male" = 24L, "Female" = 22L)

  ggplot(plot_df, aes(x = est, y = fc_label)) +
    facet_wrap(~ diagnosis, ncol = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(
      aes(xmin = ci_low, xmax = ci_high, colour = sex),
      height   = 0.25, linewidth = 0.65,
      position = position_dodge(width = 0.5)
    ) +
    geom_point(
      aes(colour = sex, fill = pt_fill, shape = sex),
      size     = 3.2, stroke = 0.9,
      position = position_dodge(width = 0.5)
    ) +
    geom_text(
      data = ~ subset(.x, stars != ""),
      aes(x = est, label = stars, colour = sex),
      size = 5, hjust = 0.5, vjust = -0.7,
      position = position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    geom_text(
      data = ~ subset(.x, degen_ci),
      aes(x = est, label = paste0("n=", n_used, "\nCI n/a"), colour = sex),
      size = 2.5, hjust = -0.05,
      position = position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    scale_colour_manual(name = "Sex", values = pal_s) +
    scale_fill_identity() +
    scale_shape_manual(name = "Sex", values = shp_s) +
    guides(
      colour = guide_legend(override.aes = list(
        shape = unname(shp_s), fill = unname(pal_s), size = 3.5
      )),
      shape = "none", fill = "none"
    ) +
    labs(
      title    = title_str,
      subtitle = paste0("Outcome: ", outcome_labels[outcome],
                        "; 95% CI (t-dist); BH-FDR: * p<.05  ** p<.01  *** p<.001"),
      x        = "β (standardised), 95% CI",
      y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = 12, face = "bold")
    )
}

# --- Produce all forest plots for each outcome --------------------------------

for (yv in outcomes) {

  ylab <- gsub("[^A-Za-z0-9]", "_", outcome_labels[yv])

  # 8a: Sex subgroups (Full / Male / Female)
  ggsave(
    file.path(out_dir, paste0("lm_forest_sex_", ylab, ".png")),
    plot = make_forest(
      results_list     = list(results_full, results_male, results_female),
      subgroup_levels  = c("Full", "Male", "Female"),
      outcome          = yv,
      pal_col          = pal_sex,
      pal_shp          = shp_sex,
      title_str        = paste0("Linear regression (FC → ", outcome_labels[yv],
                                "): Full / Male / Female")
    ),
    width = 14, height = 7, dpi = 300
  )
  cat("Saved: lm_forest_sex_", ylab, ".png\n", sep = "")

  # 8b: Diagnostic subgroups (ASD / NAC)
  ggsave(
    file.path(out_dir, paste0("lm_forest_diag_", ylab, ".png")),
    plot = make_forest(
      results_list     = list(results_asd, results_nac),
      subgroup_levels  = c("ASD", "NAC"),
      outcome          = yv,
      pal_col          = pal_diag,
      pal_shp          = shp_diag,
      title_str        = paste0("Linear regression (FC → ", outcome_labels[yv],
                                "): ASD vs NAC")
    ),
    width = 14, height = 7, dpi = 300
  )
  cat("Saved: lm_forest_diag_", ylab, ".png\n", sep = "")

  # 8c: Within-ASD sex (ASD_Male / ASD_Female)
  ggsave(
    file.path(out_dir, paste0("lm_forest_asd_sex_", ylab, ".png")),
    plot = make_forest(
      results_list     = list(results_asd_male, results_asd_female),
      subgroup_levels  = c("ASD_Male", "ASD_Female"),
      outcome          = yv,
      pal_col          = pal_asd_sex,
      pal_shp          = shp_asd_sex,
      title_str        = paste0("Linear regression (FC → ", outcome_labels[yv],
                                "): ASD by sex")
    ),
    width = 14, height = 7, dpi = 300
  )
  cat("Saved: lm_forest_asd_sex_", ylab, ".png\n", sep = "")

  # 8d: Within-NAC sex (NAC_Male / NAC_Female)
  ggsave(
    file.path(out_dir, paste0("lm_forest_nac_sex_", ylab, ".png")),
    plot = make_forest(
      results_list     = list(results_nac_male, results_nac_female),
      subgroup_levels  = c("NAC_Male", "NAC_Female"),
      outcome          = yv,
      pal_col          = pal_nac_sex,
      pal_shp          = shp_nac_sex,
      title_str        = paste0("Linear regression (FC → ", outcome_labels[yv],
                                "): NAC by sex")
    ),
    width = 14, height = 7, dpi = 300
  )
  cat("Saved: lm_forest_nac_sex_", ylab, ".png\n", sep = "")

  # 8e: 4-group faceted (ASD/NAC × Male/Female)
  ggsave(
    file.path(out_dir, paste0("lm_forest_4group_", ylab, ".png")),
    plot = make_forest_4grp(
      results_list = list(results_asd_male, results_asd_female,
                          results_nac_male, results_nac_female),
      outcome      = yv,
      title_str    = paste0("Linear regression (FC → ", outcome_labels[yv],
                            "): diagnosis × sex")
    ),
    width = 18, height = 7, dpi = 300
  )
  cat("Saved: lm_forest_4group_", ylab, ".png\n", sep = "")
}

scatter_y_labels <- c(
  auditory_c = "Auditory subscore",
  fso_com_z  = "Social communication subscore"
)

# ---- Section 9: Scatter plots (FC → outcome): Full / Male / Female overlay --
# Partial regression plots: both FC and outcome are residualised against the
# same covariates used in run_lm(), so the plotted slope matches the forest β.

# Helper: residualise FC cols and outcome against covariates; return long format.
partial_long <- function(data, yv, covariates, grp_var, grp_label,
                         extra_cols = character(0)) {
  req <- c(fc_c_cols, yv, covariates)
  d   <- data[complete.cases(data[, req, drop = FALSE]), ]
  if (nrow(d) == 0) return(NULL)
  cov_rhs <- paste(covariates, collapse = " + ")
  y_res   <- residuals(lm(as.formula(paste(yv, "~", cov_rhs)), data = d))
  rows <- lapply(fc_c_cols, function(fc) {
    fc_res <- residuals(lm(as.formula(paste(fc, "~", cov_rhs)), data = d))
    base   <- data.frame(fc_var  = fc,
                         fc_val  = as.numeric(fc_res),
                         y_resid = as.numeric(y_res),
                         stringsAsFactors = FALSE)
    for (ec in extra_cols) base[[ec]] <- d[[ec]]
    base
  })
  out <- do.call(rbind, rows)
  out[[grp_var]] <- grp_label
  out$fc_label   <- factor(fc_label_map[out$fc_var], levels = fc_label_map)
  out
}

pal_sex_scatter <- c("Full" = "#0072B2", "Male" = "#009E73", "Female" = "#D55E00")
shp_sex_scatter <- c("Full" = 21L,       "Male" = 24L,       "Female" = 22L)

for (yv in outcomes) {

  scatter_long <- bind_rows(
    partial_long(df,        yv, c("age_c", "meanFD_c", "sex_F"), "sex_grp", "Full"),
    partial_long(df_male,   yv, c("age_c", "meanFD_c"),          "sex_grp", "Male"),
    partial_long(df_female, yv, c("age_c", "meanFD_c"),          "sex_grp", "Female")
  ) |>
    mutate(sex_grp = factor(sex_grp, levels = c("Full", "Male", "Female")))

  annot_sex_all <- do.call(rbind, lapply(
    list(list(tbl = results_full,   grp = "Full"),
         list(tbl = results_male,   grp = "Male"),
         list(tbl = results_female, grp = "Female")),
    function(x) {
      sub <- x$tbl[x$tbl$y_var == yv, ]
      data.frame(
        fc_label  = factor(fc_label_map[sub$x_var], levels = fc_label_map),
        sex_grp   = factor(x$grp, levels = c("Full", "Male", "Female")),
        stars     = sig_stars(sub$pFDR),
        ann_label = sprintf("β=%.3f  pFDR=%.3f", sub$est, sub$pFDR),
        stringsAsFactors = FALSE
      )
    }
  ))
  annot_sex <- annot_sex_all[annot_sex_all$stars != "", ]

  p <- ggplot(scatter_long,
              aes(x = fc_val, y = y_resid,
                  colour = sex_grp, fill = sex_grp)) +
    facet_grid(sex_grp ~ fc_label, scales = "free") +
    geom_point(alpha = 0.5, size = 1.8, stroke = 0.6) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    geom_text(
      data        = annot_sex,
      aes(x = Inf, y = Inf, label = stars, colour = sex_grp),
      vjust = 1.3, hjust = 1.05,
      size = 5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = annot_sex_all,
      aes(x = -Inf, y = Inf, label = ann_label, colour = sex_grp),
      vjust = 1.3, hjust = -0.05,
      size = 3.2,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(name = "Subgroup", values = pal_sex_scatter) +
    scale_fill_manual(name = "Subgroup",   values = pal_sex_scatter) +
    guides(
      colour = guide_legend(override.aes = list(
        size = 3, linetype = 1, linewidth = 0.9
      )),
      fill = "none"
    ) +
    labs(
      title    = sprintf("FC → %s: Full / Male / Female (partial regression)",
                         outcome_labels[yv]),
      subtitle = "Each row: within-group covariate-adjusted residuals (age, meanFD [+ sex for Full]); slopes match forest plot β;  * p<.05   ** p<.01   *** p<.001 (BH-FDR)",
      x        = "FC (residualised)",
      y        = paste0(scatter_y_labels[yv], " (residualised)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x     = element_text(size = 8, face = "bold"),
      strip.text.y     = element_text(size = 10, face = "bold")
    )

  fname <- file.path(out_dir,
                     paste0("scatter_sex_",
                            gsub("[^A-Za-z0-9]", "_", outcome_labels[yv]),
                            ".png"))
  ggsave(fname, plot = p, width = 18, height = 12, dpi = 300)
  cat("Saved:", basename(fname), "\n")
}

# ---- Section 10: Scatter plots (FC → outcome): Full / ASD / NAC overlay ----

pal_diag_scatter <- c("Full" = "#2CA02C", "ASD" = "#D62728", "NAC" = "#1F77B4")
shp_diag_scatter <- c("Full" = 21L,       "ASD" = 22L,       "NAC" = 24L)

for (yv in outcomes) {

  scatter_long_diag <- bind_rows(
    partial_long(df,     yv, c("age_c", "meanFD_c", "sex_F"), "diag_grp", "Full", extra_cols = "group"),
    partial_long(df_asd, yv, c("age_c", "meanFD_c", "sex_F"), "diag_grp", "ASD",  extra_cols = "group"),
    partial_long(df_nac, yv, c("age_c", "meanFD_c", "sex_F"), "diag_grp", "NAC",  extra_cols = "group")
  ) |>
    mutate(diag_grp = factor(diag_grp, levels = c("Full", "ASD", "NAC")),
           group    = factor(group,    levels = c("ASD", "NAC")))

  annot_diag_all <- do.call(rbind, lapply(
    list(list(tbl = results_full, grp = "Full"),
         list(tbl = results_asd,  grp = "ASD"),
         list(tbl = results_nac,  grp = "NAC")),
    function(x) {
      sub <- x$tbl[x$tbl$y_var == yv, ]
      data.frame(
        fc_label  = factor(fc_label_map[sub$x_var], levels = fc_label_map),
        diag_grp  = factor(x$grp, levels = c("Full", "ASD", "NAC")),
        stars     = sig_stars(sub$pFDR),
        ann_label = sprintf("β=%.3f  pFDR=%.3f", sub$est, sub$pFDR),
        stringsAsFactors = FALSE
      )
    }
  ))
  annot_diag <- annot_diag_all[annot_diag_all$stars != "", ]

  p <- ggplot(scatter_long_diag, aes(x = fc_val, y = y_resid)) +
    facet_grid(diag_grp ~ fc_label, scales = "free") +
    geom_point(aes(colour = group), alpha = 0.5, size = 1.8, stroke = 0.6) +
    geom_smooth(aes(colour = diag_grp, fill = diag_grp),
                method = "lm", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    geom_text(
      data        = annot_diag,
      aes(x = Inf, y = Inf, label = stars, colour = diag_grp),
      vjust = 1.3, hjust = 1.05,
      size = 5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = annot_diag_all,
      aes(x = -Inf, y = Inf, label = ann_label, colour = diag_grp),
      vjust = 1.3, hjust = -0.05,
      size = 3.2,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(name = "Diagnosis", values = pal_diag_scatter) +
    scale_fill_manual(name = "Diagnosis",   values = pal_diag_scatter, guide = "none") +
    guides(
      colour = guide_legend(override.aes = list(
        shape = 16, size = 3, linetype = 1, linewidth = 0.9
      ))
    ) +
    labs(
      title    = sprintf("FC → %s: Full / ASD / NAC (partial regression)",
                         outcome_labels[yv]),
      subtitle = "Each row: within-group covariate-adjusted residuals (age, meanFD, sex); slopes match forest plot β;  * p<.05   ** p<.01   *** p<.001 (BH-FDR)",
      x        = "FC (residualised)",
      y        = paste0(scatter_y_labels[yv], " (residualised)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x     = element_text(size = 8, face = "bold"),
      strip.text.y     = element_text(size = 10, face = "bold")
    )

  fname <- file.path(out_dir,
                     paste0("scatter_diag_",
                            gsub("[^A-Za-z0-9]", "_", outcome_labels[yv]),
                            ".png"))
  ggsave(fname, plot = p, width = 18, height = 12, dpi = 300)
  cat("Saved:", basename(fname), "\n")
}

# ---- Section 12: Scatter (FC → outcome): Full/ASD/NAC — male-only & female-only --
# Partial regression version: residualise against age_c + meanFD_c (sex is
# constant within each sex stratum so it is not added as a covariate).
# Stars in the top-right corner come from the covariate-adjusted forest results.

make_diag_sex_scatter <- function(df_full_sex, df_asd_sex, df_nac_sex,
                                  res_full_sex, res_asd_sex, res_nac_sex,
                                  sex_label, yv) {

  cov <- c("age_c", "meanFD_c")

  plot_data <- bind_rows(
    partial_long(df_full_sex, yv, cov, "diag_grp", "Full", extra_cols = "group"),
    partial_long(df_asd_sex,  yv, cov, "diag_grp", "ASD",  extra_cols = "group"),
    partial_long(df_nac_sex,  yv, cov, "diag_grp", "NAC",  extra_cols = "group")
  ) |>
    mutate(diag_grp = factor(diag_grp, levels = c("Full", "ASD", "NAC")),
           group    = factor(group,    levels = c("ASD", "NAC")))

  annot_df2_all <- do.call(rbind, lapply(
    list(list(tbl = res_full_sex, grp = "Full"),
         list(tbl = res_asd_sex,  grp = "ASD"),
         list(tbl = res_nac_sex,  grp = "NAC")),
    function(x) {
      sub <- x$tbl[x$tbl$y_var == yv, ]
      data.frame(
        fc_label  = factor(fc_label_map[sub$x_var], levels = fc_label_map),
        diag_grp  = factor(x$grp, levels = c("Full", "ASD", "NAC")),
        stars     = sig_stars(sub$pFDR),
        ann_label = sprintf("β=%.3f  pFDR=%.3f", sub$est, sub$pFDR),
        stringsAsFactors = FALSE
      )
    }
  ))
  annot_df2 <- annot_df2_all |> filter(stars != "")

  ggplot(plot_data, aes(x = fc_val, y = y_resid)) +
    facet_grid(diag_grp ~ fc_label, scales = "free") +
    geom_point(aes(colour = group), alpha = 0.5, size = 1.8, stroke = 0.6) +
    geom_smooth(aes(colour = diag_grp, fill = diag_grp),
                method = "lm", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    geom_text(
      data        = annot_df2,
      aes(x = Inf, y = Inf, label = stars, colour = diag_grp),
      vjust = 1.3, hjust = 1.05,
      size = 5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = annot_df2_all,
      aes(x = -Inf, y = Inf, label = ann_label, colour = diag_grp),
      vjust = 1.3, hjust = -0.05,
      size = 3.2,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(name = "Diagnosis", values = pal_diag_scatter) +
    scale_fill_manual(name = "Diagnosis",   values = pal_diag_scatter, guide = "none") +
    guides(
      colour = guide_legend(override.aes = list(
        shape = 16, size = 3, linetype = 1, linewidth = 0.9
      ))
    ) +
    labs(
      title    = sprintf("FC → %s: Full / ASD / NAC  [%s only] (partial regression)",
                         outcome_labels[yv], sex_label),
      subtitle = "Each row: within-group covariate-adjusted residuals (age, meanFD);  * p<.05   ** p<.01   *** p<.001 (BH-FDR)",
      x        = "FC (residualised)",
      y        = paste0(scatter_y_labels[yv], " (residualised)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x     = element_text(size = 8, face = "bold"),
      strip.text.y     = element_text(size = 10, face = "bold")
    )
}

for (yv in outcomes) {

  ylab <- gsub("[^A-Za-z0-9]", "_", outcome_labels[yv])

  # Male-only
  p_male <- make_diag_sex_scatter(
    df_male, df_asd_male, df_nac_male,
    results_male, results_asd_male, results_nac_male,
    "Male", yv
  )
  fname_m <- file.path(out_dir, paste0("scatter_diag_male_", ylab, ".png"))
  ggsave(fname_m, plot = p_male, width = 18, height = 12, dpi = 300)
  cat("Saved:", basename(fname_m), "\n")

  # Female-only
  p_female <- make_diag_sex_scatter(
    df_female, df_asd_female, df_nac_female,
    results_female, results_asd_female, results_nac_female,
    "Female", yv
  )
  fname_f <- file.path(out_dir, paste0("scatter_diag_female_", ylab, ".png"))
  ggsave(fname_f, plot = p_female, width = 18, height = 12, dpi = 300)
  cat("Saved:", basename(fname_f), "\n")
}

# ---- Section 13: Export beta / p / pFDR for Section 12 subgroups ----------

export_src <- list(
  list(tbl = results_male,       sex = "Male",   grp = "Full"),
  list(tbl = results_asd_male,   sex = "Male",   grp = "ASD"),
  list(tbl = results_nac_male,   sex = "Male",   grp = "NAC"),
  list(tbl = results_female,     sex = "Female", grp = "Full"),
  list(tbl = results_asd_female, sex = "Female", grp = "ASD"),
  list(tbl = results_nac_female, sex = "Female", grp = "NAC")
)

export_df <- do.call(rbind, lapply(export_src, function(x) {
  x$tbl |>
    mutate(
      sex       = x$sex,
      diagnosis = x$grp,
      fc_label  = fc_label_map[x_var],
      stars_fdr = sig_stars(pFDR)
    ) |>
    select(sex, diagnosis, fc_label, y_var,
           beta = est, SE = se, t_val,
           p_original = pval, pFDR, stars_fdr,
           n_used)
})) |>
  mutate(
    outcome   = outcome_labels[y_var],
    sex       = factor(sex,       levels = c("Male", "Female")),
    diagnosis = factor(diagnosis, levels = c("Full", "ASD", "NAC"))
  ) |>
  arrange(sex, diagnosis, outcome, fc_label) |>
  select(sex, diagnosis, outcome, fc_label,
         beta, SE, t_val, p_original, pFDR, stars_fdr, n_used)

write_xlsx(
  list(
    all        = export_df,
    male_full  = export_df[export_df$sex == "Male"   & export_df$diagnosis == "Full", ],
    male_ASD   = export_df[export_df$sex == "Male"   & export_df$diagnosis == "ASD",  ],
    male_NAC   = export_df[export_df$sex == "Male"   & export_df$diagnosis == "NAC",  ],
    female_full= export_df[export_df$sex == "Female" & export_df$diagnosis == "Full", ],
    female_ASD = export_df[export_df$sex == "Female" & export_df$diagnosis == "ASD",  ],
    female_NAC = export_df[export_df$sex == "Female" & export_df$diagnosis == "NAC",  ]
  ),
  file.path(out_dir, "lm_FC_sex_diag_results.xlsx")
)
cat("Saved: lm_FC_sex_diag_results.xlsx\n")

# ---- Section 14: Forest plots — Male vs Female only (new folder) ------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal_mf <- c("Male" = "#7B79CC", "Female" = "#E69F00")
shp_mf <- c("Male" = 24L,       "Female" = 22L)

for (yv in outcomes) {
  ylab <- gsub("[^A-Za-z0-9]", "_", outcome_labels[yv])

  p <- make_forest(
    results_list    = list(results_male, results_female),
    subgroup_levels = c("Male", "Female"),
    outcome         = yv,
    pal_col         = pal_mf,
    pal_shp         = shp_mf,
    title_str       = paste0("Linear regression (FC → ", outcome_labels[yv],
                             "): Male vs Female")
  )

  fname <- file.path(out_dir,
                     paste0("lm_forest_MaleVsFemale_", ylab, ".png"))
  ggsave(fname, plot = p, width = 14, height = 7, dpi = 300)
  cat("Saved:", basename(fname), "\n")
}

# ---- Section 15: Individual scatter plots for significant results (pFDR < .05) ----
# For every FC × outcome × subgroup combination that passes BH-FDR, saves a
# separate partial-regression scatter plot to sig_scatter/.

sig_out_dir <- file.path(out_dir, "sig_scatter")
dir.create(sig_out_dir, recursive = TRUE, showWarnings = FALSE)

pal_subgroup <- c(
  Full       = "#2CA02C",
  Male       = "#2CA02C",
  Female     = "#2CA02C",
  ASD        = "#D62728",
  NAC        = "#1F77B4",
  ASD_Male   = "#D62728",
  ASD_Female = "#D62728",
  NAC_Male   = "#1F77B4",
  NAC_Female = "#1F77B4"
)

subgroup_meta <- list(
  list(tbl = results_full,       data = df,            covs = c("age_c", "meanFD_c", "sex_F"), label = "Full"),
  list(tbl = results_male,       data = df_male,        covs = c("age_c", "meanFD_c"),          label = "Male"),
  list(tbl = results_female,     data = df_female,      covs = c("age_c", "meanFD_c"),          label = "Female"),
  list(tbl = results_asd,        data = df_asd,         covs = c("age_c", "meanFD_c", "sex_F"), label = "ASD"),
  list(tbl = results_nac,        data = df_nac,         covs = c("age_c", "meanFD_c", "sex_F"), label = "NAC"),
  list(tbl = results_asd_male,   data = df_asd_male,    covs = c("age_c", "meanFD_c"),          label = "ASD_Male"),
  list(tbl = results_asd_female, data = df_asd_female,  covs = c("age_c", "meanFD_c"),          label = "ASD_Female"),
  list(tbl = results_nac_male,   data = df_nac_male,    covs = c("age_c", "meanFD_c"),          label = "NAC_Male"),
  list(tbl = results_nac_female, data = df_nac_female,  covs = c("age_c", "meanFD_c"),          label = "NAC_Female")
)

for (sg in subgroup_meta) {
  sig_rows <- sg$tbl[!is.na(sg$tbl$pFDR) & sg$tbl$pFDR < 0.05, ]
  if (nrow(sig_rows) == 0) next

  for (i in seq_len(nrow(sig_rows))) {
    r       <- sig_rows[i, ]
    fc      <- r$x_var
    yv      <- r$y_var
    covs    <- sg$covs
    d_full  <- sg$data
    lbl     <- sg$label

    req <- c(fc, yv, covs)
    d   <- d_full[complete.cases(d_full[, req, drop = FALSE]), ]

    cov_rhs <- paste(covs, collapse = " + ")
    fc_res  <- residuals(lm(as.formula(paste(fc, "~", cov_rhs)), data = d))
    y_res   <- residuals(lm(as.formula(paste(yv, "~", cov_rhs)), data = d))
    plot_d  <- data.frame(fc_res    = as.numeric(fc_res),
                          y_res     = as.numeric(y_res),
                          diagnosis = factor(d$group, levels = c("ASD", "NAC")))

    col     <- pal_subgroup[lbl]
    fc_nice <- fc_label_map[fc]
    yv_nice <- outcome_labels[yv]
    stars   <- sig_stars(r$pFDR)

    p <- ggplot(plot_d, aes(x = fc_res, y = y_res)) +
      geom_point(aes(colour = diagnosis), alpha = 0.65, size = 2.2, stroke = 0.5) +
      geom_smooth(method = "lm", colour = col, fill = col,
                  alpha = 0.15, linewidth = 1.1, se = TRUE) +
      annotate("text", x = Inf, y = Inf, label = stars,
               hjust = 1.15, vjust = 1.5,
               size = 7, fontface = "bold", colour = col) +
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("β = %.3f  pFDR = %.4f", r$est, r$pFDR),
               hjust = -0.05, vjust = 1.7,
               size = 4.5, colour = col) +
      scale_colour_manual(name   = "Diagnosis",
                          values = c("ASD" = "#D62728", "NAC" = "#1F77B4")) +
      labs(
        title    = sprintf("%s  →  %s   [%s]", fc_nice, yv_nice, lbl),
        subtitle = sprintf("β = %.3f, 95%% CI [%.3f, %.3f], p = %.4f, pFDR = %.4f, n = %d",
                           r$est, r$ci_low, r$ci_high, r$pval, r$pFDR, r$n_used),
        x        = paste0(fc_nice, " (residualised)"),
        y        = paste0(scatter_y_labels[yv], " (residualised)")
      ) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

    fname <- file.path(
      sig_out_dir,
      sprintf("sig_%s_%s_%s.png",
              lbl, fc, gsub("[^A-Za-z0-9]", "_", yv))
    )
    ggsave(fname, plot = p, width = 6, height = 5, dpi = 300)
    cat("Saved:", basename(fname), "\n")
  }
}

cat("\n====== Done. All outputs written to:\n", out_dir, "\n")
cat("Sex-only forest plots written to:\n", out_dir, "\n")
sessionInfo()
