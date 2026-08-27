# =============================================================================
# Linear Regression: AASP Auditory score → SRS fso_com (direct path)
# Subgroups: full | males | females | ASD | NAC | ASD×sex | NAC×sex
# Covariates: age_c + meanFD_c + sex_F (mixed-sex groups)
#             age_c + meanFD_c only    (within-sex groups)
# Outputs: scatter plots with regression lines + forest plot of beta across subgroups
# =============================================================================

# ---- Section 1: Setup -------------------------------------------------------

pkgs <- c("readxl", "dplyr", "ggplot2", "writexl", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

set.seed(42)

data_file <- "/Users/james/Desktop/CYL/ASD_ThaCor/Score_sheet/drop_nosig_raw/Combine_All_cleaned_0.01-0.08_drop_nosig_raw_AuTha.xlsx"
out_dir   <- "/Users/james/Desktop/CYL/ASD_ThaCor/Score_sheet/drop_nosig_raw/ForGraph_lm_fsocom_auditory"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sig_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ""))))
}

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
    ),
    auditory_c = as.numeric(scale(auditory)),
    fso_com_z  = as.numeric(scale(fso_com)),
    age_c      = as.numeric(scale(age)),
    meanFD_c   = as.numeric(scale(meanFD))
  )

# ---- Section 3: Linear regression function ----------------------------------

run_lm_direct <- function(data, covariates = c("age_c", "meanFD_c", "sex_F")) {
  all_vars <- c("fso_com_z", "auditory_c", covariates)
  d <- data[complete.cases(data[, all_vars, drop = FALSE]), ]
  n_used <- nrow(d)

  if (n_used < length(all_vars) + 2) {
    message("  Skipping: n too small (", n_used, ")")
    return(NULL)
  }

  rhs <- paste(c("auditory_c", covariates), collapse = " + ")
  fit <- tryCatch(
    lm(as.formula(paste0("fso_com_z ~ ", rhs)), data = d),
    error = function(e) { message("  lm error: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)

  s  <- summary(fit)
  ci <- confint(fit, level = 0.95)

  if (!("auditory_c" %in% rownames(s$coefficients))) return(NULL)

  data.frame(
    est     = s$coefficients["auditory_c", "Estimate"],
    se      = s$coefficients["auditory_c", "Std. Error"],
    t_val   = s$coefficients["auditory_c", "t value"],
    pval    = s$coefficients["auditory_c", "Pr(>|t|)"],
    ci_low  = ci["auditory_c", 1],
    ci_high = ci["auditory_c", 2],
    r2      = s$r.squared,
    adj_r2  = s$adj.r.squared,
    n_used  = n_used,
    stringsAsFactors = FALSE
  )
}

# ---- Section 4: Run all subgroups -------------------------------------------

subgroup_defs <- list(
  Full       = list(data = df,                                                             covs = c("age_c", "meanFD_c", "sex_F")),
  Male       = list(data = df[df$sex_F == 0, ],                                            covs = c("age_c", "meanFD_c")),
  Female     = list(data = df[df$sex_F == 1, ],                                            covs = c("age_c", "meanFD_c")),
  ASD        = list(data = df[!is.na(df$group) & df$group == "ASD", ],                    covs = c("age_c", "meanFD_c", "sex_F")),
  NAC        = list(data = df[!is.na(df$group) & df$group == "NAC", ],                    covs = c("age_c", "meanFD_c", "sex_F")),
  ASD_Male   = list(data = df[!is.na(df$group) & df$group == "ASD" & df$sex_F == 0, ],   covs = c("age_c", "meanFD_c")),
  ASD_Female = list(data = df[!is.na(df$group) & df$group == "ASD" & df$sex_F == 1, ],   covs = c("age_c", "meanFD_c")),
  NAC_Male   = list(data = df[!is.na(df$group) & df$group == "NAC" & df$sex_F == 0, ],   covs = c("age_c", "meanFD_c")),
  NAC_Female = list(data = df[!is.na(df$group) & df$group == "NAC" & df$sex_F == 1, ],   covs = c("age_c", "meanFD_c"))
)

results_list <- lapply(names(subgroup_defs), function(nm) {
  cat(sprintf("\n===== %s (n = %d) =====\n", nm, nrow(subgroup_defs[[nm]]$data)))
  r <- run_lm_direct(subgroup_defs[[nm]]$data, subgroup_defs[[nm]]$covs)
  if (!is.null(r)) r$subgroup <- nm
  r
})
names(results_list) <- names(subgroup_defs)

all_results <- do.call(rbind, Filter(Negate(is.null), results_list))
all_results$subgroup <- factor(all_results$subgroup, levels = names(subgroup_defs))
all_results$pFDR    <- p.adjust(all_results$pval, method = "BH")

# ---- Section 5: Console summary ---------------------------------------------

cat("\n\n===== Summary: auditory_c → fso_com_z =====\n")
cat(sprintf("  %-12s  %7s  %7s  %9s  %9s  %7s  %7s  %6s\n",
            "Subgroup", "est", "SE", "CI_lo", "CI_hi", "p", "R2", "n"))
for (i in seq_len(nrow(all_results))) {
  r <- all_results[i, ]
  cat(sprintf("  %-12s  %7.4f  %7.4f  %9.4f  %9.4f  %7.4f  %6.4f  %6d\n",
              r$subgroup, r$est, r$se, r$ci_low, r$ci_high, r$pval, r$r2, r$n_used))
}

# ---- Section 6: Export results ----------------------------------------------

write_xlsx(
  list(all_subgroups = all_results),
  file.path(out_dir, "lm_fsocom_auditory_results.xlsx")
)
cat("\nResults saved: lm_fsocom_auditory_results.xlsx\n")

# ---- Section 7: Scatter plots per subgroup ----------------------------------
# Wong (2011) colorblind-safe palette
pal_all <- c(
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

make_scatter <- function(data, subgroup_name, colour, beta_row, covariates) {
  all_vars <- c("fso_com_z", "auditory_c", covariates)
  d <- data[complete.cases(data[, all_vars, drop = FALSE]), ]

  cov_rhs      <- paste(covariates, collapse = " + ")
  auditory_res <- as.numeric(residuals(lm(as.formula(paste("auditory_c ~", cov_rhs)), data = d)))
  fsocom_res   <- as.numeric(residuals(lm(as.formula(paste("fso_com_z ~", cov_rhs)),  data = d)))

  plot_d <- data.frame(
    auditory_res = auditory_res,
    fsocom_res   = fsocom_res,
    diagnosis    = factor(d$group, levels = c("ASD", "NAC"))
  )

  pval_str <- if (beta_row$pval < 0.001) "p < .001" else sprintf("p = %.3f", beta_row$pval)
  subtitle_str <- sprintf("β = %.3f, 95%% CI [%.3f, %.3f], %s, n = %d",
                          beta_row$est, beta_row$ci_low, beta_row$ci_high,
                          pval_str, beta_row$n_used)
  stars <- sig_stars(beta_row$pval)

  p <- ggplot(plot_d, aes(x = auditory_res, y = fsocom_res)) +
    geom_point(aes(colour = diagnosis), alpha = 0.65, size = 2) +
    geom_smooth(method = "lm", colour = colour, fill = colour, alpha = 0.15,
                linewidth = 1.1, se = TRUE) +
    scale_colour_manual(name   = "Diagnosis",
                        values = c("ASD" = "#D62728", "NAC" = "#1F77B4")) +
    labs(
      title    = sprintf("Auditory subscore → Social communication subscore: %s", subgroup_name),
      subtitle = subtitle_str,
      x        = "Auditory subscore (residualised)",
      y        = "Social communication subscore (residualised)"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  if (nchar(stars) > 0) {
    p <- p + annotate("text", x = Inf, y = Inf, label = stars,
                      hjust = 1.15, vjust = 1.5,
                      size = 7, fontface = "bold", colour = colour)
  }
  p <- p + annotate("text", x = -Inf, y = Inf,
                    label = sprintf("β = %.3f  pFDR = %.4f", beta_row$est, beta_row$pFDR),
                    hjust = -0.05, vjust = 1.7,
                    size = 4.5, colour = colour)
  p
}

scatter_plots <- lapply(names(subgroup_defs), function(nm) {
  row <- all_results[all_results$subgroup == nm, ]
  if (nrow(row) == 0) return(NULL)
  make_scatter(subgroup_defs[[nm]]$data, nm, pal_all[nm], row, subgroup_defs[[nm]]$covs)
})
names(scatter_plots) <- names(subgroup_defs)

# Save individual scatter plots
for (nm in names(scatter_plots)) {
  if (is.null(scatter_plots[[nm]])) next
  ggsave(
    file.path(out_dir, sprintf("scatter_%s.png", nm)),
    plot = scatter_plots[[nm]],
    width = 7, height = 5.5, dpi = 300
  )
  cat(sprintf("Saved: scatter_%s.png\n", nm))
}

# Combined panel: 3×3 grid (Full + sex + diag + diag×sex)
valid_plots <- Filter(Negate(is.null), scatter_plots)
combined <- wrap_plots(valid_plots, ncol = 3) +
  plot_annotation(
    title   = "Auditory subscore → Social communication subscore: all subgroups",
    caption = "Regression line with 95% CI band; β from lm(fso_com_z ~ auditory_c + covariates)"
  )

ggsave(
  file.path(out_dir, "scatter_all_subgroups.png"),
  plot = combined,
  width = 21, height = 15, dpi = 300
)
cat("Saved: scatter_all_subgroups.png\n")

# ---- Section 8: Forest plot of beta across subgroups -----------------------

forest_df <- all_results |>
  mutate(
    subgroup  = factor(subgroup, levels = rev(names(subgroup_defs))),
    sig_point = pval < 0.05,
    pt_fill   = ifelse(sig_point, pal_all[as.character(subgroup)], "white"),
    stars     = sig_stars(pval)
  )

forest_plot <- ggplot(forest_df, aes(x = est, y = subgroup)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high, colour = subgroup),
    height = 0.35, linewidth = 0.7
  ) +
  geom_point(
    aes(colour = subgroup, fill = pt_fill),
    shape = 21, size = 4, stroke = 1
  ) +
  geom_text(
    data = ~ subset(.x, stars != ""),
    aes(label = stars, colour = subgroup),
    size = 5, hjust = 0.5, vjust = -0.7,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = sprintf("β=%.3f, p=%.3f, n=%d", est, pval, n_used)),
    hjust = -0.08, size = 3, colour = "grey30"
  ) +
  scale_colour_manual(values = pal_all, guide = "none") +
  scale_fill_identity() +
  labs(
    title    = "Auditory subscore → Social communication subscore: β across subgroups",
    subtitle = "95% CI (t-distribution); filled point = p < .05 (uncorrected);  * p<.05   ** p<.01   *** p<.001",
    x        = "β (standardised), 95% CI",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(out_dir, "forest_beta_subgroups.png"),
  plot  = forest_plot,
  width = 14, height = 7, dpi = 300
)
cat("Saved: forest_beta_subgroups.png\n")

cat("\n====== Done. All outputs written to:\n", out_dir, "\n")
sessionInfo()
