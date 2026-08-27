install.packages(c("MatchIt", "readxl", "openxlsx"))
library(MatchIt)
library(readxl)
library(openxlsx)

# ── Load & clean ─────────────────────────────────────────────────────────────
abide_raw   <- read_excel("../profile_ABIDE+NTUH.xlsx")
abide_clean <- subset(abide_raw, Sex %in% c(1, 2) & (!is.na(FSIQ) | Site == "NTUH") &
                        !(HighMotion %in% "x") & !(HighMeanFD %in% "x"))
# GROUP: 1 = ASD, 2 = Control  |  Sex: 1 = Male, 2 = Female
abide_clean$is_ASD    <- as.integer(abide_clean$GROUP == 1)
abide_clean$is_Female <- as.integer(abide_clean$Sex   == 2)

set.seed(123)

# ── Shared MatchIt settings ───────────────────────────────────────────────────
match_covs    <- ~ Age + MeanFD   # RHS covariates used in all three matchings
match_caliper <- c(Age = 3, MeanFD = 0.1)
match_exact   <- ~Site                 # set to NULL or ~Site to disable or use exact site matching
match_std_cal <- FALSE

# ── Anchor: ALL female ASD — never filtered ───────────────────────────────────
all_female_asd <- subset(abide_clean, Sex == 2 & GROUP == 1)
female_asd_ids <- all_female_asd$ID
cat(sprintf("\nAnchor: %d female ASD subjects (all retained)\n", nrow(all_female_asd)))

# ── Matching 1: Trim female TD pool to match female ASD ───────────────────────
# Flip treatment: female TD = treated (gets trimmed), female ASD = control pool (kept)
females       <- subset(abide_clean, Sex == 2)
females$is_TD <- as.integer(females$GROUP == 2)

m1 <- matchit(
  formula     = update(match_covs, is_TD ~ .),
  data        = females,
  method      = "nearest",
  distance    = "glm",
  ratio       = 1,
  exact       = match_exact,
  caliper     = match_caliper,
  std.caliper = match_std_cal
)
cat("\n=== Matching 1: Female TD trimmed to match Female ASD (1:1) ===\n")
print(summary(m1, un = FALSE))

# Only extract the matched female TD; female ASD anchor is taken in full separately
female_ctrl_ids <- match.data(m1, group = "treated")$ID

# ── Matching 2: Select male ASD matched to female ASD (3:1) ──────────────────
# Female ASD = treated (anchor, ratio drives 3 males per female)
# Male ASD   = control pool (gets trimmed)
# No exact = ~ Site: sex-balancing within diagnosis across sites is acceptable,
# and avoids dropping female ASD anchors at sites with too few males.
asd_data            <- subset(abide_clean, GROUP == 1)
asd_data$is_Female  <- as.integer(asd_data$Sex == 2)

m2 <- matchit(
  formula     = update(match_covs, is_Female ~ .),
  data        = asd_data,
  method      = "nearest",
  distance    = "glm",
  ratio       = 3,
  exact       = match_exact,
  caliper     = match_caliper,
  std.caliper = match_std_cal
)
cat("\n=== Matching 2: Male ASD selected to match Female ASD (3:1) ===\n")
print(summary(m2, un = FALSE))

male_asd_ids <- match.data(m2, group = "control")$ID

# ── Matching 3: Select male TD matched to female TD (3:1) ────────────────────
# Female TD  = treated (anchor); male TD = control pool (gets trimmed)
# No exact = ~ Site for the same reason as Matching 2.
ctrl_data           <- subset(abide_clean, GROUP == 2 &
                                (Sex == 1 | ID %in% female_ctrl_ids))
ctrl_data$is_Female <- as.integer(ctrl_data$Sex == 2)

m3 <- matchit(
  formula     = update(match_covs, is_Female ~ .),
  data        = ctrl_data,
  method      = "nearest",
  distance    = "glm",
  ratio       = 3,
  exact       = match_exact,
  caliper     = match_caliper,
  std.caliper = match_std_cal
)
cat("\n=== Matching 3: Male TD selected to match Female TD (3:1) ===\n")
print(summary(m3, un = FALSE))

male_ctrl_ids <- match.data(m3, group = "control")$ID

# ── Combine into final dataset ────────────────────────────────────────────────
all_ids      <- unique(c(female_asd_ids, female_ctrl_ids, male_asd_ids, male_ctrl_ids))
matched_data <- subset(abide_clean, ID %in% all_ids)

cat("\n=== Final matched sample ===\n")
print(table(Sex = matched_data$Sex, GROUP = matched_data$GROUP))

# ── Sorted subject list: GROUP → Sex → Site ──────────────────────────────────
matched_sorted <- matched_data[order(matched_data$GROUP, matched_data$Sex, matched_data$Site), ]

# ── Summary + pairwise t-tests across sex × group subgroups ──────────────────
sg_def <- list(
  list(label = "Female ASD",     sex = 2, group = 1),
  list(label = "Female Control", sex = 2, group = 2),
  list(label = "Male ASD",       sex = 1, group = 1),
  list(label = "Male Control",   sex = 1, group = 2)
)
vars <- intersect(c("MeanFD", "Age", "FSIQ"), names(matched_data))

# Only same-sex or same-diagnosis comparisons (exclude cross-sex cross-diagnosis)
# Kept: Female ASD vs Female Control, Male ASD vs Male Control,
#        Female ASD vs Male ASD, Female Control vs Male Control
pairs <- list(c(1,2), c(3,4), c(1,3), c(2,4))

summary_rows <- do.call(rbind, lapply(pairs, function(idx) {
  sg_a <- sg_def[[idx[1]]]; sg_b <- sg_def[[idx[2]]]
  da   <- subset(matched_data, Sex == sg_a$sex & GROUP == sg_a$group)
  db   <- subset(matched_data, Sex == sg_b$sex & GROUP == sg_b$group)
  do.call(rbind, lapply(vars, function(v) {
    tt <- t.test(da[[v]], db[[v]])
    data.frame(
      Group_A     = sg_a$label,
      N_A         = nrow(da),
      Mean_A      = round(mean(da[[v]], na.rm = TRUE), 3),
      SD_A        = round(sd(da[[v]],   na.rm = TRUE), 3),
      Group_B     = sg_b$label,
      N_B         = nrow(db),
      Mean_B      = round(mean(db[[v]], na.rm = TRUE), 3),
      SD_B        = round(sd(db[[v]],   na.rm = TRUE), 3),
      Variable    = v,
      t_statistic = round(tt$statistic, 3),
      df          = round(tt$parameter, 1),
      p_value     = round(tt$p.value,   4),
      stringsAsFactors = FALSE
    )
  }))
}))

cat("\n=== Pairwise t-tests ===\n")
print(summary_rows)

# ── NTUH case counts in matched sample ───────────────────────────────────────
ntuh_matched <- subset(matched_data, Site == "NTUH")
ntuh_summary <- data.frame(
  NTUH  = c("fa", "fc", "ma", "mc"),
  n     = c(
    nrow(subset(ntuh_matched, Sex == 2 & GROUP == 1)),
    nrow(subset(ntuh_matched, Sex == 2 & GROUP == 2)),
    nrow(subset(ntuh_matched, Sex == 1 & GROUP == 1)),
    nrow(subset(ntuh_matched, Sex == 1 & GROUP == 2))
  ),
  stringsAsFactors = FALSE
)
cat("\n=== NTUH subjects in matched sample ===\n")
print(ntuh_summary)

# ── Recode Sex/GROUP for export: Female=0, Male=1; ASD=1, Control=0 ──────────
recode_for_export <- function(df) {
  df$Sex   <- as.integer(df$Sex   == 1)   # Male=1, Female=0
  df$GROUP <- as.integer(df$GROUP == 1)   # ASD=1,  Control=0
  df
}
matched_data_out   <- recode_for_export(matched_data)
matched_sorted_out <- recode_for_export(matched_sorted)

# ── Export to Excel ───────────────────────────────────────────────────────────
wb <- createWorkbook()

addWorksheet(wb, "matched_all")
writeData(wb, "matched_all", matched_data_out)

addWorksheet(wb, "by_group")
writeData(wb, "by_group", matched_sorted_out)

addWorksheet(wb, "summary")

# Row 1: MatchIt settings above the column-header row (derived from shared settings)
# Placed at Group_A (col 1), Group_B (col 5), Variable (col 9) to align with data below
match_vars_str <- sprintf("[%s]", deparse(match_covs[[2]]))
caliper_str    <- sprintf("[%s]", paste(names(match_caliper), match_caliper, sep = " = ", collapse = ", "))
exact_str      <- if (is.null(match_exact)) "[no exact site]" else
                    sprintf("[exact: %s]", paste(all.vars(match_exact), collapse = " + "))

writeData(wb, "summary", match_vars_str, startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, "summary", caliper_str,    startRow = 1, startCol = 5, colNames = FALSE)
writeData(wb, "summary", exact_str,      startRow = 1, startCol = 9, colNames = FALSE)

# Row 2+: summary t-test table (headers at row 2, data from row 3)
writeData(wb, "summary", summary_rows, startRow = 2, startCol = 1)

# NTUH summary: 2-column table to the right (col 14-15, col 13 left as gap)
writeData(wb, "summary", ntuh_summary, startRow = 2, startCol = 14)

# ── Build filename suffix from MatchIt settings ───────────────────────────────
covariate_abbrev <- c(Age = "A", FSIQ = "F", MeanFD = "M")
cov_vars   <- all.vars(match_covs)
cov_str    <- paste(covariate_abbrev[cov_vars[cov_vars %in% names(covariate_abbrev)]], collapse = "")
age_cal    <- if ("Age"    %in% names(match_caliper)) match_caliper["Age"]    else NA
meanfd_cal <- if ("MeanFD" %in% names(match_caliper)) match_caliper["MeanFD"] else NA
exact_flag <- if (is.null(match_exact)) "n" else "y"
file_suffix <- sprintf("_%s_%s_%s_%s", cov_str,
                       ifelse(is.na(age_cal),    "NA", format(age_cal,    scientific = FALSE)),
                       ifelse(is.na(meanfd_cal), "NA", format(meanfd_cal, scientific = FALSE)),
                       exact_flag)
out_file <- sprintf("../profile_ABIDE+fullNTUH_matched%s.xlsx", file_suffix)

saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("\nSaved: %s\n", out_file))
