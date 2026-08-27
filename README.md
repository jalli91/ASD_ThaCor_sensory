# ASD_ThaCor_sensory
neuroimaging and behavioral analysis scripts for manuscript: Sexual Dimorphism in Thalamocortical Connectivity Underlies Auditory-Mediated Social Difficulties in Autistic Adults

> Ta-Cheng Lin, Tsung-Ren Huang, Min Liu, Li-Ying Fan, Susan Shur-Fen Gau, Yi-Ling Chien.
> *Sexual Dimorphism in Thalamocortical Connectivity Underlies Auditory-Mediated Social
> Difficulties in Autistic Adults.*
> Correspondence: Yi-Ling Chien, MD, PhD (ylchien71@ntu.edu.tw), Department of Psychiatry,
> National Taiwan University Hospital and College of Medicine, Taipei, Taiwan.

This repository is provided for **reviewer reference**: it documents exactly how the reported
results were produced. It is a read-only mirror of the analysis scripts — see
[Reproducibility notes](#reproducibility-notes) before attempting to run anything.

---

## Study in brief

Resting-state fMRI seed-based functional connectivity (FC) between the bilateral thalamus and
sensory-regulating cortical ROIs (auditory, somatosensory, visual, gustatory), tested for
diagnosis, sex, and diagnosis × sex effects in autistic adults.

Two samples:

| Sample | N | Description |
|---|---|---|
| **NTUH cohort** (primary) | 96 (47 ASD / 49 non-autistic controls; all Asian) | Recruited at the Adult Autism Clinic, National Taiwan University Hospital. Single Siemens MAGNETOM Prisma 3T scanner. |
| **ABIDE combined sample** (corroboration) | 189 (103 ASD / 86 NAC; 65 female) | NTUH cohort combined with adults from ABIDE I & II, assembled by propensity-score matching. |

Headline findings reproduced by this code: thalamus–superior temporal gyrus (STG)
hyperconnectivity in ASD, a group × sex interaction localising that effect to females,
female-specific brain–behaviour correlations with AASP auditory and SRS social communication
scores, and an exploratory mediation of FC → auditory sensory features → social communication
difficulties.

## Repository contents

```
.
└── code/
    ├── ABIDE/              ABIDE arm + the combined ABIDE+NTUH models
    │   ├── S1_ABIDE_ASD_preprocess_multisite.m
    │   ├── S2_ABIDE_motion_profiles_report.m
    │   ├── S3_ABIDE_MatchIt_NTUH_skipIQ.R
    │   ├── S4_ABIDE_ThaCor_GLM_NTUH_loop_no_site.m
    │   └── rescue/         Re-run scripts for subjects that failed the main pass
    └── NTUH/               NTUH-only arm
        ├── S1_CYL_batch_denoising_analysis_loop_for_ASD.m
        ├── S2_CYL_FC_difference_GLM_fitlm.m
        ├── S3_CYL_FC_difference_GLM_fitlm_sex.m
        ├── S4_CYL_FC_difference_raw.m
        ├── S5_ForGraph_lm_FC.R
        ├── S6_ForGraph_lm_fsocom_auditory.R
        ├── S7_ForGraph_sim_AuM_lm_lavaan.R
        ├── sub_extract.m       FDR correction + heat-map helper
        ├── motion_analx.m      Motion QC
        └── outlier_summary.m   ART scrubbing summary
```

Scripts are prefixed `S1…S7` in intended execution order within each arm. The `S`-numbering is
a reading order for reviewers, not a build system — each script is run on its own.

---

## Pipeline and mapping to the manuscript

### NTUH arm (primary analyses)

| Script | What it does | Manuscript element |
|---|---|---|
| `NTUH/S1_CYL_batch_denoising_analysis_loop_for_ASD.m` | CONN denoising + first-level ROI-to-ROI analysis per subject (band-pass 0.01–0.08 Hz; WM/CSF principal components, six realignment parameters, ART scrubbing regressors, rest-condition effects). | Methods → *Image acquisition and preprocessing* |
| `NTUH/motion_analx.m`, `NTUH/outlier_summary.m` | Mean framewise displacement (rotations converted to mm on a 50-mm sphere) and scrubbed-frame proportions; used for the mean-FD exclusion (> 2 SD above sample mean) and the < 10 % scrubbing check. | Methods → *Participants*; Table 1, Table S1 |
| `NTUH/S2_CYL_FC_difference_GLM_fitlm.m` | Whole-sample GLM per ROI pair: `FC ~ age + meanFD + sex + group + sex×group`. Extracts β, 95 % CI and p for the group, sex, and interaction terms across thalamus × {auditory, somatosensory, visual, insular} ROI sets. | Figure 2; Tables S3–S5 |
| `NTUH/S3_CYL_FC_difference_GLM_fitlm_sex.m` | The same GLM re-fitted within each sex subgroup: `FC ~ age + meanFD + group`. | Figure 3; Tables S6–S7 |
| `NTUH/S4_CYL_FC_difference_raw.m` | Exports raw (Fisher-z) FC values per subject per ROI pair — the input to all brain–behaviour and mediation models. | Input to S5–S7 |
| `NTUH/S5_ForGraph_lm_FC.R` | Standardised linear regressions of FC on AASP auditory subscore and SRS Social Communication subscore, across subgroups (full / male / female / ASD / NAC / ASD×sex / NAC×sex). Covariates: age, mean FD, and sex in mixed-sex groups; age and mean FD only within sex. BH-FDR across the eight thalamus–STG pairs. Emits forest and scatter plots. | Figures 4 and 5 |
| `NTUH/S6_ForGraph_lm_fsocom_auditory.R` | The direct AASP auditory → SRS Social Communication path, same subgroup and covariate scheme. | Figure 6 |
| `NTUH/S7_ForGraph_sim_AuM_lm_lavaan.R` | Simple mediation (Hayes Model 4) in the female subgroup: X = FC, M = AASP auditory, Y = SRS Social Communication. Paths `a`, `b`, `c'`, `c` by OLS with t-distribution CIs; indirect effect by `lavaan::sem(se = "bootstrap", bootstrap = 5000)` with percentile CIs. BH-FDR applied to all four paths across the eight FC pairs. | Table 2 |

### ABIDE arm (corroboration in the combined sample)

| Script | What it does | Manuscript element |
|---|---|---|
| `ABIDE/S1_ABIDE_ASD_preprocess_multisite.m` | Full CONN pipeline for ABIDE I/II subjects: setup + preprocessing (`default_mni`: realign → coregister → segment → normalise to 2 mm MNI → 6 mm FWHM smoothing), denoising (0.01–0.08 Hz), and first-level bivariate-correlation analysis. Contains a per-site parameter library (TR, TE, slice count, slice order, multiband factor) transcribed from each site's published protocol, applied individually during slice-timing correction and CONN project setup. | Methods → *Image acquisition and preprocessing* (ABIDE paragraph) |
| `ABIDE/rescue/*.m` | Targeted re-runs for subjects that were missing or produced empty first-level output in the main pass (`ABIDE_OILH_preprocess_rescue.m`, `ABIDE_empty_analysis_rescue.m`, `ABIDE_IP_rerun.m`). Same pipeline, restricted subject lists. | — |
| `ABIDE/S2_ABIDE_motion_profiles_report.m` | Builds the motion/phenotype profile table (mean FD, scrubbed frames, scrub proportion, age, FIQ, site one-hot columns) that feeds the matching step. | Input to S3 |
| `ABIDE/S3_ABIDE_MatchIt_NTUH_skipIQ.R` | Propensity-score matching that constructs the ABIDE combined sample (`MatchIt`, `method = "nearest"`, `distance = "glm"`). See [Matching design](#matching-design) below. | Supplementary → *Selection of ABIDE combined sample*; Table S2 |
| `ABIDE/S4_ABIDE_ThaCor_GLM_NTUH_loop_no_site.m` | Re-runs the full GLM pipeline on the matched combined sample: whole-sample group / sex / group×sex effects plus within-sex subgroup group effects, for every matched profile file found. Exports one Excel workbook per profile with `All` / `Female` / `Male` sheets. | Figure S1; Tables S8–S9 |

### Matching design

Because autistic females are scarce across ABIDE sites, **all 40 female ASD participants are
retained as anchors** and every other cell is matched to them (`ABIDE/S3_…skipIQ.R`):

1. Female NAC matched 1:1 to female ASD.
2. Male ASD selected 3:1 against female ASD.
3. Male NAC selected 3:1 against female NAC.

Settings: covariates `~ Age + MeanFD`; calipers 3 years and 0.1 mm (unstandardised);
`exact = ~ Site`; `set.seed(123)`. The 3:1 ratio reflects the male:female proportion in the
NTUH cohort. Full-scale IQ is deliberately **not** a matching covariate — hence the
`skipIQ` suffix — because it is missing for part of the ABIDE sample. Participants with
FSIQ < 70, age < 18, or high-motion flags are removed beforehand.

Resulting sample: 40 female ASD, 25 female NAC, 63 male ASD, 61 male NAC (N = 189; 82 from
ABIDE I, 40 from ABIDE II, 67 from NTUH).

The script writes `profile_ABIDE+fullNTUH_matched<suffix>.xlsx`, where the suffix encodes the
matching parameters as `_<covariate initials>_<age caliper>_<meanFD caliper>_<y|n exact site>`
(A = Age, F = FSIQ, M = MeanFD). The configuration above yields **`_AM_3_0.1_y`**. `S4` globs
these files and fits every profile it finds, which is how matching-parameter sensitivity checks
were run.

---

## Regions of interest

Thalamic seeds are the bilateral thalamus from the FSL Harvard–Oxford subcortical atlas
(CONN's built-in `Default_atlas`). Cortical ROIs are drawn from the Harvard–Oxford cortical
atlas and CONN's HCP ICA-derived networks, and grouped by sensory modality in the analysis
scripts:

| Modality | ROIs |
|---|---|
| Auditory | aSTG L/R, pSTG L/R, Heschl's gyrus L/R |
| Somatosensory | Postcentral gyrus, superior parietal lobule, anterior/posterior supramarginal gyrus, angular gyrus, parietal operculum (L/R) |
| Visual | Occipital, medial and lateral visual networks |
| Gustatory | Insular cortex L/R |

ROI indices are recovered by name-matching against the `names` field of each subject's CONN
`resultsROI_Subject001_Condition001.mat`, then relabelled to the short names used in the
figures and tables.

## Statistical settings

- Bivariate correlations, Fisher z-transformed, taken from CONN first-level output.
- GLMs fitted with MATLAB `fitlm`; 95 % confidence intervals from `coefCI(model, 0.05)`.
- Multiple comparisons: Benjamini–Hochberg FDR at **q = 0.05**, applied **within each sensory
  modality** for the FC GLMs (`sub_extract.m`), and across the eight thalamus–STG pairs for the
  brain–behaviour and mediation analyses.
- Brain–behaviour models use z-standardised FC, age, mean FD, and clinical scores, so reported
  β are standardised coefficients.
- Mediation: 5,000 bootstrap resamples, percentile CIs; an indirect-effect CI excluding zero is
  the criterion for significance.
- Seeds: `set.seed(123)` in the matching script, `set.seed(42)` in the R analysis/plotting
  scripts. The MATLAB scripts are deterministic.

## Variable coding

Group and sex coding **changes at the R → MATLAB boundary** (function `recode_for_export()` in
`ABIDE/S3_ABIDE_MatchIt_NTUH_skipIQ.R`). Reading a β with the wrong convention inverts its sign:

| Stage | GROUP | Sex |
|---|---|---|
| Raw profile sheets and R scripts | 1 = ASD, 2 = Control | 1 = Male, 2 = Female |
| Exported `matched_all` sheet and all MATLAB scripts | **1 = ASD, 0 = Control** | **1 = Male, 0 = Female** |

The MATLAB convention (ASD = 1, NAC = 0; Male = 1, Female = 0) is the one stated in the figure
legends. Under it, β<sub>group</sub> > 0 means hyperconnectivity in ASD, and
β<sub>group×sex</sub> < 0 means that hyperconnectivity is attenuated in males.

Subject folder names carry group and sex:

- **ABIDE arm** — `[m|f][a|c]<numeric ID>_session_<n>`; `m`/`f` = male/female, `a`/`c` =
  ASD/control (e.g. `ma0050952`). Where a subject has more than one session, only session 1 is
  used by `S4`.
- **NTUH arm** — CONN output folders `conn_1*` = control, `conn_2*` = ASD, with sex read from
  a separate `gender_index.xlsx` (M = 1, F = 0).

---

## Software requirements

**MATLAB** (Statistics and Machine Learning Toolbox required for `fitlm` / `coefCI`), with
**SPM12** ([RRID:SCR_002592](https://www.nitrc.org/projects/spm)) and the
**CONN toolbox** v21+ ([RRID:SCR_009550](https://www.conn-toolbox.org/)) on the path. CONN v21+
is recommended for native `.nii.gz` support; earlier versions need the ABIDE images
decompressed first.

**R 4.5.3**, with `MatchIt` 4.7.2 ([RRID:SCR_025618](https://CRAN.R-project.org/package=MatchIt)),
`lavaan` ([RRID:SCR_027663](https://cran.r-project.org/package=lavaan)), `readxl`, `openxlsx`,
`writexl`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`.

`sub_extract.m` calls `fdr_BH`, a Benjamini–Hochberg implementation that must also be on the
MATLAB path.

## Data availability

No participant data are included in this repository.

- **NTUH cohort** — data supporting the findings are available on reasonable request from the
  corresponding author. Some data cannot be made public because of privacy and ethical
  restrictions.
- **ABIDE I / II** — publicly available from
  <https://fcon_1000.projects.nitrc.org/indi/abide/>, under those repositories' own terms.

Ethics: approved by the Research Ethics Committee of National Taiwan University Hospital
(REC no. 201903126RINA); conducted in accordance with the Declaration of Helsinki; written
informed consent obtained from all participants.

---

## Reproducibility notes

Please read these before judging the code as runnable — several are deliberate properties of
the working environment rather than oversights.

1. **Paths are absolute and hardcoded** to the external volumes the analyses were run from
   (`/Volumes/Milk/CYL/ASD_ThaCor/…` for the ABIDE arm, `/Volumes/Sea/CYL/ASD_ThaCor/NTUH/…`
   for the NTUH arm; some earlier NTUH scripts still point at `Milk`, and the R plotting
   scripts read from a local `~/Desktop/CYL/…` copy). They are preserved verbatim so that the
   provenance of every input and output file is visible. Running the pipeline elsewhere
   requires editing the `dataDir` / `data_file` / `out_dir` lines at the top of each script.

2. **Scripts are run interactively** — `.m` files from the MATLAB desktop, `.R` files from
   RStudio. This is why `S3_ABIDE_MatchIt_NTUH_skipIQ.R` opens with an unconditional
   `install.packages()`. There is no build, test, or driver script, and none is implied.

3. **Preprocessing is long-running and not resumable.** `S1_ABIDE_ASD_preprocess_multisite.m`
   runs hours to days per site, writes gigabytes of CONN output, and begins with `clear all`
   without resume logic. The `rescue/` scripts and the hand-edited restart index in
   `NTUH/S1_…_loop_for_ASD.m` (`for subNum = 71:length(useable_data)`) exist for exactly that
   reason and record where the original runs were resumed.

4. **`sub_extract.m` exists in two incompatible signatures.** The copy in this repository
   (`NTUH/sub_extract.m`) takes 5 arguments and is the one called by `NTUH/S2` and `NTUH/S3`.
   `ABIDE/S4` calls a 6-argument variant with a leading table-name argument, which lived
   alongside the ABIDE scripts on the analysis volume and is not mirrored here. The two differ
   only in that added label argument; the FDR logic is identical.

5. **`S4` is the no-site variant.** Its covariate matrix is allocated with 19 columns, but only
   columns 1–5 (`age`, `meanFD`, `sex`, `group`, `sex×group`) are populated; the site-dummy
   columns are left as zeros. This keeps `fitlm` coefficient indices aligned with the
   site-included variants of the model while excluding site from the fit. Coefficient 4 is sex,
   5 is group, 6 is the interaction.

6. **A stale comment in the matching script.** The inline comments above matchings 2 and 3 in
   `S3_ABIDE_MatchIt_NTUH_skipIQ.R` read "No exact = ~ Site", but both calls pass
   `exact = match_exact`, i.e. `~ Site`. The executed behaviour — exact matching within
   acquisition site throughout — is what the Supplementary Materials describe; the comment is
   left over from an earlier configuration and does not reflect the code.

7. **Traditional-Chinese identifiers.** The R analysis scripts read the Excel sheet named
   `工作表1` (the default Chinese sheet name), and some source filenames on the analysis volume
   contain Chinese characters.

## Funding

Supported by grants from the Ministry of Science and Technology
(113-2314-B-002-186-MY2, 114-2410-H-002-117-MY3), National Taiwan University Hospital
(VN113-01, 113-S0148, 114-CTC0011), and the National Health Research Institute
(NHRI-EX110-110008PC, NHRI-EX111-110008PC, NHRI-EX112-110008PC, NHRI-EX113-110008PC), Taiwan.

We thank the participants and their families for their contribution to this study, and the
ABIDE consortium for making their data publicly available.
