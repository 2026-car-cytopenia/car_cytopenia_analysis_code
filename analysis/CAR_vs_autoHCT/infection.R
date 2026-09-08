#!/usr/bin/env Rscript
# Run: Rscript this_script.R /path/to/private/settings.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 "..", "_shared", "settings.R"))
configure_study("infection")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(MatchIt)
  library(cobalt)
  library(survival)
  library(survminer)
  library(cmprsk)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(janitor)
  library(broom)
})

`%notin%` <- Negate(`%in%`)

# =========================================================
# CONFIG
# =========================================================
PATH_INF    <- study_path("PATH_INF")
PATH_MATCH  <- study_path("PATH_MATCH")
PATH_CYTO   <- study_path("PATH_CYTO")

out_dir <- study_path("out_dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(...) file.path(out_dir, ...)

figshare_dir <- file.path(out_dir, "figshare_datasets")
dir.create(figshare_dir, showWarnings = FALSE, recursive = TRUE)
figshare_file <- function(...) file.path(figshare_dir, ...)

theme_set(theme_classic(base_family = "Arial", base_size = 16))

pal_treatment <- c("HCT" = "#377EB8", "CAR" = "#E41A1C")

manifest <- tibble::tibble(
  file = character(),
  rows = integer(),
  columns = integer(),
  description = character()
)

# =========================================================
# HELPERS
# =========================================================
safe_read_delim <- function(path, delim = NULL, ...) {
  if (!file.exists(path)) stop("Missing required input file: ", path)
  if (is.null(delim)) {
    readr::read_delim(path, show_col_types = FALSE, guess_max = 100000, progress = FALSE, ...)
  } else {
    readr::read_delim(path, delim = delim, show_col_types = FALSE, guess_max = 100000, progress = FALSE, ...)
  }
}

as_num_days <- function(x) {
  if (inherits(x, "difftime")) return(as.numeric(x, units = "days"))
  if (inherits(x, "Date")) stop("Date provided; subtract index date first.")
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) return(readr::parse_number(x))
  as.numeric(x)
}

fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

first_existing <- function(candidates, nm) {
  hit <- intersect(candidates, nm)
  if (length(hit) > 0) return(hit[[1]])
  hit_clean <- intersect(janitor::make_clean_names(candidates), janitor::make_clean_names(nm))
  if (length(hit_clean) > 0) return(nm[match(hit_clean[[1]], janitor::make_clean_names(nm))])
  NA_character_
}

extract_lymid <- function(x) {
  stringr::str_extract(as.character(x), "LYM[0-9]+")
}

make_id_map <- function(raw_ids, prefix = "Subject") {
  raw_chr <- as.character(raw_ids)
  raw_chr <- raw_chr[!is.na(raw_chr) & raw_chr != ""]
  raw_chr <- sort(unique(raw_chr))
  tibble::tibble(raw_id = raw_chr) %>%
    dplyr::mutate(
      lymid = extract_lymid(raw_id),
      public_id = ifelse(
        !is.na(lymid) & lymid != "",
        lymid,
        paste0(prefix, "_", sprintf("%03d", dplyr::row_number()))
      )
    ) %>%
    dplyr::select(raw_id, public_id)
}

lookup_public_id <- function(raw_ids, id_map, prefix = "Subject") {
  raw_chr <- as.character(raw_ids)
  out <- id_map$public_id[match(raw_chr, id_map$raw_id)]
  missing <- is.na(out) & !is.na(raw_chr) & raw_chr != ""
  if (any(missing)) {
    # Deterministic fallback; should rarely be used if maps are built correctly.
    tmp <- make_id_map(raw_chr[missing], prefix = prefix)
    out[missing] <- tmp$public_id[match(raw_chr[missing], tmp$raw_id)]
  }
  out
}

remove_internal_id_columns <- function(df) {
  bad_exact <- c(
    "MRN", "mrn", "SID", "sid",
    "Sample_ID", "sample_id", "SAMPLE_ID", "SampleID", "sampleID",
    "CCT_ID", "cct_id", "CCTN", "cctn",
    "Trial_ID", "trial_id", "UNIQUE.ID",
    "Specimen_ID", "Specimen", "raw_sample", "raw_subject_id", "raw_pair_id",
    "raw_cct_id", "raw_id", "PID", "pid"
  )
  df %>% dplyr::select(-dplyr::any_of(bad_exact))
}

audit_no_disallowed_ids <- function(df, label) {
  if (nrow(df) == 0 || ncol(df) == 0) return(invisible(TRUE))

  chr_df <- df %>% dplyr::select(where(~is.character(.x) || is.factor(.x)))
  if (ncol(chr_df) == 0) return(invisible(TRUE))

  vals <- unlist(lapply(chr_df, as.character), use.names = FALSE)
  vals <- vals[!is.na(vals)]

  # LYM/AAID identifiers are intentionally allowed.
  bad_pattern <- study_setting("study_value_001")
  bad_vals <- unique(vals[stringr::str_detect(vals, bad_pattern)])

  if (length(bad_vals) > 0) {
    stop(
      "Disallowed raw identifier detected in ", label, ": ",
      paste(utils::head(bad_vals, 20), collapse = ", "),
      ifelse(length(bad_vals) > 20, " ...", "")
    )
  }

  invisible(TRUE)
}

write_figshare_csv <- function(df, filename, description) {
  df_out <- df %>% remove_internal_id_columns()
  audit_no_disallowed_ids(df_out, filename)
  readr::write_csv(df_out, figshare_file(filename), na = "")

  manifest <<- dplyr::bind_rows(
    manifest,
    tibble::tibble(
      file = filename,
      rows = nrow(df_out),
      columns = ncol(df_out),
      description = description
    )
  )

  invisible(df_out)
}

save_plot_all <- function(p, basename, width = 7, height = 5, dpi = 400) {
  ggsave(out_file(paste0(basename, ".pdf")), plot = p, width = width, height = height, device = cairo_pdf)
  ggsave(out_file(paste0(basename, ".png")), plot = p, width = width, height = height, dpi = dpi)
  try(ggsave(out_file(paste0(basename, ".eps")), plot = p, width = width, height = height, device = cairo_ps), silent = TRUE)
  invisible(p)
}

coerce_treatment <- function(x) {
  z <- stringr::str_to_upper(stringr::str_trim(as.character(x)))
  dplyr::case_when(
    z %in% c("1", "CAR", "CART", "CAR T", "CAR-T", "CCT") ~ "CAR",
    z %in% c("0", "HCT", "AUTO", "AUTOHCT", "AUTO-HCT", "TRANSPLANT") ~ "HCT",
    TRUE ~ as.character(x)
  )
}

# =========================================================
# COMPETING RISK HELPERS
# =========================================================
make_ci_dataset <- function(df,
                            event_type = c("severe", "all"),
                            exclude_covid = FALSE,
                            cap_pre_landmark = 1230,
                            landmark_day = 30,
                            id_col = "CCT_ID") {
  event_type <- match.arg(event_type)

  out <- df %>%
    dplyr::mutate(
      Grade_num        = suppressWarnings(as.numeric(Grade)),
      Days_Infection   = as_num_days(Days_Infection),
      Days_Censor      = as_num_days(Days_Censor),
      Infection_Censor = suppressWarnings(as.integer(Infection_Censor)),
      Censor_Type      = suppressWarnings(as.integer(Censor_Type)),
      Treatment        = factor(Treatment, levels = c("HCT", "CAR"))
    )

  if (exclude_covid && "Type" %in% names(out)) {
    out <- out %>%
      dplyr::mutate(Days_Infection = if_else(Type == "COVID", NA_real_, Days_Infection))
  }

  out <- out %>%
    dplyr::mutate(
      is_event = dplyr::case_when(
        event_type == "severe" ~ !is.na(Grade_num) & Grade_num >= 3 & !is.na(Days_Infection),
        event_type == "all"    ~ !is.na(Days_Infection),
        TRUE ~ FALSE
      ),
      time_raw = dplyr::if_else(is_event, Days_Infection, Days_Censor),
      status_raw = dplyr::case_when(
        is_event & Days_Infection <= Days_Censor ~ 1L,
        TRUE ~ dplyr::coalesce(Censor_Type, 0L)
      ),
      time_cap = pmin(time_raw, cap_pre_landmark, na.rm = TRUE),
      status_cap = if_else(time_raw > cap_pre_landmark, 0L, status_raw)
    )

  out %>%
    dplyr::filter(!is.na(time_cap)) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::slice_min(order_by = time_cap, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::filter(time_cap >= landmark_day) %>%
    dplyr::mutate(
      landmark_day = landmark_day,
      ftime = time_cap - landmark_day,
      fstatus = status_cap,
      event_type_analysis = event_type,
      covid_excluded_as_index = exclude_covid,
      cap_pre_landmark = cap_pre_landmark
    )
}

tidy_cr_stats <- function(dat, label) {
  if (nrow(dat) == 0 || length(unique(dat$Treatment)) < 2) {
    return(tibble::tibble(
      analysis = label,
      statistic = c("Gray test", "Fine-Gray HR CAR vs HCT", "Cause-specific Cox HR CAR vs HCT"),
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p_value = NA_real_,
      p_label = NA_character_
    ))
  }

  ci <- cmprsk::cuminc(ftime = dat$ftime, fstatus = dat$fstatus, group = dat$Treatment, cencode = 0)

  gray_p <- tryCatch({
    as.data.frame(ci$Tests)$pv[1]
  }, error = function(e) NA_real_)

  cov1 <- stats::model.matrix(~ Treatment, dat)[, -1, drop = FALSE]

  fg_row <- tryCatch({
    fg <- cmprsk::crr(ftime = dat$ftime, fstatus = dat$fstatus, cov1 = cov1, failcode = 1, cencode = 0)
    beta <- as.numeric(fg$coef[1])
    se <- sqrt(diag(fg$var))[1]
    p <- 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
    tibble::tibble(
      statistic = "Fine-Gray HR CAR vs HCT",
      estimate = exp(beta),
      conf.low = exp(beta - 1.96 * se),
      conf.high = exp(beta + 1.96 * se),
      p_value = p
    )
  }, error = function(e) {
    tibble::tibble(statistic = "Fine-Gray HR CAR vs HCT", estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p_value = NA_real_)
  })

  cs_row <- tryCatch({
    cs <- survival::coxph(Surv(ftime, as.integer(fstatus == 1)) ~ Treatment, data = dat)
    td <- broom::tidy(cs, exponentiate = TRUE, conf.int = TRUE)
    td %>%
      dplyr::filter(stringr::str_detect(term, "Treatment")) %>%
      dplyr::transmute(
        statistic = "Cause-specific Cox HR CAR vs HCT",
        estimate = estimate,
        conf.low = conf.low,
        conf.high = conf.high,
        p_value = p.value
      )
  }, error = function(e) {
    tibble::tibble(statistic = "Cause-specific Cox HR CAR vs HCT", estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p_value = NA_real_)
  })

  dplyr::bind_rows(
    tibble::tibble(statistic = "Gray test", estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p_value = gray_p),
    fg_row,
    cs_row
  ) %>%
    dplyr::mutate(analysis = label, .before = 1, p_label = fmt_p(p_value))
}

tidy_cuminc_curve <- function(dat, label) {
  if (nrow(dat) == 0 || length(unique(dat$Treatment)) < 1) {
    return(tibble::tibble())
  }

  ci <- cmprsk::cuminc(ftime = dat$ftime, fstatus = dat$fstatus, group = dat$Treatment, cencode = 0)

  purrr::imap_dfr(ci, function(x, nm) {
    if (nm == "Tests" || is.null(x$time) || is.null(x$est)) return(NULL)
    # Keep event-of-interest curve only, i.e. failcode 1.
    if (!stringr::str_detect(nm, " 1$")) return(NULL)
    tibble::tibble(
      analysis = label,
      treatment = stringr::str_remove(nm, " 1$"),
      time_days_post_landmark = x$time,
      cumulative_incidence = x$est,
      variance = x$var
    )
  })
}

plot_cif_curve <- function(curve_df, label, out_base) {
  if (nrow(curve_df) == 0) return(invisible(NULL))

  p <- ggplot(curve_df, aes(x = time_days_post_landmark, y = cumulative_incidence, color = treatment)) +
    geom_step(linewidth = 1.2) +
    scale_color_manual(values = pal_treatment, name = "") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    theme_classic(base_size = 16, base_family = "Arial") +
    xlab("Days post-landmark") +
    ylab("Cumulative incidence") +
    ggtitle(label)

  save_plot_all(p, out_base, width = 7, height = 5)
}

# =========================================================
# PART 1. INPUTS + MATCHED COHORT
# =========================================================
infxn_text <- readr::read_delim(
  PATH_INF,
  delim = "\t",
  quote = "",
  escape_backslash = FALSE,
  escape_double = FALSE,
  trim_ws = TRUE,
  locale = readr::locale(encoding = "windows-1252"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  progress = FALSE
)

match_raw <- readr::read_tsv(PATH_MATCH, show_col_types = FALSE, guess_max = 100000, progress = FALSE) %>%
  janitor::clean_names()

required_match_cols <- c("cct", "age", "sex_m", "lines", "time_from_infusion", "cct_id")
missing_match <- setdiff(required_match_cols, names(match_raw))
if (length(missing_match) > 0) {
  stop("The match sheet is missing required columns after clean_names(): ", paste(missing_match, collapse = ", "))
}

match <- match_raw %>%
  dplyr::mutate(
    age = as.numeric(age),
    lines = as.numeric(lines),
    time_from_infusion = as.numeric(time_from_infusion),
    treatment_group = coerce_treatment(cct)
  )

vars_tbl1 <- intersect(c("age", "sex_m", "lines", "time_from_infusion", "pre_wbc", "pre_anc", "pre_alc", "pre_hg", "pre_plt", "auto_regimen"), names(match))

m_mod <- MatchIt::matchit(
  cct ~ age + sex_m + lines + time_from_infusion,
  data = match,
  method = "nearest",
  caliper = 0.25
)

dta_m <- MatchIt::match.data(m_mod)
keep_ids <- as.character(dta_m$cct_id)

id_map <- make_id_map(keep_ids, prefix = "InfectionSubject")

# Create optional columns up front so transmute() is stable across source-file versions.
for (nm in c("subclass", "weights", "pre_wbc", "pre_anc", "pre_alc", "pre_hg", "pre_plt", "auto_regimen")) {
  if (!nm %in% names(dta_m)) dta_m[[nm]] <- NA
}

matched_covariates <- dta_m %>%
  dplyr::mutate(
    figure_patient_id = lookup_public_id(cct_id, id_map, prefix = "InfectionSubject"),
    treatment_group = coerce_treatment(cct)
  ) %>%
  dplyr::transmute(
    figure_patient_id,
    treatment_group,
    matched_subclass = as.character(subclass),
    matching_weight = suppressWarnings(as.numeric(weights)),
    age,
    sex_m,
    lines,
    time_from_infusion,
    pre_wbc,
    pre_anc,
    pre_alc,
    pre_hg,
    pre_plt,
    auto_regimen = as.character(auto_regimen)
  )

write_figshare_csv(
  matched_covariates,
  "figshare_infection_matched_baseline_covariates.csv",
  "Matched CAR versus autoHCT cohort baseline covariates used for the infection analyses."
)

balance_tbl <- tryCatch({
  bal <- cobalt::bal.tab(m_mod, un = TRUE)
  as.data.frame(bal$Balance) %>%
    tibble::rownames_to_column("covariate") %>%
    dplyr::rename_with(~str_replace_all(.x, "[^A-Za-z0-9]+", "_"))
}, error = function(e) {
  tibble::tibble(covariate = character(), note = character())
})

write_figshare_csv(
  balance_tbl,
  "figshare_infection_matching_balance_statistics.csv",
  "Standardized mean-difference balance statistics from the propensity-matched cohort."
)

# Save love plot for record; not included in manuscript upload if not needed.
try({
  p_love <- cobalt::love.plot(cobalt::bal.tab(m_mod, un = TRUE), abs = TRUE, thresholds = c(m = 0.1)) +
    theme_classic(base_size = 14, base_family = "Arial")
  save_plot_all(p_love, "infection_matching_love_plot", width = 7, height = 5)
}, silent = TRUE)

# Restrict infection records to matched cohort.
if (!"CCT_ID" %in% names(infxn_text)) stop("Infection text is missing CCT_ID.")
keep_match_infxn <- infxn_text %>%
  dplyr::filter(CCT_ID %in% keep_ids) %>%
  dplyr::mutate(
    Treatment = factor(Treatment, levels = c("HCT", "CAR")),
    figure_patient_id = lookup_public_id(CCT_ID, id_map, prefix = "InfectionSubject")
  )

for (nm in c("Type", "Grade", "Days_Infection", "Days_Censor", "Censor_Type")) {
  if (!nm %in% names(keep_match_infxn)) keep_match_infxn[[nm]] <- NA
}

infection_event_rows <- keep_match_infxn %>%
  dplyr::transmute(
    figure_patient_id,
    Treatment,
    infection_type = Type,
    infection_grade = suppressWarnings(as.numeric(Grade)),
    days_infection = as_num_days(Days_Infection),
    days_censor = as_num_days(Days_Censor),
    censor_type = suppressWarnings(as.integer(Censor_Type))
  )

write_figshare_csv(
  infection_event_rows,
  "figshare_infection_all_matched_event_rows.csv",
  "Deidentified infection event/censor rows for the matched cohort prior to landmark/event filtering."
)

# =========================================================
# PART 2. COMPETING RISKS: INFECTION ENDPOINTS
# =========================================================
ci_specs <- tibble::tribble(
  ~analysis_id,              ~event_type, ~exclude_covid, ~cap_pre_landmark, ~landmark_day, ~title,
  "severe_postD30",          "severe",    FALSE,          1230,              30,            "Severe infection (Grade >=3), post-Day 30",
  "all_postD30",             "all",       FALSE,          1230,              30,            "All infections, post-Day 30",
  "severe_noCOVID_postD30",  "severe",    TRUE,           1230,              30,            "Severe infection (Grade >=3), no COVID as index, post-Day 30",
  "severe_postD365",         "severe",    FALSE,          1500,              365,           "Severe infection (Grade >=3), post-Day 365"
)

ci_patient_all <- purrr::pmap_dfr(ci_specs, function(analysis_id, event_type, exclude_covid, cap_pre_landmark, landmark_day, title) {
  dat <- make_ci_dataset(
    keep_match_infxn,
    event_type = event_type,
    exclude_covid = exclude_covid,
    cap_pre_landmark = cap_pre_landmark,
    landmark_day = landmark_day,
    id_col = "CCT_ID"
  ) %>%
    dplyr::mutate(
      analysis_id = analysis_id,
      analysis_title = title,
      figure_patient_id = lookup_public_id(CCT_ID, id_map, prefix = "InfectionSubject")
    )

  dat %>%
    dplyr::transmute(
      analysis_id,
      analysis_title,
      figure_patient_id,
      Treatment,
      landmark_day,
      event_type_analysis,
      covid_excluded_as_index,
      cap_pre_landmark,
      ftime,
      fstatus,
      event_of_interest = as.integer(fstatus == 1),
      infection_type = if ("Type" %in% names(.)) Type else NA_character_,
      infection_grade = if ("Grade_num" %in% names(.)) Grade_num else NA_real_,
      days_infection = Days_Infection,
      days_censor = Days_Censor
    )
})

write_figshare_csv(
  ci_patient_all,
  "figshare_infection_competing_risk_patient_level_data.csv",
  "Patient-level landmarked competing-risk datasets for severe/all infection endpoints."
)

ci_stats_all <- purrr::pmap_dfr(ci_specs, function(analysis_id, event_type, exclude_covid, cap_pre_landmark, landmark_day, title) {
  dat <- make_ci_dataset(
    keep_match_infxn,
    event_type = event_type,
    exclude_covid = exclude_covid,
    cap_pre_landmark = cap_pre_landmark,
    landmark_day = landmark_day,
    id_col = "CCT_ID"
  )
  tidy_cr_stats(dat, analysis_id)
})

write_figshare_csv(
  ci_stats_all,
  "figshare_infection_competing_risk_statistics.csv",
  "Gray test, Fine-Gray, and cause-specific Cox statistics for infection competing-risk analyses."
)

ci_curves_all <- purrr::pmap_dfr(ci_specs, function(analysis_id, event_type, exclude_covid, cap_pre_landmark, landmark_day, title) {
  dat <- make_ci_dataset(
    keep_match_infxn,
    event_type = event_type,
    exclude_covid = exclude_covid,
    cap_pre_landmark = cap_pre_landmark,
    landmark_day = landmark_day,
    id_col = "CCT_ID"
  )
  tidy_cuminc_curve(dat, analysis_id)
})

write_figshare_csv(
  ci_curves_all,
  "figshare_infection_competing_risk_cif_curve_data.csv",
  "Step-function cumulative-incidence curve coordinates generated from cmprsk::cuminc for infection endpoints."
)

purrr::pwalk(ci_specs, function(analysis_id, event_type, exclude_covid, cap_pre_landmark, landmark_day, title) {
  curve_df <- ci_curves_all %>% dplyr::filter(analysis == analysis_id)
  plot_cif_curve(curve_df, title, paste0("infection_CIF_", analysis_id))
})

# =========================================================
# PART 3. POST-DAY 30 OVERALL SURVIVAL / NON-RELAPSE SURVIVAL
# =========================================================
if (all(c("OS_Days", "OS_Event") %in% names(keep_match_infxn))) {
  df_os <- keep_match_infxn %>%
    dplyr::mutate(
      OS_Days = as_num_days(OS_Days),
      OS_Event = suppressWarnings(as.integer(OS_Event))
    ) %>%
    dplyr::group_by(CCT_ID) %>%
    dplyr::slice_min(order_by = OS_Days, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::filter(OS_Days >= 30) %>%
    dplyr::mutate(
      OS_time = OS_Days - 30,
      figure_patient_id = lookup_public_id(CCT_ID, id_map, prefix = "InfectionSubject")
    )

  os_figshare <- df_os %>%
    dplyr::transmute(
      figure_patient_id,
      Treatment,
      OS_time_post_day30 = OS_time,
      OS_Days,
      OS_Event
    )

  write_figshare_csv(
    os_figshare,
    "figshare_infection_postD30_survival_plot_data.csv",
    "Patient-level post-Day 30 survival data for matched CAR versus autoHCT cohort."
  )

  os_stats <- tryCatch({
    sf <- survival::survdiff(Surv(OS_time, OS_Event) ~ Treatment, data = df_os)
    p <- stats::pchisq(sf$chisq, df = length(sf$n) - 1, lower.tail = FALSE)
    tibble::tibble(
      analysis = "postD30_survival",
      statistic = "log-rank",
      p_value = p,
      p_label = fmt_p(p)
    )
  }, error = function(e) tibble::tibble(analysis = "postD30_survival", statistic = "log-rank", p_value = NA_real_, p_label = NA_character_))

  write_figshare_csv(
    os_stats,
    "figshare_infection_postD30_survival_statistics.csv",
    "Log-rank statistics for post-Day 30 survival."
  )

  p_os <- survminer::ggsurvplot(
    fit = survival::survfit(Surv(OS_time, OS_Event) ~ Treatment, data = df_os),
    xlab = "Days Post Infusion",
    ylab = "Non-Relapse Survival",
    pval = TRUE, conf.int = FALSE, surv.median.line = "hv",
    palette = pal_treatment,
    risk.table = TRUE, break.time.by = 200, xlim = c(0, 1400),
    risk.table.height = 0.20, risk.table.y.text = FALSE,
    ggtheme = survminer::theme_survminer(text = element_text(size = 14))
  )
  p_os$plot <- p_os$plot +
    scale_x_continuous(
      breaks = c(0, 170, 370, 570, 770, 970, 1170, 1370),
      labels = c(30, 200, 400, 600, 800, 1000, 1200, 1400)
    )

  ggsave(out_file("infection_postD30_survival.pdf"), plot = p_os$plot, width = 7, height = 5, device = cairo_pdf)
  ggsave(out_file("infection_postD30_survival.png"), plot = p_os$plot, width = 7, height = 5, dpi = 400)
}

# =========================================================
# PART 4. INFECTION CATEGORY / GRADE BAR PLOT
# =========================================================
if (all(c("Type", "Grade", "Treatment") %in% names(keep_match_infxn))) {
  infxn_text_noCMV <- keep_match_infxn %>% dplyr::filter(Type != "CMV Reactivation")

  infection_bar_data <- infxn_text_noCMV %>%
    dplyr::count(Treatment, Type, Grade, name = "n_infections") %>%
    dplyr::group_by(Treatment) %>%
    dplyr::mutate(
      pct_within_treatment = n_infections / sum(n_infections),
      pct_within_treatment_percent = 100 * pct_within_treatment
    ) %>%
    dplyr::ungroup()

  infection_type_order <- infection_bar_data %>%
    dplyr::group_by(Type) %>%
    dplyr::summarise(total_pct = sum(pct_within_treatment), .groups = "drop") %>%
    dplyr::arrange(total_pct)

  infection_bar_data <- infection_bar_data %>%
    dplyr::mutate(
      Type = factor(Type, levels = infection_type_order$Type),
      Grade = factor(Grade, levels = sort(unique(Grade)))
    )

  write_figshare_csv(
    infection_bar_data %>% dplyr::rename(infection_type = Type, infection_grade = Grade),
    "figshare_infection_category_grade_bar_plot_data.csv",
    "Infection type and grade frequencies used for treatment-faceted stacked bar plot."
  )

  grade_cols <- c("#4DAF4A", "#377EB8", "#984EA3", "#FF7F00", "#E41A1C")
  names(grade_cols) <- levels(infection_bar_data$Grade)[seq_along(levels(infection_bar_data$Grade))]

  p_bar <- ggplot(infection_bar_data, aes(x = pct_within_treatment_percent, y = Type, fill = Grade)) +
    geom_vline(xintercept = c(5, 10, 15, 20), linetype = "dotted") +
    geom_col(position = "stack") +
    scale_fill_manual(name = "Grade", values = grade_cols, drop = FALSE) +
    facet_grid(~Treatment) +
    coord_cartesian(xlim = c(0, 20)) +
    ylab("Category of Infection") +
    xlab("Percent of Infections by Treatment") +
    theme_classic(base_size = 14, base_family = "Arial")

  save_plot_all(p_bar, "infection_category_grade_bar_plot", width = 9, height = 7)
}

# =========================================================
# PART 5. MATCHED CYTOPENIA D28-D800 LMM DATASETS
# =========================================================
if (file.exists(PATH_CYTO)) {
  match_cyto <- readr::read_delim(PATH_CYTO, show_col_types = FALSE, guess_max = 100000, progress = FALSE)

  if (!"Systemic_Lines" %in% names(match_cyto)) {
    if ("Prior_Lines" %in% names(match_cyto)) match_cyto <- dplyr::rename(match_cyto, Systemic_Lines = Prior_Lines)
    if ("Lines" %in% names(match_cyto)) match_cyto <- dplyr::rename(match_cyto, Systemic_Lines = Lines)
  }

  cyto_id_col <- first_existing(c("PID", "CCT_ID", "cct_id", "ID", "Patient_ID", "LYMID", "AAID"), names(match_cyto))
  if (is.na(cyto_id_col)) {
    match_cyto$raw_cyto_id <- seq_len(nrow(match_cyto))
    cyto_id_col <- "raw_cyto_id"
  }

  # Filter to the matched cohort only when the cytopenia ID column overlaps with match CCT IDs.
  if (any(as.character(match_cyto[[cyto_id_col]]) %in% keep_ids)) {
    match_cyto <- match_cyto %>% dplyr::filter(as.character(.data[[cyto_id_col]]) %in% keep_ids)
  } else {
    message("Cytopenia ID column did not overlap with matched CCT IDs; retaining all cytopenia rows and assigning figure IDs from ", cyto_id_col, ".")
  }

  cyto_id_map <- make_id_map(match_cyto[[cyto_id_col]], prefix = "CytoSubject")
  day_levels_post <- c("28", "60", "90", "180", "365", "500", "800")
  outcomes <- c("WBC", "ANC", "ALC", "HG", "PLT")
  outcomes <- outcomes[outcomes %in% names(match_cyto)]

  cyto_match_post <- match_cyto %>%
    dplyr::filter(!as.character(Day) %in% c("-6", "0", "7", "14", "21")) %>%
    dplyr::mutate(
      figure_patient_id = lookup_public_id(.data[[cyto_id_col]], cyto_id_map, prefix = "CytoSubject"),
      random_effect_id = factor(as.character(.data[[cyto_id_col]])),
      Day = factor(as.character(Day), levels = day_levels_post),
      CCT_Type = factor(coerce_treatment(CCT_Type), levels = c("HCT", "CAR")),
      Sex = factor(as.character(Sex)),
      Sex = if ("M" %in% levels(Sex)) stats::relevel(Sex, ref = "M") else Sex,
      Age = suppressWarnings(as.numeric(Age)),
      Systemic_Lines = suppressWarnings(as.numeric(Systemic_Lines))
    ) %>%
    dplyr::filter(!is.na(Day), !is.na(CCT_Type))

  cyto_figshare_long <- cyto_match_post %>%
    dplyr::select(figure_patient_id, Day, CCT_Type, Age, Sex, Systemic_Lines, dplyr::all_of(outcomes)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(outcomes), names_to = "outcome", values_to = "value") %>%
    dplyr::filter(!is.na(value))

  write_figshare_csv(
    cyto_figshare_long,
    "figshare_matched_cytopenia_longitudinal_values_D28toD800.csv",
    "Long-format matched cytopenia values used for D28-D800 median plots and LMMs."
  )

  cyto_summary <- cyto_figshare_long %>%
    dplyr::group_by(outcome, Day, CCT_Type) %>%
    dplyr::summarise(
      median = median(value, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      n = dplyr::n(),
      se = sd / sqrt(n),
      .groups = "drop"
    )

  write_figshare_csv(
    cyto_summary,
    "figshare_matched_cytopenia_median_summary_D28toD800.csv",
    "Median/SE cytopenia summaries by treatment group and day."
  )

  fit_lmm <- function(y) {
    dat <- cyto_match_post %>% dplyr::filter(!is.na(.data[[y]]))
    fml <- stats::as.formula(paste0(y, " ~ Day * (CCT_Type + Age + Sex + Systemic_Lines) + (1|random_effect_id)"))
    lmerTest::lmer(fml, data = dat, REML = TRUE)
  }

  models_matched <- purrr::set_names(purrr::map(outcomes, ~tryCatch(fit_lmm(.x), error = function(e) e)), outcomes)

  type3_stats <- purrr::imap_dfr(models_matched, function(mod, outcome) {
    if (inherits(mod, "error")) {
      return(tibble::tibble(outcome = outcome, term = NA_character_, statistic = NA_real_, p_value = NA_real_, model_error = mod$message))
    }
    a <- as.data.frame(stats::anova(mod, type = 3))
    a %>%
      tibble::rownames_to_column("term") %>%
      dplyr::transmute(
        outcome = outcome,
        term,
        numDF = if ("NumDF" %in% names(.)) NumDF else NA_real_,
        denDF = if ("DenDF" %in% names(.)) DenDF else NA_real_,
        statistic = if ("F.value" %in% names(.)) F.value else if ("F value" %in% names(.)) `F value` else NA_real_,
        p_value = if ("Pr(>F)" %in% names(.)) `Pr(>F)` else NA_real_,
        p_label = fmt_p(p_value),
        model_error = NA_character_
      )
  })

  write_figshare_csv(
    type3_stats,
    "figshare_matched_cytopenia_LMM_type3_statistics_D28toD800.csv",
    "Type III ANOVA statistics for matched cytopenia LMMs."
  )

  emm_byday <- purrr::imap_dfr(models_matched, function(mod, outcome) {
    if (inherits(mod, "error")) return(NULL)
    out <- tryCatch({
      emmeans::contrast(
        emmeans::emmeans(mod, ~ CCT_Type | Day),
        method = list("CAR - HCT" = c(-1, +1)),
        adjust = "holm"
      ) %>%
        summary(infer = TRUE) %>%
        as.data.frame() %>%
        dplyr::transmute(
          outcome = outcome,
          Day = as.character(Day),
          contrast,
          estimate,
          SE,
          df,
          lower.CL,
          upper.CL,
          p_value = p.value,
          p_label = fmt_p(p_value)
        )
    }, error = function(e) NULL)
    out
  })

  write_figshare_csv(
    emm_byday,
    "figshare_matched_cytopenia_LMM_daywise_contrasts_D28toD800.csv",
    "Day-wise adjusted CAR-HCT contrasts from matched cytopenia LMMs, Holm-adjusted."
  )

  emm_avg <- purrr::imap_dfr(models_matched, function(mod, outcome) {
    if (inherits(mod, "error")) return(NULL)
    out <- tryCatch({
      emmeans::contrast(
        emmeans::emmeans(mod, ~ CCT_Type, weights = "equal"),
        method = list("CAR - HCT" = c(-1, +1)),
        adjust = "none"
      ) %>%
        summary(infer = TRUE) %>%
        as.data.frame() %>%
        dplyr::transmute(
          outcome = outcome,
          contrast,
          estimate,
          SE,
          df,
          lower.CL,
          upper.CL,
          p_value = p.value,
          p_label = fmt_p(p_value)
        )
    }, error = function(e) NULL)
    out
  })

  write_figshare_csv(
    emm_avg,
    "figshare_matched_cytopenia_LMM_averaged_contrasts_D28toD800.csv",
    "Equal-weight averaged CAR-HCT contrasts across D28-D800 from matched cytopenia LMMs."
  )

  emm_adjusted_means <- purrr::imap_dfr(models_matched, function(mod, outcome) {
    if (inherits(mod, "error")) return(NULL)
    tryCatch({
      as.data.frame(emmeans::emmeans(mod, ~ CCT_Type | Day)) %>%
        dplyr::transmute(
          outcome = outcome,
          Day = as.character(Day),
          CCT_Type,
          emmean,
          SE,
          df,
          lower.CL,
          upper.CL
        )
    }, error = function(e) NULL)
  })

  write_figshare_csv(
    emm_adjusted_means,
    "figshare_matched_cytopenia_LMM_adjusted_mean_plot_data_D28toD800.csv",
    "Adjusted mean estimates used for matched cytopenia EMM plots."
  )

  # Optional time-weighted estimates from the by-day adjusted contrast table.
  time_weights <- tibble::tibble(
    Day = day_levels_post,
    day_num = as.numeric(day_levels_post)
  ) %>%
    dplyr::arrange(day_num) %>%
    dplyr::mutate(
      weight = dplyr::case_when(
        dplyr::row_number() == 1 ~ (dplyr::lead(day_num) - day_num) / 2,
        dplyr::row_number() == dplyr::n() ~ (day_num - dplyr::lag(day_num)) / 2,
        TRUE ~ ((day_num - dplyr::lag(day_num)) + (dplyr::lead(day_num) - day_num)) / 2
      )
    )

  emm_timeweighted <- emm_byday %>%
    dplyr::left_join(time_weights, by = "Day") %>%
    dplyr::group_by(outcome, contrast) %>%
    dplyr::summarise(
      weighted_average_estimate = sum(estimate * weight, na.rm = TRUE) / sum(weight[!is.na(estimate)], na.rm = TRUE),
      AUC_difference = sum(estimate * weight, na.rm = TRUE),
      total_weight_days = sum(weight[!is.na(estimate)], na.rm = TRUE),
      .groups = "drop"
    )

  write_figshare_csv(
    emm_timeweighted,
    "figshare_matched_cytopenia_LMM_timeweighted_estimates_D28toD800.csv",
    "Descriptive time-weighted CAR-HCT differences calculated from day-wise LMM contrasts."
  )

  # Plot adjusted means and raw medians for each endpoint.
  purrr::walk(outcomes, function(y) {
    adj_dat <- emm_adjusted_means %>% dplyr::filter(outcome == y)
    if (nrow(adj_dat) > 0) {
      p_adj <- ggplot(adj_dat, aes(x = Day, y = emmean, color = CCT_Type, group = CCT_Type)) +
        geom_point(position = position_dodge(width = 0.12)) +
        geom_line(position = position_dodge(width = 0.12), linewidth = 1) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = position_dodge(width = 0.12)) +
        scale_color_manual(values = pal_treatment) +
        theme_classic(base_size = 14, base_family = "Arial") +
        xlab("Day") + ylab(paste0("Adjusted ", y))
      save_plot_all(p_adj, paste0("matched_cytopenia_adjusted_mean_", y, "_D28toD800"), width = 7, height = 5)
    }

    med_dat <- cyto_summary %>% dplyr::filter(outcome == y)
    if (nrow(med_dat) > 0) {
      p_med <- ggplot(med_dat, aes(x = Day, y = median, color = CCT_Type, group = CCT_Type)) +
        geom_point() +
        geom_line(linewidth = 1.1) +
        geom_errorbar(aes(ymin = median - se, ymax = median + se), width = 0.1) +
        scale_color_manual(values = pal_treatment) +
        theme_classic(base_size = 14, base_family = "Arial") +
        xlab("Day") + ylab(paste0("Median ", y))
      save_plot_all(p_med, paste0("matched_cytopenia_median_", y, "_D28toD800"), width = 7, height = 5)
    }
  })
} else {
  message("Cytopenia input file not found; skipping matched cytopenia LMM exports.")
}

# =========================================================
# MANIFEST + README
# =========================================================
readr::write_csv(manifest, figshare_file("figshare_infection_CAR_vs_autoHCT_dataset_manifest.csv"), na = "")

readme <- c(
  "CAR vs autoHCT infection and cytopenia Figshare datasets",
  "=======================================================",
  "",
  "This directory contains figure-specific datasets generated from the CAR T-cell versus autologous HCT infection analysis script.",
  "LYMID/AAID-style identifiers are retained if present. CCT IDs, S IDs, MRNs, SID, Sample_ID, Trial_ID, and internal sample columns are not exported.",
  "When the source data only contained CCT IDs, the script generated figure_patient_id labels that preserve grouping but are not keys back to clinical identifiers.",
  "",
  "Competing-risk datasets use landmarking as specified in the analysis_id fields.",
  "For fstatus: 0 = censored, 1 = infection event of interest, and other values retain competing-event codes from the source censor field.",
  "",
  "Datasets include only rows retained by the analysis filters used for figure generation. Values excluded in code are not retained in the final upload files.",
  "",
  "Generated files are listed in figshare_infection_CAR_vs_autoHCT_dataset_manifest.csv."
)

writeLines(readme, con = figshare_file("figshare_infection_CAR_vs_autoHCT_README.txt"))
writeLines(capture.output(sessionInfo()), con = out_file("infection_CAR_vs_autoHCT_figshare_sessionInfo.txt"))

cat("\nSaved infection/CAR-vs-autoHCT figure datasets to:\n", figshare_dir, "\n", sep = "")
cat("Saved figure outputs to:\n", out_dir, "\n", sep = "")
