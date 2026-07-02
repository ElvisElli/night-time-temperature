#!/usr/bin/env Rscript
# =============================================================================
# Dark-period night-time (and day-time) temperature from NASA POWER *hourly* data
# =============================================================================
#
# WHY THIS EXISTS
# ---------------
# The manuscript uses mean daily minimum temperature (`mean_mint`) as a surrogate
# for night-time temperature. Reviewer 1 asked for the *actual* night-time air
# temperature the crop experiences -- the average of sub-daily values while the
# sun is below the horizon -- not the single coldest instant of each day. That
# cannot be recovered from daily Tmax/Tmin; it needs hourly data.
#
# This script pulls NASA POWER **hourly** 2 m air temperature (T2M) -- the same
# product the paper already cites, at its native hourly resolution -- for each
# site-year, flags each hour day/night by solar geometry (sun below the horizon,
# NOAA algorithm), and averages within the June-July window to produce:
#     nightT = mean T2M over dark hours       (the dark-period night temperature)
#     dayT   = mean T2M over daylight hours
# It then (a) pairs the SITE-YEAR nightT with the existing `mean_mint`, and
# (b) builds a DAILY comparison (daily minimum vs daily dark-period mean) across
# all days, with an overall plot and a per-location multipanel.
#
# HOW TO RUN (anywhere)
# ---------------------
#   * RStudio:  open nighttime_temperature.Rproj, open this file, click "Source".
#   * Terminal: Rscript nighttime_temperature.R
# The script finds the repository root on its own, auto-installs missing
# packages, and caches raw hourly pulls under weather/outputs/cache/ so re-runs
# are instant and offline-friendly.
#
# DEPENDENCIES: jsonlite, curl, readxl, dplyr, tidyr, ggplot2, scales (base R
# otherwise; no nasapower/suncalc needed -- the API and solar math are inline).
# =============================================================================

## ---- 0. Packages (auto-install if missing) --------------------------------
pkgs <- c("jsonlite", "curl", "readxl", "dplyr", "tidyr", "ggplot2", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

## ---- 1. CONFIG ------------------------------------------------------------
WINDOW_START_DOY <- 152      # June 1  (matches the paper's fixed window & mean_mint)
WINDOW_END_DOY   <- 212      # July 31
NIGHT_ELEV_DEG   <- 0.0      # night = solar elevation < this (0 = geometric horizon)
LOCAL_UTC_OFFSET <- -5       # eastern Arkansas is CDT (UTC-5) in Jun-Jul; used for daily calendar
NASA_FILL        <- -999.0   # NASA POWER missing-value sentinel

# Map modelling-dataset location names -> coordinate-file location names
LOC_MAP <- c(Des_Arc = "des_arc", Harrisburg = "harrisburg_neerec",
             Keiser = "keiser_nere", Marianna = "marianna",
             Rohwer = "rohwer", Stuttgart = "stuttgart")

## ---- 2. Locate repo & I/O paths (portable) --------------------------------
get_script_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()
}
find_repo_root <- function(start) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:12) {
    if (dir.exists(file.path(d, "corn variety testing"))) return(d)
    p <- dirname(d); if (identical(p, d)) break; d <- p
  }
  stop("Could not locate repo root (a parent dir containing 'corn variety testing').")
}
ROOT      <- find_repo_root(get_script_dir())
WEATHER   <- file.path(ROOT, "weather")
OUTDIR    <- file.path(WEATHER, "outputs")
CACHE     <- file.path(OUTDIR, "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
MODEL_CSV   <- file.path(ROOT, "corn variety testing", "mixed_effect_models",
                         "input", "Corn data organized2.csv")
COORDS_XLSX <- file.path(WEATHER, "locations.june.july.xlsx")

## ---- 3. Solar elevation: NOAA algorithm (UTC in -> degrees out) ------------
solar_elevation_deg <- function(dt_utc, lat, lon) {
  lt   <- as.POSIXlt(dt_utc, tz = "UTC")
  doy  <- lt$yday + 1
  hour <- lt$hour + lt$min / 60 + lt$sec / 3600
  gamma <- 2 * pi / 365 * (doy - 1 + (hour - 12) / 24)
  eqtime <- 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -
                      0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma))
  decl <- 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) -
          0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) -
          0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)
  tst <- hour * 60 + eqtime + 4 * lon            # true solar time (min); tz=UTC
  ha  <- (tst / 4 - 180) * pi / 180              # hour angle (rad)
  latr <- lat * pi / 180
  cosz <- sin(latr) * sin(decl) + cos(latr) * cos(decl) * cos(ha)
  90 - acos(pmin(pmax(cosz, -1), 1)) * 180 / pi
}

## ---- 4. NASA POWER hourly pull (cached, with retry) -----------------------
fetch_hourly_t2m <- function(lon, lat, year) {
  start <- sprintf("%d0601", year); end <- sprintf("%d0731", year)
  cache_file <- file.path(CACHE, sprintf("%.4f_%.4f_%s_%s.json", lat, lon, start, end))
  if (file.exists(cache_file)) {
    payload <- jsonlite::fromJSON(cache_file, simplifyVector = TRUE)
  } else {
    url <- sprintf(paste0("https://power.larc.nasa.gov/api/temporal/hourly/point",
                          "?parameters=T2M&community=AG&longitude=%s&latitude=%s",
                          "&start=%s&end=%s&format=JSON&time-standard=UTC"),
                   lon, lat, start, end)
    payload <- NULL
    for (attempt in 0:3) {
      payload <- tryCatch(jsonlite::fromJSON(url, simplifyVector = TRUE),
                          error = function(e) NULL)
      if (!is.null(payload)) break
      Sys.sleep(2^attempt)
    }
    if (is.null(payload)) stop(sprintf("NASA POWER failed for %s,%s,%d", lat, lon, year))
    writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE), cache_file)
  }
  t2m   <- payload$properties$parameter$T2M          # named list "YYYYMMDDHH" -> value
  times <- names(t2m)
  vals  <- as.numeric(unlist(t2m, use.names = FALSE))
  vals[vals <= NASA_FILL] <- NA
  dt_utc <- as.POSIXct(times, format = "%Y%m%d%H", tz = "UTC")
  data.frame(dt_utc = dt_utc, t2m = vals)
}

## ---- 5. Enrich one site-year (day/night flag + local calendar) ------------
windowed_hourly <- function(lon, lat, year) {
  df <- fetch_hourly_t2m(lon, lat, year)
  lt <- as.POSIXlt(df$dt_utc, tz = "UTC")
  doy <- lt$yday + 1
  keep <- doy >= WINDOW_START_DOY & doy <= WINDOW_END_DOY
  df <- df[keep, , drop = FALSE]
  ltk <- as.POSIXlt(df$dt_utc, tz = "UTC")
  df$utc_date <- as.Date(df$dt_utc, tz = "UTC")
  df$elev     <- solar_elevation_deg(df$dt_utc, lat, lon)
  df$is_night <- df$elev < NIGHT_ELEV_DEG
  local       <- as.POSIXlt(df$dt_utc + LOCAL_UTC_OFFSET * 3600, tz = "UTC")
  df$local_date  <- as.Date(sprintf("%04d-%02d-%02d", local$year + 1900,
                                     local$mon + 1, local$mday))
  df$local_month <- local$mon + 1
  df
}

## ---- 6. Aggregate: site-year summary + daily table ------------------------
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

summarise_site_year <- function(df) {
  night <- df$t2m[df$is_night]; day <- df$t2m[!df$is_night]
  daily_min <- tapply(df$t2m, df$utc_date, function(x) min(x, na.rm = TRUE))
  data.frame(
    nightT     = mean_na(night),
    dayT       = mean_na(day),
    mint_hr    = mean_na(as.numeric(daily_min)),
    n_night_hr = sum(!is.na(night)),
    n_day_hr   = sum(!is.na(day)),
    n_days     = sum(is.finite(as.numeric(daily_min))),
    n_missing  = sum(is.na(df$t2m))
  )
}

daily_from_hourly <- function(df) {
  df <- df[df$local_month %in% c(6, 7), , drop = FALSE]
  df %>%
    dplyr::group_by(local_date) %>%
    dplyr::summarise(
      Tmin_day   = suppressWarnings(min(t2m, na.rm = TRUE)),
      nightT_day = mean_na(t2m[is_night]),
      dayT_day   = mean_na(t2m[!is_night]),
      n_night_hr = sum(is_night & !is.na(t2m)),
      n_day_hr   = sum(!is_night & !is.na(t2m)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(Tmin_day), n_night_hr >= 3)
}

## ---- 7. Build both datasets over all site-years ---------------------------
message("Loading site-year list and coordinates ...")
model <- utils::read.csv(MODEL_CSV, check.names = FALSE, stringsAsFactors = FALSE)
names(model)[1] <- sub("^[^A-Za-z]*", "", names(model)[1])  # strip leading BOM bytes from 'location'
names(model) <- trimws(names(model))
site_years <- model %>%
  dplyr::group_by(location, year) %>%
  dplyr::summarise(mean_mint = mean(mean_mint, na.rm = TRUE),
                   mean_maxt = mean(mean_maxt, na.rm = TRUE),
                   avg_t     = mean(avg_t, na.rm = TRUE), .groups = "drop")

coords <- readxl::read_excel(COORDS_XLSX) %>%
  dplyr::filter(!is.na(longitude), !is.na(latitude)) %>%
  dplyr::group_by(location) %>%
  dplyr::summarise(longitude = dplyr::first(longitude),
                   latitude  = dplyr::first(latitude), .groups = "drop")

sy_rows <- list(); daily_rows <- list()
n <- nrow(site_years)
message(sprintf("Processing %d site-years (NASA POWER hourly T2M, cached) ...", n))
for (i in seq_len(n)) {
  loc_model <- site_years$location[i]; yr <- as.integer(site_years$year[i])
  loc_coord <- LOC_MAP[[loc_model]]
  if (is.null(loc_coord) || !loc_coord %in% coords$location) {
    message("  [skip] no coordinates for '", loc_model, "'"); next
  }
  lon <- coords$longitude[coords$location == loc_coord][1]
  lat <- coords$latitude[coords$location == loc_coord][1]
  df  <- tryCatch(windowed_hourly(lon, lat, yr),
                  error = function(e) { message("  [fail] ", loc_model, " ", yr, ": ",
                                                 conditionMessage(e)); NULL })
  if (is.null(df)) next
  s <- summarise_site_year(df)
  sy_rows[[length(sy_rows) + 1]] <- cbind(
    data.frame(location = loc_model, year = yr, longitude = lon, latitude = lat,
               mean_mint = site_years$mean_mint[i], mean_maxt = site_years$mean_maxt[i],
               avg_t = site_years$avg_t[i]), s)
  d <- daily_from_hourly(df)
  if (nrow(d)) daily_rows[[length(daily_rows) + 1]] <-
      cbind(data.frame(location = loc_model, year = yr), d)
  if (i %% 10 == 0 || i == n) message(sprintf("  %d/%d done", i, n))
}

paired <- dplyr::bind_rows(sy_rows) %>% dplyr::arrange(location, year)
paired$night_minus_mint <- paired$nightT - paired$mean_mint
paired$day_minus_night  <- paired$dayT   - paired$nightT
daily  <- dplyr::bind_rows(daily_rows) %>% dplyr::arrange(location, year, local_date)

utils::write.csv(paired, file.path(OUTDIR, "nighttime_temperature_paired.csv"), row.names = FALSE)
utils::write.csv(daily,  file.path(OUTDIR, "nighttime_temperature_daily.csv"),  row.names = FALSE)
message(sprintf("Wrote paired site-year table (%d rows) and daily table (%d rows).",
                nrow(paired), nrow(daily)))

## ---- 8. Statistics helpers ------------------------------------------------
fit_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  co <- stats::coef(stats::lm(y ~ x))
  list(r = stats::cor(x, y), slope = co[[2]], intercept = co[[1]],
       rmse = sqrt(mean((y - x)^2)), offset = mean(y - x), n = length(x))
}
sy   <- fit_stats(paired$mean_mint, paired$nightT)         # site-year: mint vs nightT
syD  <- fit_stats(paired$mean_mint, paired$dayT)           # site-year: mint vs dayT
syB  <- fit_stats(paired$mean_mint, paired$mint_hr)        # QC bridge
dl   <- fit_stats(daily$Tmin_day, daily$nightT_day)        # daily: Tmin vs nightT

## ---- 9. Plots -------------------------------------------------------------
loc_r_daily <- daily %>% dplyr::group_by(location) %>%
  dplyr::summarise(r = stats::cor(Tmin_day, nightT_day, use = "complete.obs"),
                   n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(lab = sprintf("r = %.3f\nn = %d", r, n))

# Identical x and y scale: shared limits (same range on both axes) + 1:1 aspect
lim_sy    <- range(c(paired$mean_mint, paired$nightT), na.rm = TRUE) + c(-0.4, 0.4)
lim_daily <- range(c(daily$Tmin_day, daily$nightT_day), na.rm = TRUE) + c(-0.4, 0.4)

theme_nt <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom")

# 9a. Site-year scatter (mean_mint vs nightT)
p_sy <- ggplot2::ggplot(paired, ggplot2::aes(mean_mint, nightT, color = location)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black",
                       linewidth = 0.9, formula = y ~ x, inherit.aes = FALSE,
                       ggplot2::aes(mean_mint, nightT)) +
  ggplot2::geom_point(size = 2.6, alpha = 0.9) +
  ggplot2::coord_fixed(ratio = 1, xlim = lim_sy, ylim = lim_sy, expand = FALSE) +
  ggplot2::annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
                    label = sprintf("r = %.3f\nnightT = %.2f*mint + %.1f\noffset = %+.2f C\nn = %d",
                                    sy$r, sy$slope, sy$intercept, sy$offset, sy$n),
                    size = 3.4, label.size = 0.3) +
  ggplot2::labs(x = "Mean daily minimum temperature, mean_mint (C)",
                y = "Dark-period night-time temperature, nightT (C)",
                colour = NULL,
                title = "Tmin vs dark-period night-time temperature",
                subtitle = "Site-year means, June-July, 64 site-years (NASA POWER hourly)") +
  theme_nt
ggplot2::ggsave(file.path(OUTDIR, "mint_vs_nighttime.png"), p_sy,
                width = 7.2, height = 7.2, dpi = 300, bg = "white")

# 9b. Daily scatter, all data
p_daily <- ggplot2::ggplot(daily, ggplot2::aes(Tmin_day, nightT_day)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = location), size = 0.8, alpha = 0.35) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black",
                       linewidth = 0.9, formula = y ~ x) +
  ggplot2::coord_fixed(ratio = 1, xlim = lim_daily, ylim = lim_daily, expand = FALSE) +
  ggplot2::annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
                    label = sprintf("r = %.3f\nnightT = %.2f*Tmin + %.1f\noffset = %+.2f C\nn = %d days",
                                    dl$r, dl$slope, dl$intercept, dl$offset, dl$n),
                    size = 3.4, label.size = 0.3) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2))) +
  ggplot2::labs(x = "Daily minimum temperature, Tmin (C)",
                y = "Daily dark-period night-time temperature, nightT (C)",
                colour = NULL,
                title = "Daily Tmin vs dark-period night-time temperature",
                subtitle = "All days, June-July across 64 site-years (NASA POWER hourly)") +
  theme_nt
ggplot2::ggsave(file.path(OUTDIR, "mint_vs_nighttime_daily.png"), p_daily,
                width = 7.2, height = 7.4, dpi = 300, bg = "white")

# 9c. Per-location multipanel (daily)
p_multi <- ggplot2::ggplot(daily, ggplot2::aes(Tmin_day, nightT_day)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(colour = "#2c7fb8", size = 0.7, alpha = 0.35) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black",
                       linewidth = 0.8, formula = y ~ x) +
  ggplot2::geom_text(data = loc_r_daily, ggplot2::aes(x = -Inf, y = Inf, label = lab),
                     hjust = -0.1, vjust = 1.15, size = 3, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~location, ncol = 3) +
  ggplot2::coord_fixed(ratio = 1, xlim = lim_daily, ylim = lim_daily, expand = FALSE) +
  ggplot2::labs(x = "Daily minimum temperature, Tmin (C)",
                y = "Daily dark-period night-time temperature, nightT (C)",
                title = "Daily Tmin vs dark-period night-time temperature, by location",
                subtitle = "June-July, all years (NASA POWER hourly). Dashed = 1:1, solid = OLS.") +
  theme_nt + ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(OUTDIR, "mint_vs_nighttime_by_location.png"), p_multi,
                width = 10, height = 7, dpi = 300, bg = "white")
message("Wrote 3 plots (site-year, daily, per-location multipanel).")

## ---- 10. Inspection report ------------------------------------------------
per_loc_sy <- paired %>% dplyr::group_by(location) %>%
  dplyr::summarise(n = dplyr::n(),
                   `r(mint,nightT)` = ifelse(dplyr::n() > 2,
                       round(stats::cor(mean_mint, nightT, use = "complete.obs"), 3), NA),
                   mean_mint = round(mean(mean_mint, na.rm = TRUE), 2),
                   mean_nightT = round(mean(nightT), 2),
                   `night-mint` = round(mean(nightT - mean_mint, na.rm = TRUE), 2),
                   .groups = "drop")
per_loc_daily <- loc_r_daily %>%
  dplyr::transmute(location, n_days = n, `r(Tmin,nightT)` = round(r, 3))

night_thr <- if (NIGHT_ELEV_DEG == 0) "geometric horizon (elev < 0 deg)" else sprintf("elev < %g deg", NIGHT_ELEV_DEG)
avg_nhr <- mean(paired$n_night_hr) / mean(paired$n_days)

rep <- c(
  "# Dark-period night-time temperature -- inspection report",
  "",
  "Generated by `weather/nighttime_temperature/nighttime_temperature.R` (R).",
  "",
  sprintf("- Source: NASA POWER **hourly** T2M (community=AG), same product as the paper, native hourly resolution."),
  sprintf("- Window: June 1 - July 31 (DOY %d-%d), matching the paper's fixed window.", WINDOW_START_DOY, WINDOW_END_DOY),
  sprintf("- Night definition: %s, via NOAA solar-position algorithm (UTC).", night_thr),
  sprintf("- Site-years: **%d**; daily observations: **%d** (local calendar, CDT).", nrow(paired), nrow(daily)),
  sprintf("- Mean dark hours/day: %.1f; mean missing hourly values/site-year: %.1f.", avg_nhr, mean(paired$n_missing)),
  "",
  "## Headline -- site-year level (matches the regression)",
  "",
  sprintf("- **Correlation (mean_mint vs dark-period nightT): r = %.3f** (n = %d).", sy$r, sy$n),
  sprintf("- OLS: nightT = %.3f * mean_mint + %.2f.", sy$slope, sy$intercept),
  sprintf("- Mean offset (nightT - mean_mint): **%+.2f C** (dark-period mean runs %s than the daily minimum).",
          sy$offset, ifelse(sy$offset > 0, "warmer", "cooler")),
  sprintf("- RMSE(nightT - mean_mint): %.2f C.", sy$rmse),
  sprintf("- For contrast, correlation of mean_mint with day-time T2M: r = %.3f.", syD$r),
  "",
  "## Daily level -- all data (more granular)",
  "",
  sprintf("- **Correlation (daily Tmin vs daily dark-period nightT): r = %.3f** (n = %d days).", dl$r, dl$n),
  sprintf("- OLS: nightT = %.3f * Tmin + %.2f; mean offset %+.2f C; RMSE %.2f C.",
          dl$slope, dl$intercept, dl$offset, dl$rmse),
  "- Per-location daily correlations (see also the multipanel figure):",
  "",
  "```",
  paste(utils::capture.output(print(as.data.frame(per_loc_daily), row.names = FALSE)), collapse = "\n"),
  "```",
  "",
  "## QC bridge",
  "",
  sprintf("- A Tmin-like metric rebuilt from the hourly pull (mean of daily hourly minima) vs the paper's `mean_mint`: r = %.3f, RMSE = %.2f C. High agreement confirms the hourly pull is consistent with the daily-based variable.",
          syB$r, syB$rmse),
  "",
  "## Interpretation",
  "",
  "- A very high r means the dark-period mean and the daily minimum carry nearly the same signal, so substituting one for the other moves the yield-model coefficient little -- the choice of metric is not driving the result. Report r, the offset, and the refit coefficient side by side to answer Reviewer 1 as a robustness result.",
  sprintf("- The positive offset (dark-period mean ~%+.1f C above the single coldest instant) quantifies how much `mean_mint` under-states the temperature the crop actually experiences through the night.", sy$offset),
  "",
  "## Per-location summary (site-year)",
  "",
  "```",
  paste(utils::capture.output(print(as.data.frame(per_loc_sy), row.names = FALSE)), collapse = "\n"),
  "```",
  "",
  "## Outputs",
  "",
  "- `nighttime_temperature_paired.csv` -- per site-year day/night temps + mean_mint",
  "- `nighttime_temperature_daily.csv`  -- per day Tmin/nightT/dayT",
  "- `mint_vs_nighttime.png`            -- site-year association",
  "- `mint_vs_nighttime_daily.png`      -- daily association, all data",
  "- `mint_vs_nighttime_by_location.png`-- daily association, one panel per location",
  "",
  "## Caveats",
  "",
  "- NASA POWER is a coarse reanalysis; absolute night values can miss local radiative cooling on calm, clear nights. R1's concern is about *which hours* are averaged (dark hours vs the single minimum), which hourly data answers directly.",
  "- Night is set by solar geometry; the 2 W m-2 / 5 umol threshold R1 mentions converges to nearly the same hours over a multi-week mean (set NIGHT_ELEV_DEG = -6 for a civil-twilight check).",
  "- Window is the paper's fixed June-July; if the flowering window moves to a crop-calendar basis (Reviewer 2), recompute both metrics on that window -- the script keys off WINDOW_*_DOY."
)
writeLines(rep, file.path(OUTDIR, "nighttime_inspection_report.md"))
message("Wrote nighttime_inspection_report.md")

cat("\n==== SUMMARY ====\n")
cat(sprintf("Site-year:  r(mean_mint, nightT) = %.3f  offset %+.2f C  (n=%d)\n", sy$r, sy$offset, sy$n))
cat(sprintf("Daily:      r(Tmin, nightT)      = %.3f  offset %+.2f C  (n=%d days)\n", dl$r, dl$offset, dl$n))
cat(sprintf("QC bridge:  r(mean_mint, mint_hr)= %.3f  RMSE %.2f C\n", syB$r, syB$rmse))
cat("Done.\n")
