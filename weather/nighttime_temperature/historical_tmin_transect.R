#!/usr/bin/env Rscript
# =============================================================================
# Historical June-July MINIMUM temperature trends across the US maize belt
# Longer-record companion to historical_nt_transect.R.
# =============================================================================
#
# This version uses daily minimum temperature (NASA POWER DAILY T2M_MIN), which
# reaches back to 1981 -- ~45 years vs the 25-year hourly record used for the
# dark-period night-time analysis. Tmin is the paper's actual variable ("mean
# daily minimum temperature"), so this is the long-record trend behind it.
#
# For each site (south->north transect) we compute June-July mean Tmin per year,
# then the trend as:
#   * ABSOLUTE slope  (C per decade)          -- OLS and robust Theil-Sen
#   * RELATIVE slope  (% per decade)          -- absolute / mean(Tmin) * 100
# The relative slope makes warming comparable across sites with very different
# baseline night temperatures (TX ~24 C vs MN ~16 C).
#
# Run:  Rscript historical_tmin_transect.R
# Depends: jsonlite, dplyr, ggplot2, scales.
# =============================================================================

pkgs <- c("jsonlite", "dplyr", "ggplot2", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

## ---- CONFIG ---------------------------------------------------------------
YEAR_START       <- 1981     # NASA POWER daily record begins 1981-01-01
YEAR_END         <- 2025
WINDOW_START_DOY <- 152      # Jun 1
WINDOW_END_DOY   <- 212      # Jul 31
NASA_FILL        <- -999.0
PARAM            <- "T2M_MIN"
CACHE_VER        <- "tmin"

SITES <- data.frame(
  state = c("TX", "AR", "MO", "IA", "MN"),
  site  = c("College Station", "Keiser", "Columbia", "Ames", "Waseca"),
  lon   = c(-96.334, -90.081, -92.334, -93.632, -93.507),
  lat   = c( 30.628,  35.676,  38.952,  42.030,  44.073),
  stringsAsFactors = FALSE)
SITES <- SITES[order(SITES$lat), ]

## ---- paths ----------------------------------------------------------------
get_script_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()
}
find_repo_root <- function(start) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:12) { if (dir.exists(file.path(d, "corn variety testing"))) return(d)
    p <- dirname(d); if (identical(p, d)) break; d <- p }
  stop("Could not locate repo root.")
}
ROOT   <- find_repo_root(get_script_dir())
OUTDIR <- file.path(ROOT, "weather", "outputs")
CACHE  <- file.path(OUTDIR, "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

## ---- NASA POWER daily pull (one call per site, full range, cached) --------
fetch_daily_tmin <- function(lon, lat) {
  start <- sprintf("%d0101", YEAR_START); end <- sprintf("%d1231", YEAR_END)
  cf <- file.path(CACHE, sprintf("%s_%.4f_%.4f_%s_%s.json", CACHE_VER, lat, lon, start, end))
  if (file.exists(cf)) {
    payload <- jsonlite::fromJSON(cf, simplifyVector = TRUE)
  } else {
    url <- sprintf(paste0("https://power.larc.nasa.gov/api/temporal/daily/point",
                          "?parameters=%s&community=AG&longitude=%s&latitude=%s",
                          "&start=%s&end=%s&format=JSON"), PARAM, lon, lat, start, end)
    payload <- NULL
    for (attempt in 0:3) {
      payload <- tryCatch(jsonlite::fromJSON(url, simplifyVector = TRUE), error = function(e) NULL)
      if (!is.null(payload)) break
      Sys.sleep(2^attempt)
    }
    if (is.null(payload)) stop(sprintf("POWER daily failed %s,%s", lat, lon))
    writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE), cf)
  }
  v <- payload$properties$parameter$T2M_MIN         # named: "YYYYMMDD" -> Tmin
  d <- data.frame(date = names(v), tmin = as.numeric(unlist(v, use.names = FALSE)))
  d$tmin[d$tmin <= NASA_FILL] <- NA
  d$year <- as.integer(substr(d$date, 1, 4))
  d$doy  <- as.integer(format(as.Date(d$date, "%Y%m%d"), "%j"))
  d
}

# June-July mean Tmin per year for a site
jj_tmin_by_year <- function(lon, lat) {
  d <- fetch_daily_tmin(lon, lat)
  d <- d[d$doy >= WINDOW_START_DOY & d$doy <= WINDOW_END_DOY, ]
  d %>% dplyr::group_by(year) %>%
    dplyr::summarise(tmin = mean(tmin, na.rm = TRUE), .groups = "drop")
}

## ---- Build dataset --------------------------------------------------------
message(sprintf("Pulling %d sites, daily Tmin %d-%d (one call/site, cached) ...",
                nrow(SITES), YEAR_START, YEAR_END))
rows <- list()
for (i in seq_len(nrow(SITES))) {
  yy <- tryCatch(jj_tmin_by_year(SITES$lon[i], SITES$lat[i]),
                 error = function(e) { message("  [fail] ", SITES$site[i], ": ", conditionMessage(e)); NULL })
  if (is.null(yy)) next
  yy$state <- SITES$state[i]; yy$site <- SITES$site[i]; yy$lat <- SITES$lat[i]
  rows[[length(rows) + 1]] <- yy
  message(sprintf("  %s done (%d years)", SITES$site[i], nrow(yy)))
}
dat <- dplyr::bind_rows(rows)
dat$label <- factor(sprintf("%s - %s (%.1f N)", dat$state, dat$site, dat$lat),
                    levels = sprintf("%s - %s (%.1f N)", SITES$state, SITES$site, SITES$lat))
utils::write.csv(dat, file.path(OUTDIR, "historical_tmin_transect.csv"), row.names = FALSE)

## ---- Trends: OLS + Theil-Sen; ABSOLUTE (C/dec) and RELATIVE (%/dec) --------
theil_sen <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x); s <- c()
  for (i in 1:(n - 1)) for (j in (i + 1):n) s <- c(s, (y[j] - y[i]) / (x[j] - x[i]))
  stats::median(s)
}
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

slopes <- dat %>% dplyr::group_by(state, site, lat, label) %>%
  dplyr::group_modify(function(d, ...) {
    fit <- stats::lm(tmin ~ year, data = d); sm <- summary(fit)$coefficients["year", ]
    ci <- stats::confint(fit)["year", ]; mean_tmin <- mean(d$tmin, na.rm = TRUE)
    abs_dec <- unname(sm["Estimate"]) * 10
    data.frame(
      mean_tmin   = mean_tmin,
      abs_dec     = abs_dec,                                   # C / decade
      abs_ci_lo   = ci[1] * 10, abs_ci_hi = ci[2] * 10,
      p_value     = unname(sm["Pr(>|t|)"]),
      sen_dec     = theil_sen(d$year, d$tmin) * 10,            # robust C / decade
      rel_dec     = abs_dec / mean_tmin * 100,                 # % / decade
      rel_ci_lo   = ci[1] * 10 / mean_tmin * 100,
      rel_ci_hi   = ci[2] * 10 / mean_tmin * 100,
      r2          = summary(fit)$r.squared,
      n_years     = sum(is.finite(d$tmin)))
  }) %>% dplyr::ungroup() %>% dplyr::arrange(lat)
slopes$sig <- stars(slopes$p_value)
slopes$legend <- sprintf("%s  %+.2f C/dec (%+.1f%%/dec) %s",
                         slopes$label, slopes$abs_dec, slopes$rel_dec, slopes$sig)
utils::write.csv(slopes, file.path(OUTDIR, "historical_tmin_transect_slopes.csv"), row.names = FALSE)
mean_abs <- mean(slopes$abs_dec); mean_rel <- mean(slopes$rel_dec); n_sig <- sum(slopes$p_value < 0.05)

## ---- MAIN figure: ladder with trend lines ---------------------------------
dat2 <- dat %>% dplyr::left_join(slopes[, c("label", "legend")], by = "label")
dat2$legend <- factor(dat2$legend, levels = slopes$legend)
legcols <- setNames(c("#d73027", "#fc8d59", "#4a7f3f", "#4575b4", "#542788"), slopes$legend)

p_main <- ggplot2::ggplot(dat2, ggplot2::aes(year, tmin, colour = legend)) +
  ggplot2::geom_point(size = 1.5, alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.12, formula = y ~ x) +
  ggplot2::scale_colour_manual(values = legcols, name = "State - site (lat)   absolute (relative) trend") +
  ggplot2::scale_x_continuous(breaks = seq(1985, 2025, 10)) +
  ggplot2::labs(
    x = NULL, y = expression("June-July minimum temperature ("*degree*"C)"),
    title = "June-July minimum temperature is rising across the US maize belt",
    subtitle = sprintf("Daily Tmin, NASA POWER, %d-%d (45-yr record). Mean trend %+.2f C/decade (%+.1f%%/decade).",
                       YEAR_START, YEAR_END, mean_abs, mean_rel),
    caption = "Trend = OLS slope; significance *** p<0.001, ** p<0.01, * p<0.05, ns. Relative = absolute/mean(Tmin)x100. Bands = 95% CI.") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = c(0.015, 0.985), legend.justification = c(0, 1),
                 legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.78), colour = "grey70"),
                 legend.title = ggplot2::element_text(size = 8.5, face = "bold"),
                 legend.text = ggplot2::element_text(size = 8),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(OUTDIR, "historical_tmin_transect.png"), p_main,
                width = 10.5, height = 6.5, dpi = 300, bg = "white")

## ---- FACET figure ---------------------------------------------------------
dat$label_ns <- factor(dat$label, levels = rev(levels(dat$label)))
slab <- slopes %>% dplyr::mutate(
  label_ns = factor(label, levels = rev(levels(dat$label))),
  txt = sprintf("%+.2f C/dec %s (95%% CI %.2f, %.2f)\nrelative %+.1f%%/dec  |  Theil-Sen %+.2f  R2=%.2f",
                abs_dec, sig, abs_ci_lo, abs_ci_hi, rel_dec, sen_dec, r2))
p_facet <- ggplot2::ggplot(dat, ggplot2::aes(year, tmin)) +
  ggplot2::geom_point(size = 1.2, colour = "grey30") +
  ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#c0392b", fill = "#c0392b",
                       alpha = 0.15, linewidth = 0.9, formula = y ~ x) +
  ggplot2::geom_text(data = slab, ggplot2::aes(x = -Inf, y = Inf, label = txt),
                     hjust = -0.03, vjust = 1.15, size = 3, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~label_ns, ncol = 1, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = seq(1985, 2025, 10)) +
  ggplot2::labs(x = NULL, y = expression("June-July minimum temperature ("*degree*"C)"),
                title = sprintf("June-July minimum temperature trend by site (%d-%d)", YEAR_START, YEAR_END),
                subtitle = "Daily Tmin, NASA POWER. Absolute (C/decade) and relative (%/decade) slopes shown.") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(OUTDIR, "historical_tmin_transect_facet.png"), p_facet,
                width = 8, height = 11, dpi = 300, bg = "white")
message("Wrote Tmin transect figures.")

## ---- Report ---------------------------------------------------------------
fmt <- function(v, d = 2) formatC(v, format = "f", digits = d)
rep <- c(
  "# June-July MINIMUM temperature trends across the US maize belt (long record)",
  "",
  sprintf("Daily minimum temperature (NASA POWER DAILY T2M_MIN), June-July, %d-%d (~%d years),",
          YEAR_START, YEAR_END, YEAR_END - YEAR_START + 1),
  "one maize site per state, south-to-north. Companion to the 25-year hourly dark-period version.",
  "",
  "## Trends (absolute C/decade and relative %/decade)",
  "",
  "| State - site (lat) | mean Tmin (C) | abs C/dec | 95% CI | rel %/dec | p | Theil-Sen C/dec | R2 |",
  "|---|---|---|---|---|---|---|---|",
  paste(sprintf("| %s | %s | %+.2f %s | %.2f, %.2f | %+.1f%% | %.4f | %+.2f | %.2f |",
                slopes$label, fmt(slopes$mean_tmin, 1), slopes$abs_dec, slopes$sig,
                slopes$abs_ci_lo, slopes$abs_ci_hi, slopes$rel_dec, slopes$p_value,
                slopes$sen_dec, slopes$r2), collapse = "\n"),
  "",
  sprintf("Mean across sites: **%+.2f C/decade (%+.1f%%/decade)**; %d of 5 sites significant at p<0.05.",
          mean_abs, mean_rel, n_sig),
  "",
  "## Relative slope -- why it is added",
  "",
  "- The absolute slope (C/decade) is the physical warming rate. The **relative slope = absolute / mean(Tmin) x 100** expresses that warming as a percent of each site's baseline night temperature per decade, which makes sites with very different baselines comparable (e.g. a +0.3 C/decade rise is a larger *fractional* change on a cool 16 C Minnesota night than on a warm 24 C Texas night).",
  "- Reported per decade here; divide by 10 for %/year.",
  "",
  "## Critical analysis",
  "",
  sprintf("- **Longer record, more power.** Over ~%d years, %d of 5 sites reach p<0.05 -- a marked gain over the 25-year hourly dark-period version, where none did. The extra two decades pull the trend out of the interannual noise. This is the key reason to show the Tmin long-record version alongside the physically-preferred hourly one.",
          YEAR_END - YEAR_START + 1, n_sig),
  "- **What Tmin is (and isn't).** Tmin is the single coldest instant of each day (near sunrise), not the dark-period mean the crop experiences; it is the paper's variable and the only one with a multi-decade record. The hourly analysis showed dark-period NT is ~1:1 collinear with Tmin (r~0.98), so the Tmin trend is a good proxy for the night-temperature trend while extending the record back to 1981.",
  "- **Absolute vs relative gradient.** Compare the two columns: because northern sites have lower baseline Tmin, a similar absolute slope translates into a larger relative (%) trend there; note whether the warming is spatially uniform in absolute or in relative terms.",
  "- **Robustness.** Theil-Sen slopes accompany OLS; agreement in sign and rough magnitude indicates trends are not driven by a few anomalous years (e.g. the 2011-2012 heat/drought).",
  "- **Data caveats.** NASA POWER daily is MERRA-2-based reanalysis at ~0.5 deg; a grid cell represents the landscape, not a specific field, and absolute values carry reanalysis bias (trends are less sensitive to a stable bias). OLS treats years as independent; mild autocorrelation would modestly inflate significance. Single site per state is illustrative of the latitudinal transect, not a production-weighted state value.",
  "",
  "## Outputs",
  "- `historical_tmin_transect.png`        -- main figure (absolute + relative trends in legend)",
  "- `historical_tmin_transect_facet.png`  -- per-site detail panels",
  "- `historical_tmin_transect.csv`        -- site x year June-July mean Tmin",
  "- `historical_tmin_transect_slopes.csv` -- per-site absolute & relative slope statistics")
writeLines(rep, file.path(OUTDIR, "historical_tmin_transect_report.md"))

cat("\n==== Tmin TRANSECT SLOPES ====\n")
print(slopes[, c("state", "site", "mean_tmin", "abs_dec", "rel_dec", "p_value", "sen_dec")], row.names = FALSE)
cat(sprintf("\nMean: %+.2f C/decade  (%+.1f%%/decade);  %d/5 significant\n", mean_abs, mean_rel, n_sig))
cat("Done.\n")
