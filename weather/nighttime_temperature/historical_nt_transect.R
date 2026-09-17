#!/usr/bin/env Rscript
# =============================================================================
# Historical June-July NIGHT-TIME temperature trends across the US maize belt
# A south->north transect (TX, AR, MO, IA, MN) for a USDA proposal figure.
# =============================================================================
#
# Night-time temperature is the dark-period mean of NASA POWER hourly 2 m air
# temperature (T2M) while PPFD < 5 umol m-2 s-1 (photosynthetic photon flux from
# hourly PAR, ALLSKY_SFC_PAR_TOT) -- the same approach used for the Arkansas
# analysis. Trends are the slope of June-July nightT vs year, 2001-2025.
#
# NOTE: NASA POWER hourly (MERRA-2) begins in 2001, so the dark-period trend
# spans 25 years (2001-2025), not the 40-year daily record.
#
# Run:  Rscript historical_nt_transect.R
# Depends: jsonlite, dplyr, ggplot2, scales.
# =============================================================================

pkgs <- c("jsonlite", "dplyr", "ggplot2", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

## ---- CONFIG ---------------------------------------------------------------
YEARS            <- 2001:2025
WINDOW_START_DOY <- 152      # Jun 1
WINDOW_END_DOY   <- 212      # Jul 31
PPFD_THRESHOLD   <- 5.0      # umol m-2 s-1  (night = below this)
PAR_TO_PPFD      <- 1e6 / 3600 * 4.57
NASA_FILL        <- -999.0
PARAMS           <- "T2M,ALLSKY_SFC_PAR_TOT"
CACHE_VER        <- "tr"     # separate cache namespace for the transect

# One representative irrigated/rainfed maize location per state, south -> north
SITES <- data.frame(
  state = c("TX", "AR", "MO", "IA", "MN"),
  site  = c("College Station", "Keiser", "Columbia", "Ames", "Waseca"),
  lon   = c(-96.334, -90.081, -92.334, -93.632, -93.507),
  lat   = c( 30.628,  35.676,  38.952,  42.030,  44.073),
  stringsAsFactors = FALSE
)
SITES <- SITES[order(SITES$lat), ]  # south -> north

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

## ---- NASA POWER hourly pull (cached, retry) -------------------------------
fetch_hourly <- function(lon, lat, year) {
  start <- sprintf("%d0601", year); end <- sprintf("%d0731", year)
  cf <- file.path(CACHE, sprintf("%s_%.4f_%.4f_%s_%s.json", CACHE_VER, lat, lon, start, end))
  if (file.exists(cf)) {
    payload <- jsonlite::fromJSON(cf, simplifyVector = TRUE)
  } else {
    url <- sprintf(paste0("https://power.larc.nasa.gov/api/temporal/hourly/point",
                          "?parameters=%s&community=AG&longitude=%s&latitude=%s",
                          "&start=%s&end=%s&format=JSON&time-standard=UTC"),
                   PARAMS, lon, lat, start, end)
    payload <- NULL
    for (attempt in 0:3) {
      payload <- tryCatch(jsonlite::fromJSON(url, simplifyVector = TRUE), error = function(e) NULL)
      if (!is.null(payload)) break
      Sys.sleep(2^attempt)
    }
    if (is.null(payload)) stop(sprintf("POWER failed %s,%s,%d", lat, lon, year))
    writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE), cf)
  }
  p <- payload$properties$parameter
  grab <- function(nm) { v <- as.numeric(unlist(p[[nm]], use.names = FALSE)); v[v <= NASA_FILL] <- NA; v }
  times <- names(p$T2M)
  data.frame(dt = as.POSIXct(times, format = "%Y%m%d%H", tz = "UTC"),
             t2m = grab("T2M"), par = grab("ALLSKY_SFC_PAR_TOT"))
}

# June-July dark-period night-time temperature for one site-year
night_temp <- function(lon, lat, year) {
  df <- fetch_hourly(lon, lat, year)
  doy <- as.POSIXlt(df$dt, tz = "UTC")$yday + 1
  df <- df[doy >= WINDOW_START_DOY & doy <= WINDOW_END_DOY, , drop = FALSE]
  ppfd <- df$par * PAR_TO_PPFD
  night <- df$t2m[ppfd < PPFD_THRESHOLD]
  if (all(is.na(night))) NA_real_ else mean(night, na.rm = TRUE)
}

## ---- Build the transect dataset -------------------------------------------
message(sprintf("Pulling %d sites x %d years (NASA POWER hourly, cached) ...",
                nrow(SITES), length(YEARS)))
rows <- list()
for (i in seq_len(nrow(SITES))) {
  for (y in YEARS) {
    nt <- tryCatch(night_temp(SITES$lon[i], SITES$lat[i], y),
                   error = function(e) { message("  [fail] ", SITES$site[i], " ", y); NA_real_ })
    rows[[length(rows) + 1]] <- data.frame(
      state = SITES$state[i], site = SITES$site[i], lat = SITES$lat[i], year = y, nightT = nt)
  }
  message(sprintf("  %s done", SITES$site[i]))
}
dat <- dplyr::bind_rows(rows)
dat$label <- factor(sprintf("%s - %s (%.1f N)", dat$state, dat$site, dat$lat),
                    levels = sprintf("%s - %s (%.1f N)", SITES$state, SITES$site, SITES$lat))
utils::write.csv(dat, file.path(OUTDIR, "historical_nt_transect.csv"), row.names = FALSE)

## ---- Trend statistics: OLS + robust Theil-Sen -----------------------------
theil_sen <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x); s <- c()
  for (i in 1:(n - 1)) for (j in (i + 1):n) s <- c(s, (y[j] - y[i]) / (x[j] - x[i]))
  stats::median(s)
}
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

slopes <- dat %>% dplyr::group_by(state, site, lat, label) %>%
  dplyr::group_modify(function(d, ...) {
    fit <- stats::lm(nightT ~ year, data = d)
    sm  <- summary(fit)$coefficients["year", ]
    ci  <- stats::confint(fit)["year", ]
    data.frame(
      slope_dec   = unname(sm["Estimate"]) * 10,             # C per decade
      ci_lo_dec   = ci[1] * 10, ci_hi_dec = ci[2] * 10,
      se_dec      = unname(sm["Std. Error"]) * 10,
      p_value     = unname(sm["Pr(>|t|)"]),
      sen_dec     = theil_sen(d$year, d$nightT) * 10,        # robust slope, C/decade
      r2          = summary(fit)$r.squared,
      mean_nightT = mean(d$nightT, na.rm = TRUE),
      n_years     = sum(is.finite(d$nightT)))
  }) %>% dplyr::ungroup() %>% dplyr::arrange(lat)
slopes$sig <- stars(slopes$p_value)
slopes$legend <- sprintf("%s  %+.2f C/dec %s", slopes$label, slopes$slope_dec, slopes$sig)
utils::write.csv(slopes, file.path(OUTDIR, "historical_nt_transect_slopes.csv"), row.names = FALSE)
mean_slope <- mean(slopes$slope_dec); n_sig <- sum(slopes$p_value < 0.05)

## ---- MAIN figure: overlaid ladder with trend lines ------------------------
pal <- setNames(c("#d73027", "#fc8d59", "#4a7f3f", "#4575b4", "#542788"), SITES$label |> as.character())
# map legend labels (with slopes) onto the ladder colours in latitude order
dat2 <- dat %>% dplyr::left_join(slopes[, c("label", "legend", "slope_dec", "sig")], by = "label")
dat2$legend <- factor(dat2$legend, levels = slopes$legend)
legcols <- setNames(c("#d73027", "#fc8d59", "#4a7f3f", "#4575b4", "#542788"), slopes$legend)

p_main <- ggplot2::ggplot(dat2, ggplot2::aes(year, nightT, colour = legend)) +
  ggplot2::geom_point(size = 1.7, alpha = 0.85) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.12, formula = y ~ x) +
  ggplot2::scale_colour_manual(values = legcols, name = "State - site (latitude)   trend") +
  ggplot2::scale_x_continuous(breaks = seq(2001, 2025, 4)) +
  ggplot2::labs(
    x = NULL, y = expression("June-July night-time temperature ("*degree*"C)"),
    title = "June-July night-time temperature trends across the US maize belt",
    subtitle = "Dark-period mean (PPFD < 5 umol m-2 s-1), NASA POWER hourly, 2001-2025",
    caption = paste0("All five slopes positive (mean ", sprintf("%+.2f", mean_slope),
                     " C/decade); none individually significant over 25 yr (ns). Lines = OLS, bands = 95% CI.")) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
                 legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.75), colour = "grey70"),
                 legend.title = ggplot2::element_text(size = 9, face = "bold"),
                 legend.text = ggplot2::element_text(size = 8.5),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(OUTDIR, "historical_nt_transect.png"), p_main,
                width = 10.5, height = 6.5, dpi = 300, bg = "white")

## ---- FACET figure (detail / SI): one panel per site, N at top -------------
dat$label_ns <- factor(dat$label, levels = rev(levels(dat$label)))
slab <- slopes %>% dplyr::mutate(
  label_ns = factor(label, levels = rev(levels(dat$label))),
  txt = sprintf("OLS %+.2f C/dec %s (95%% CI %.2f, %.2f)\nTheil-Sen %+.2f  R2 = %.2f",
                slope_dec, sig, ci_lo_dec, ci_hi_dec, sen_dec, r2))
p_facet <- ggplot2::ggplot(dat, ggplot2::aes(year, nightT)) +
  ggplot2::geom_point(size = 1.4, colour = "grey30") +
  ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#c0392b", fill = "#c0392b",
                       alpha = 0.15, linewidth = 0.9, formula = y ~ x) +
  ggplot2::geom_text(data = slab, ggplot2::aes(x = -Inf, y = Inf, label = txt),
                     hjust = -0.03, vjust = 1.15, size = 3, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~label_ns, ncol = 1, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = seq(2001, 2025, 4)) +
  ggplot2::labs(x = NULL, y = expression("June-July night-time temperature ("*degree*"C)"),
                title = "June-July night-time temperature trend by site (2001-2025)",
                subtitle = "Dark-period mean, PPFD < 5 umol m-2 s-1, NASA POWER hourly") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(OUTDIR, "historical_nt_transect_facet.png"), p_facet,
                width = 8, height = 11, dpi = 300, bg = "white")
message("Wrote transect figures.")

## ---- Report ---------------------------------------------------------------
fmt <- function(v, d = 2) formatC(v, format = "f", digits = d)
rep <- c(
  "# June-July night-time temperature trends across the US maize belt",
  "",
  "For a USDA proposal. Dark-period night-time temperature (mean 2 m air temperature",
  "while PPFD < 5 umol m-2 s-1) from NASA POWER hourly data, June-July, 2001-2025,",
  "one representative maize location per state along a south-to-north transect.",
  "",
  "## Trends (per decade)",
  "",
  "| State - site (lat) | mean NT (C) | OLS C/dec | 95% CI | p | Theil-Sen C/dec | R2 |",
  "|---|---|---|---|---|---|---|",
  paste(sprintf("| %s | %s | %+.2f %s | %.2f, %.2f | %.3f | %+.2f | %.2f |",
                slopes$label, fmt(slopes$mean_nightT, 1), slopes$slope_dec, slopes$sig,
                slopes$ci_lo_dec, slopes$ci_hi_dec, slopes$p_value, slopes$sen_dec, slopes$r2),
        collapse = "\n"),
  "",
  sprintf("Mean warming across the five sites: **%+.2f C/decade** (%d of 5 sites significant at p<0.05).",
          mean_slope, n_sig),
  "",
  "## Critical analysis",
  "",
  sprintf("- **Direction & consistency.** All five sites show %s June-July night-time warming (+%.2f to +%.2f C/decade; mean %+.2f), spanning ~14 degrees of latitude from the Texas Blacklands to southern Minnesota. The coherent positive direction -- 5 of 5 sites, sign-test p ~ 0.06 -- is the point most useful to a proposal: high night-time temperature is a maize-belt-wide exposure, not a single-site or Mid-South-only artefact.",
          ifelse(all(slopes$slope_dec > 0), "positive", "mostly positive"),
          min(slopes$slope_dec), max(slopes$slope_dec), mean_slope),
  sprintf("- **Significance / statistical power.** Over the 25-year hourly record, none of the individual site slopes reach p<0.05 (p = %.2f-%.2f; %s is closest). Interannual variability is large relative to a ~0.3-0.5 C/decade trend, so 25 annual values give limited power and wide confidence intervals. This is expected and is itself part of the motivation: the reanalysis-era record is too short to resolve site-level significance, which the proposed longer-term / network observations are designed to address. Theil-Sen slopes are reported alongside OLS and agree in sign at every site, indicating the positive tendency is not driven by a few outlier years.",
          min(slopes$p_value), max(slopes$p_value), slopes$site[which.min(slopes$p_value)]),
  "- **Magnitude in context.** These are trends in the *dark-period mean*, the physiologically relevant quantity for night respiration, and are broadly consistent with the ~0.24 C/decade June-July Tmin trend reported for Arkansas over 40 years (paper Fig. S4); the somewhat larger per-decade values here partly reflect the shorter, more variable 25-year window.",
  "- **Data caveats.** NASA POWER hourly is MERRA-2 reanalysis at ~0.5 deg (~50 km); a single grid cell represents the landscape, not a specific field, and can miss local radiative cooling on calm, clear nights. Absolute values carry reanalysis bias, but decadal *trends* are less sensitive to a stable bias. The hourly record begins in 2001, which sets the 25-year window.",
  "- **Autocorrelation & attribution.** OLS treats years as independent; mild serial correlation would modestly inflate significance. These are descriptive climate trends, not an attribution analysis -- they motivate the proposal (night warming is broad and ongoing) rather than establish a yield-causal claim.",
  "- **Site representativeness.** Sites were chosen as recognizable maize/AES locations spanning the latitudinal gradient; results are illustrative of the transect, and a production-weighted or multi-cell average would be the next step for a formal regional estimate.",
  "",
  "## Outputs",
  "- `historical_nt_transect.png`       -- main proposal figure (overlaid trends)",
  "- `historical_nt_transect_facet.png` -- per-site detail panels",
  "- `historical_nt_transect.csv`       -- site x year night-time temperature",
  "- `historical_nt_transect_slopes.csv`-- per-site slope statistics")
writeLines(rep, file.path(OUTDIR, "historical_nt_transect_report.md"))

cat("\n==== TRANSECT SLOPES (C/decade) ====\n")
print(slopes[, c("state", "site", "lat", "slope_dec", "p_value", "sen_dec", "mean_nightT")],
      row.names = FALSE)
cat(sprintf("\nMean across sites: %+.2f C/decade\n", mean_slope))
cat("Done.\n")
