# Dark-period night-time temperature (R project)

Computes the **actual dark-period night-time temperature** the crop experiences —
the mean of hourly 2 m air temperature while the sun is below the horizon — from
**NASA POWER hourly** data, and compares it against the mean daily minimum
temperature (`mean_mint`) used in the manuscript. This directly addresses
Reviewer 1's request to use dark-hour values rather than the single daily minimum.

## What it produces

Written to `../outputs/` (i.e. `weather/outputs/`):

| File | Contents |
|------|----------|
| `nighttime_temperature_paired.csv` | Per **site-year**: `nightT`, `dayT`, `mean_mint`, offsets, hour counts |
| `nighttime_temperature_daily.csv`  | Per **day**: `Tmin_day`, `nightT_day`, `dayT_day` |
| `mint_vs_nighttime.png`            | Site-year association (`mean_mint` vs `nightT`) |
| `mint_vs_nighttime_daily.png`      | Daily association, all days pooled |
| `mint_vs_nighttime_by_location.png`| Daily association, one panel per location |
| `nighttime_inspection_report.md`   | Correlations, offset, QC bridge, per-location tables, caveats |

## How to run (anywhere)

**RStudio** — open `nighttime_temperature.Rproj`, open `nighttime_temperature.R`,
click **Source**.

**Terminal**

```bash
Rscript nighttime_temperature.R
```

The script locates the repository root on its own (it looks upward for the
`corn variety testing` folder), so it works regardless of the working directory
as long as it lives inside the repo.

## Inputs (read from the repo)

- `weather/locations.june.july.xlsx` — site coordinates and the June–July window.
- `corn variety testing/mixed_effect_models/input/Corn data organized2.csv` —
  the modelling dataset, source of the 64 site-years and `mean_mint`.

## Dependencies

`jsonlite`, `curl`, `readxl`, `dplyr`, `tidyr`, `ggplot2`, `scales`. Missing
packages are installed automatically on first run. No `nasapower`/`suncalc`
needed — the POWER REST API call and the NOAA solar-position math are inline.

## How it works

1. For each site-year, pull hourly `T2M` for 1 Jun–31 Jul (cached under
   `../outputs/cache/`, so re-runs are instant and offline-friendly).
2. Flag each hour day/night by solar elevation (`< 0°` = below the horizon),
   using the standard NOAA solar-position algorithm.
3. Average dark hours → `nightT`; also compute a daily table for the granular
   view. Pair `nightT` with the existing `mean_mint` and plot/report.

## Configuration

All knobs are at the top of the script: `WINDOW_START_DOY`/`WINDOW_END_DOY`
(exposure window), `NIGHT_ELEV_DEG` (`0` = horizon, `-6` = civil twilight for a
sensitivity check), `LOCAL_UTC_OFFSET`, and the location-name map.
