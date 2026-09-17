# June-July MINIMUM temperature trends across the US maize belt (long record)

Daily minimum temperature (NASA POWER DAILY T2M_MIN), June-July, 1981-2025 (~45 years),
one maize site per state, south-to-north. Companion to the 25-year hourly dark-period version.

## Trends (absolute C/decade and relative %/decade)

| State - site (lat) | mean Tmin (C) | abs C/dec | 95% CI | rel %/dec | p | Theil-Sen C/dec | R2 |
|---|---|---|---|---|---|---|---|
| TX - College Station (30.6 N) | 23.1 | +0.33 *** | 0.17, 0.49 | +1.4% | 0.0002 | +0.31 | 0.28 |
| AR - Keiser (35.7 N) | 21.1 | +0.30 ** | 0.10, 0.49 | +1.4% | 0.0039 | +0.33 | 0.18 |
| MO - Columbia (39.0 N) | 19.0 | +0.38 *** | 0.17, 0.60 | +2.0% | 0.0008 | +0.37 | 0.23 |
| IA - Ames (42.0 N) | 17.0 | +0.40 ** | 0.16, 0.65 | +2.4% | 0.0018 | +0.43 | 0.21 |
| MN - Waseca (44.1 N) | 15.6 | +0.37 ** | 0.12, 0.63 | +2.4% | 0.0050 | +0.42 | 0.17 |

Mean across sites: **+0.36 C/decade (+1.9%/decade)**; 5 of 5 sites significant at p<0.05.

## Relative slope -- why it is added

- The absolute slope (C/decade) is the physical warming rate. The **relative slope = absolute / mean(Tmin) x 100** expresses that warming as a percent of each site's baseline night temperature per decade, which makes sites with very different baselines comparable (e.g. a +0.3 C/decade rise is a larger *fractional* change on a cool 16 C Minnesota night than on a warm 24 C Texas night).
- Reported per decade here; divide by 10 for %/year.

## Critical analysis

- **Longer record, more power.** Over ~45 years, 5 of 5 sites reach p<0.05 -- a marked gain over the 25-year hourly dark-period version, where none did. The extra two decades pull the trend out of the interannual noise. This is the key reason to show the Tmin long-record version alongside the physically-preferred hourly one.
- **What Tmin is (and isn't).** Tmin is the single coldest instant of each day (near sunrise), not the dark-period mean the crop experiences; it is the paper's variable and the only one with a multi-decade record. The hourly analysis showed dark-period NT is ~1:1 collinear with Tmin (r~0.98), so the Tmin trend is a good proxy for the night-temperature trend while extending the record back to 1981.
- **Absolute vs relative gradient.** Compare the two columns: because northern sites have lower baseline Tmin, a similar absolute slope translates into a larger relative (%) trend there; note whether the warming is spatially uniform in absolute or in relative terms.
- **Robustness.** Theil-Sen slopes accompany OLS; agreement in sign and rough magnitude indicates trends are not driven by a few anomalous years (e.g. the 2011-2012 heat/drought).
- **Data caveats.** NASA POWER daily is MERRA-2-based reanalysis at ~0.5 deg; a grid cell represents the landscape, not a specific field, and absolute values carry reanalysis bias (trends are less sensitive to a stable bias). OLS treats years as independent; mild autocorrelation would modestly inflate significance. Single site per state is illustrative of the latitudinal transect, not a production-weighted state value.

## Outputs
- `historical_tmin_transect.png`        -- main figure (absolute + relative trends in legend)
- `historical_tmin_transect_facet.png`  -- per-site detail panels
- `historical_tmin_transect.csv`        -- site x year June-July mean Tmin
- `historical_tmin_transect_slopes.csv` -- per-site absolute & relative slope statistics
