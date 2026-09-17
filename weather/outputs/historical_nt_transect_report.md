# June-July night-time temperature trends across the US maize belt

For a USDA proposal. Dark-period night-time temperature (mean 2 m air temperature
while PPFD < 5 umol m-2 s-1) from NASA POWER hourly data, June-July, 2001-2025,
one representative maize location per state along a south-to-north transect.

## Trends (per decade)

| State - site (lat) | mean NT (C) | OLS C/dec | 95% CI | p | Theil-Sen C/dec | R2 |
|---|---|---|---|---|---|---|
| TX - College Station (30.6 N) | 25.2 | +0.34 ns | -0.29, 0.98 | 0.274 | +0.19 | 0.05 |
| AR - Keiser (35.7 N) | 23.4 | +0.39 ns | -0.11, 0.89 | 0.118 | +0.30 | 0.10 |
| MO - Columbia (39.0 N) | 21.6 | +0.51 ns | -0.07, 1.08 | 0.081 | +0.53 | 0.13 |
| IA - Ames (42.0 N) | 19.7 | +0.49 ns | -0.20, 1.17 | 0.154 | +0.54 | 0.09 |
| MN - Waseca (44.1 N) | 18.5 | +0.28 ns | -0.35, 0.92 | 0.368 | +0.31 | 0.04 |

Mean warming across the five sites: **+0.40 C/decade** (0 of 5 sites significant at p<0.05).

## Critical analysis

- **Direction & consistency.** All five sites show positive June-July night-time warming (+0.28 to +0.51 C/decade; mean +0.40), spanning ~14 degrees of latitude from the Texas Blacklands to southern Minnesota. The coherent positive direction -- 5 of 5 sites, sign-test p ~ 0.06 -- is the point most useful to a proposal: high night-time temperature is a maize-belt-wide exposure, not a single-site or Mid-South-only artefact.
- **Significance / statistical power.** Over the 25-year hourly record, none of the individual site slopes reach p<0.05 (p = 0.08-0.37; Columbia is closest). Interannual variability is large relative to a ~0.3-0.5 C/decade trend, so 25 annual values give limited power and wide confidence intervals. This is expected and is itself part of the motivation: the reanalysis-era record is too short to resolve site-level significance, which the proposed longer-term / network observations are designed to address. Theil-Sen slopes are reported alongside OLS and agree in sign at every site, indicating the positive tendency is not driven by a few outlier years.
- **Magnitude in context.** These are trends in the *dark-period mean*, the physiologically relevant quantity for night respiration, and are broadly consistent with the ~0.24 C/decade June-July Tmin trend reported for Arkansas over 40 years (paper Fig. S4); the somewhat larger per-decade values here partly reflect the shorter, more variable 25-year window.
- **Data caveats.** NASA POWER hourly is MERRA-2 reanalysis at ~0.5 deg (~50 km); a single grid cell represents the landscape, not a specific field, and can miss local radiative cooling on calm, clear nights. Absolute values carry reanalysis bias, but decadal *trends* are less sensitive to a stable bias. The hourly record begins in 2001, which sets the 25-year window.
- **Autocorrelation & attribution.** OLS treats years as independent; mild serial correlation would modestly inflate significance. These are descriptive climate trends, not an attribution analysis -- they motivate the proposal (night warming is broad and ongoing) rather than establish a yield-causal claim.
- **Site representativeness.** Sites were chosen as recognizable maize/AES locations spanning the latitudinal gradient; results are illustrative of the transect, and a production-weighted or multi-cell average would be the next step for a formal regional estimate.

## Outputs
- `historical_nt_transect.png`       -- main proposal figure (overlaid trends)
- `historical_nt_transect_facet.png` -- per-site detail panels
- `historical_nt_transect.csv`       -- site x year night-time temperature
- `historical_nt_transect_slopes.csv`-- per-site slope statistics
