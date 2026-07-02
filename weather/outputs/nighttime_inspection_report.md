# Dark-period night-time temperature — inspection report

- Source: NASA POWER **hourly** T2M (community=AG), same product as the paper, native hourly resolution.
- Window: June 1 – July 31 (DOY 152–212), matching the paper's fixed window.
- Night definition: geometric horizon (elev < 0°), via NOAA solar-position algorithm (UTC).
- Site-years processed: **64**
- Mean dark hours/day: 9.5; mean missing hourly values/site-year: 0.0

## Headline

- **Correlation (mean_mint vs dark-period nightT): r = 0.983** (n = 63)
- OLS: nightT = 1.010·mean_mint + 1.48
- Mean offset (nightT − mean_mint): **+1.70 °C** (dark-period mean runs warmer than the daily minimum)
- RMSE(nightT − mean_mint): 1.71 °C
- For contrast, correlation of mean_mint with **day-time** T2M: r = 0.881
- Mean day−night difference (diurnal range proxy): 5.84 °C

## QC bridge

- Reconstructing a Tmin-like metric from the hourly pull (mean of daily hourly minima) vs the paper's `mean_mint`: r = 0.998, RMSE = 0.13 °C. High agreement confirms the hourly pull is consistent with the daily-based `mean_mint`.

## Interpretation

- A very high r means the dark-period mean and the daily minimum carry nearly the same site-year signal, so substituting one for the other will move the yield-model coefficient little — i.e. the choice of metric is not driving the result. That is the robustness point Reviewer 1's concern calls for: report r, the offset, and the refit coefficient side by side.
- The positive offset (dark-period mean above the single coldest instant) is expected and quantifies how much `mean_mint` under-states the temperature the crop actually experiences through the night.

## Per-location summary

```
               n  r(mint,nightT)  mean_mint  mean_nightT  night-mint
location                                                            
Des_Arc     10.0           0.983     21.950       23.617       1.667
Harrisburg   3.0             NaN     21.300       23.321       1.798
Keiser      13.0           0.960     21.808       23.474       1.666
Marianna    12.0           0.982     22.550       24.226       1.676
Rohwer      14.0           0.987     22.743       24.461       1.718
Stuttgart   12.0           0.984     22.500       24.229       1.729
```

## Caveats

- NASA POWER is a ~0.5° reanalysis; absolute night-time values can miss local radiative cooling on calm, clear nights. The concern R1 raised is about *which hours* are averaged (dark hours vs the single minimum), which hourly data answers directly.
- Night defined by solar geometry; the 2 W m⁻² / 5 µmol threshold R1 mentions converges to nearly the same hours over a multi-week mean (set NIGHT_ELEV_DEG=-6 for a civil-twilight check).
- Window is the paper's fixed June–July; if the flowering window moves to a crop-calendar basis (Reviewer 2), recompute both metrics on that window — the script keys off WINDOW_*_DOY.