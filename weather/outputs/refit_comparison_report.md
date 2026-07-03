# Refit: dark-period night-time temperature vs Tmin

Same mixed model as the paper (`yieldkg.ha ~ temp + z_cum_radn + z_vpd + z_april + z_may + z_july + z_august + z_sept`,
random `~1 | location/hybrid`, `varIdent(~1 | location)`), refit with the temperature term swapped.
Both models fit on the same 4663 rows.

## Temperature effect and yield sensitivity

| Model | std beta | kg/ha per SD | SD (C) | kg/ha per C | % per +1C | p | AIC | RMSE (kg/ha) |
|---|---|---|---|---|---|---|---|---|
| mean_mint (Tmin, paper) | -0.327 | -545 | 1.08 | -503 | -3.66% | 1.1e-71 | 80373.3 | 1356 |
| nightT (dark-period) | -0.437 | -728 | 1.14 | -641 | -4.66% | 8e-85 | 80302.5 | 1345 |

## Sensitivity across scenarios (% yield change)

| +C | Tmin | nightT |
|---|---|---|
| +1 | -3.7% | -4.7% |
| +2 | -7.4% | -9.4% |
| +3 | -11.0% | -14.1% |
| +4 | -14.7% | -18.8% |

## Read-out

- The yield penalty per +1 C is **-3.7% using Tmin vs -4.7% using the dark-period night-time temperature** -- the negative effect is preserved and, if anything, slightly larger with the physically-correct metric.
- The dark-period model also fits marginally better (AIC 80302.5 vs 80373.3; RMSE 1345 vs 1356 kg/ha), and the temperature effect is at least as significant (p = 8e-85 vs 1e-71).
- Because nightT is ~1:1 collinear with Tmin (r ~ 0.98, see the main report), the conclusion does not depend on the metric. This is the robustness result for Reviewer 1: the negative NT-yield relationship is not an artefact of using Tmin, and holds (slightly strengthened) with the actual dark-period temperature the crop experiences.

Outputs: `refit_comparison.csv`, `refit_sensitivity.csv`, `refit_sensitivity.png`.

Note: fit with ML (not REML) so AIC is comparable across the two mean structures.
