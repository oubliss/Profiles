# Bias correction sources

The coefficients in `profiles/bias.py` were transcribed from these. They are
kept so the transcription can be checked and so the provenance recorded in
output files points somewhere real.

| file | model | registered as |
|---|---|---|
| `rh_bias_correction_poly22.txt` | `poly22` over RH and T | `rh_poly22_uncorrected_t`, `rh_poly22_corrected_t` |
| `rh_bias_correction_poly41.png` | `poly41` (`sf_H_general`) | `rh_poly41_general` |

Both are MATLAB surface fits, so a coefficient `pij` multiplies
`rh**i * temp**j` and the values can be pasted straight out of a fit report.

Neither is applied automatically. `utils.rh_calib` is still a pass-through —
turning a correction on changes published values and is a decision for the
lab, not a default.
