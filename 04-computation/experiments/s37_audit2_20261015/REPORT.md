# S37 second independent audit (2026-10-15)

Independent adversarial referee (subagent), audit of the S37 material added after the first audit, at commit
`a405c50b8a` (THM-4658 statements 6-7, the census to 2^37, the HYP-9282 refutation, the halving null, the general-N
certificate and margins, HYP-9284, HYP-9285, the companion-note analogy). Four later commits (c64f94d093, f8a478757e,
b56bef2502, e62673a9f7) were not audited; where a later commit already fixes a finding, the report says so. Saved by
the session from the referee's returned text; all referee code is in this folder (one `.out` per `.py`).

**Summary.** No mathematical error in THM-4658 statements 6-7, the 2^37 census, the HYP-9282 refutation or the five
conditional residues; eight errors of interpretation, number or wording.

## Errors

- **E1 (moderate) - the nested-interval "consistency test" cannot fail.** By Lemma B2 the finite numerator is
  `C_N = E_N + R_N`, so `C(2^37) − C(2^34)` is the data between 2^34 and 2^37; with the census energies 0.060, 0.080,
  0.037 Cauchy-Schwarz bounds it by 1.87e-4 (actual 3.75e-6), while the 2^34 interval is wider than the 2^37 one by
  1.076e-3 under (2, 0). Nesting follows from data already in hand (margin 8.9e-4); even energies exactly 2 on those
  octaves would break it only by 3.5e-7. A code check, not evidence for (9a) beyond 2^37; same for the 2^40 nesting.
- **E2 (minor) - holdout numbers.** Single-zero holdout (fit 2^22..2^34, test 2^34..2^37) is 0.2174, not 0.19
  (unchanged for any fit start 2^14..2^30); the 44-zero in-sample error is 0.0078 (0.0074 if fit through 2^37), not
  0.006; 12/28/44-zero holdouts 0.0753, 0.0205, 0.0186 reproduce exactly.
- **E3 (minor).** "About 1800 times the largest energy ever observed" should be about 1000 (1806.6/1.769 = 1021).
- **E4 (minor).** "Maximum over all dyadic X <= 2^36 is 1.769" fails at X = 2 (`E_2 = 13`, `e = 3.25`); true for
  `4 <= X <= 2^36` (maximum 1.768761 at 2^24).
- **E5 (minor) - zero sets.** "Re ρ > 1/2", "Re ρ > 0.65" omit the upper bound `Re ρ < 1`: D has 194 further zeros
  with `0 < Im < 3000`, all with `2 < Re < 3` (none with `1 <= Re <= 2`); their residues are zero, so the formula is
  unaffected. Counts 44/28/12 confirmed by the argument principle.
- **E6 (minor-moderate) - pending 2^40 results written as obtained** at the audited commit (reproduction check,
  3,000 spot checks from the earlier script, "stress-tested on 12,042 values" without a committed script). Partly
  ledgered as MISTAKE-620.
- **E7 (moderate, typing) - statement 7 used beyond its scope.** "Carry no finite-depth power-law modes and are
  remainders" (note section 0) and "no finite-depth mode exists to resist the suppression" under PROVED (HYP-9284) go
  beyond statement 7: building upward level by level gives excess eigen-distributions with eigenvalue 1, `β(Z_3) = 0`
  and any prescribed level-1 marginal for every `ρ` tried (1/2, 0.3, 0.75+5i, 1.2−2i, `s_0`); all cylinder equations
  at levels 0-6 hold to 5e-15 (`b1b`). Statement 7's finite-depth exclusion is correct; the readings are heuristic.
- **E8 (moderate) - "H is exactly what is left of S after the D-zero modes are removed".** H is `−L[S]`; `L` kills
  floor-free D-zero powers, but H is `L` applied to the remainder plus floor residues `O(N^(Re ρ − 1))` per mode.
  Ledgered as MISTAKE-620 (clock_duality_oct15) after the audited commit; confirmed independently.

## Imprecisions

- I1. On the 64-point grid every zero above Im ≈ 290 aliases (25 of the 28). A fine-grid control (1024 per octave,
  exact data 2^22..2^31) gives the same fitted residues within 3% for every zero with `Re >= 0.83` and a similar rms;
  the low-Re coefficients are unstable across grid and range, so "28 of 28 inside their intervals" says nothing for
  the 23 intervals containing 0. Adding the 16 zeros with `3000 < Im < 4600` cuts the fine-grid holdout from 2.2% to
  0.65%.
- I2. "Below every null replicate" is a rank among 12 (p ≈ 1/13 per level) and the 30-36/14-29 split was chosen
  after seeing the data; the size of the gap is the stronger evidence.
- I3. In sample, Collatz beats every null replicate on persistence (75 vs <= 49; 20 vs <= 11; 27 vs <= 15); the null
  reproduces only the out-of-sample behaviour.
- I4. "An AR(1) correction leaves 7-17 channels" mixes methods: AR(1) with ρ = 0.63-0.65 leaves 6-8.
- I5. "25-300 times" should be 25-318; "keeps 65-91% of its value" is an autocorrelation.
- I6. Pre-registration P3 quotes null window means 4.6-22.9; the final 12 replicates give 1.82-22.86.
- I7. "Certify a sign" should be "certify a nonzero residue"; the grid has 1,728 points, not 1,727.
- I8. THM-4658's AUDITED status predates statements 6-7; non-real ρ needs the conjugate-pair form in 4(b)/7;
  "nilpotent relative to the halving permutation" is loose; the committed resonance-scan JSON lacked the spectrum-check
  key; the Re ρ = 0 numerical search was not committed.
- I9. The 2^40 run that produced the evaluated output was restarted at 20:25:46 with a revised script; disclose.
- I10. Companion note: D is approximated by (not "becomes") `1 + z^q + z^(q−p)/3`; validity
  `|s| ln 2 |log_2 3 − p/q| << 1` (for `s_0`, convergent 84/53: 0.0095).
- I11. The octave-averaged ratio 2+√2 needs many octaves (2.6, 3.1, 2.0 on 2^10..2^20).
- I12. "Fitted value near the middle of the interval" compares two estimators from the same data.

## Verdicts (abridged)

CONFIRMED: statement 6 (all 24 characters mod 9 and 27, twisted splits, all N <= 2^20); statement 7 density form,
spectra (characteristic polynomial to 1e-61 for Q <= 3, spectra to 6e-12 for Q <= 5), multiplicities (geometric
multiplicity Q at ζ = ±1, no Jordan blocks), finite-depth consequence and its stated limitation; census to 2^37 by an
independent census (all 1,728 grid points, octave sums k = 20..36, all 2^20 block moments, S and H at 2^37);
HYP-9282 refutation numbers; halving null (definition, randomness, same code; own null replicates consistent);
2^37 intervals and margins (independent: 1806.86/1195.09/292.90 at 2^37, 227.60/144.53/36.895 at 2^34; the
session's are slightly more conservative); remainder identity and algebra; zero counts; 12/28/44-zero holdouts;
attribution of the jump to 0.6939+22.25i; Theorem B3's formula applies at every simple zero with Re > 1/2 (all zeros
simple; Rouché check per zero); five certified nonzero residues re-implemented; no false certificate at a zero with
Re > 1; pre-registration integrity (adcd7eed1f precedes every 2^40 output; prediction coefficients depend only on data
below 2^37, reproduced to 7e-11; criteria unedited). REFUTED as stated: the nested-interval test (E1); the readings
of statement 7 (E7); the single-zero holdout 0.19 (E2).

## Scripts in this folder

`b1_statement7_spectrum`, `b1b_infinite_depth_eigendistributions`, `b2_statement6_remainder_identity`,
`b3_census_independent` and `b3b_census_compare`, `b4_hyp9282_recompute`, `b5_halving_null_check`,
`b6_certificate_independent`, `b6b_session_script_at_re_gt_1`, `b6c_lemma_B2_direct`, `b7_D_zeros_count`,
`b8_explicit_formula_fit`, `b8b_fit_attribution_and_high_zeros` (each `.py` with its `.out`).

## Session response

All errors and imprecisions applied (MISTAKE-621 in `01-canon/MISTAKES.md`); the pre-registration file and
prediction JSON were left as committed, with an erratum section.
