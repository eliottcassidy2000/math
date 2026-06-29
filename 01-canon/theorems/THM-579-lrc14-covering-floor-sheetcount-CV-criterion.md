---
id: THM-579
title: The LRC14 covering-floor decorrelation ratio obeys the CLOSED-FORM Cauchy-Schwarz bound R' >= 1 - CV(N_R)*sqrt((1-m_Q)/m_Q), where N_R is the 14-SHEET COUNT of the small part R and m_Q = meas(lonely(Q)); hence the covering floor R'>0 holds whenever CV(N_R)^2 < m_Q/(1-m_Q). This SEPARATES the floor into an R-only quantity (sheet-count coefficient of variation, on the finite 14-grid) and a Q-only quantity, and DERIVES the HYP-3140 sheet-count framing from first principles. Verified on the binding consec family r=2..6 and on even-heavy R.
status: PROVED (the inequality, via Cauchy-Schwarz + the exact projection identity sum_{N!=0}|chat(14N)|^2 = Var(N_R)/196). VERIFIED that the sufficient condition holds (bound>0) for the binding consec rows r=2..6 (bounds 0.39-0.69) and even-heavy R. The UNIFORM bound over ALL covering (R,Q) is OPEN (reduces to: bound CV(N_R)^2 above and m_Q below over the family).
source: mac-mini-2026-06-29-S2
depends_on:
  - HYP-3415   # the critical-path reduction: LRC14 = q-witness + LRC<=13 + covering floor R'>0
  - HYP-3129   # the SPEC resonance-lattice elementary certificate (R'>=0.642 numerically)
related:
  - HYP-3140   # fiber-PGF / sheet-count Rprime route -- this DERIVES that object from Cauchy-Schwarz
  - THM-576    # m_Q >= C(14-r,2)/91 lower bound (pairwise-avoidance cap)
  - THM-578    # same decorrelation-residual machinery (doublet R-tail); the "finiteness beats sharpness" reframe
  - OPEN-Q-108
---

# THM-579 — the covering floor as a sheet-count coefficient-of-variation bound

## Setup (kps HYP-3415 critical path)
LRC(14) = [q-witness, elementary] + [LRC<=13 induction, arXiv:2604.23906] + [the covering FLOOR].
For covering `S = R u 14Q` (R = 14-free small part, |Q|=r<=6), the floor is
> `R'(S) := meas(lonely(S)) / [ meas(lonely(R)) . meas(lonely(Q)) ] > 0`  uniformly.
Writing `m_R = meas(lonely(R))`, `m_Q = meas(lonely(Q))`, and Fourier-expanding,
`meas(lonely(S)) = product + SPEC`, `product = m_R m_Q`,
`SPEC = sum_{n!=0} chat(n) conj(ghat(n))`, where `chat` = Fourier coeffs of `1_{R-safe}`
and `ghat` = Fourier coeffs of `1_{Q-lonely}(14t)`.  So `R' = 1 + SPEC/product`.
The danger-comb coefficients are `ahat(k) = -sin(pi k/7)/(pi k)`, `|ahat(k)| <= 1/(pi|k|)`.

## The 14-sheet count
> `N_R(t) := #{ a in {0,...,13} : t + a/14 is R-safe }`,  `E[N_R] = 14 m_R`.

## The theorem (PROVED)
> **`R' >= 1 - CV(N_R) . sqrt( (1 - m_Q)/m_Q )`,   `CV(N_R) = sqrt(Var N_R)/E[N_R]`.**
> Hence the covering floor holds (`R' > 0`) whenever **`CV(N_R)^2 < m_Q/(1-m_Q)`**,
> equivalently `E[N_R^2] . (1 - m_Q) < E[N_R]^2`.

### Proof
`ghat` is supported ENTIRELY on `14Z` (Q-lonely(14t) has frequencies in `14.Lat(Q) ⊆ 14Z`), so
`SPEC = sum_{N!=0} chat(14N) conj(ghat(14N))`.  Cauchy-Schwarz:
`|SPEC| <= sqrt( sum_{N!=0}|chat(14N)|^2 ) . sqrt( sum_{N!=0}|ghat(14N)|^2 )`.

Projection identity: the projection of a function `f` onto `14Z`-frequencies is
`P_14 f(t) = (1/14) sum_{a=0}^{13} f(t+a/14)`, so for `f = 1_{R-safe}`,
`P_14 f = N_R/14`.  By Parseval, `sum_{N}|chat(14N)|^2 = ||P_14 f||^2 = E[(N_R/14)^2]`, and
subtracting the DC term `|chat(0)|^2 = m_R^2 = (E[N_R]/14)^2`,
> `sum_{N!=0}|chat(14N)|^2 = (E[N_R^2] - E[N_R]^2)/196 = Var(N_R)/196`.
Since `ghat` is supported on `14Z`, `sum_{N!=0}|ghat(14N)|^2 = sum_{n!=0}|ghat(n)|^2 = m_Q - m_Q^2`
(Parseval; `1_{Q-lonely}` is an indicator so its square is itself).  Therefore
`|SPEC| <= sqrt(Var(N_R)/196) . sqrt(m_Q(1-m_Q)) = (1/14) sqrt(Var N_R) sqrt(m_Q(1-m_Q))`, and
`|SPEC|/product = |SPEC|/(m_R m_Q) = CV(N_R) . sqrt((1-m_Q)/m_Q)`.  So `R' = 1 + SPEC/product
>= 1 - |SPEC|/product >= 1 - CV(N_R) sqrt((1-m_Q)/m_Q)`. QED

## Verification (script: lrc14_floor_CV_sheetcount_bound_macmini_20260629.py)
| (R, Q) | m_R | m_Q | R' actual | CV^2 | bound | >0 |
|---|---|---|---|---|---|---|
| {1..12}, {1,2} (r=2) | .0341 | .7857 | .7015 | 1.0946 | **.4536** | YES |
| {1..11}, {1,2,3} (r=3) | .0563 | .6905 | .9673 | .8247 | **.3920** | YES |
| {1..10}, {1..4} (r=4) | .1380 | .6190 | 1.0579 | .1591 | **.6871** | YES |
| {1..9}, {1..5} (r=5) | .1811 | .5048 | 1.0396 | .1323 | **.6397** | YES |
| {1..8}, {1..6} (r=6) | .2658 | .4571 | .9587 | .1163 | **.6284** | YES |
| even-heavy R, {1,2} | .3381 | .7857 | .8118 | .3063 | **.7110** | YES |

The bound is positive on every binding consec row (so the floor is PROVED for them by this
clean criterion, with NO exact-low computation) and on even-heavy R (the S259 2-adic worry:
even-heavy R has LARGER m_R, hence SMALLER CV, hence a HIGHER bound -- the sheet-count framing
does not single out even speeds as worse).

## Significance
1. **Closed-form sufficient floor criterion.** Replaces HYP-3129's per-row "exact-low + tail"
   numerical certificate by a single inequality in two SEPARATED quantities: `CV(N_R)` (an R-only
   property of the 14-grid sheet count) and `m_Q` (a Q-only measure, `>= C(14-r,2)/91` by THM-576).
2. **Derives the sheet-count.** HYP-3140 reached for the fiber-PGF / 14-sheet-count `N_R` as the
   right object; here it FALLS OUT of one Cauchy-Schwarz step (the `14Z` projection of R-safe).
3. **Same machinery as the cap and the doublet.** `ahat(k)=-sin(pi k/7)/(pi k)` is the comb
   coefficient of THM-578; the floor SPEC, the cap inclusion-exclusion (THM-576), and the doublet
   R-tail are one family of resonance sums.
4. **Reduces the critical path** to: prove `CV(N_R)^2 < m_Q/(1-m_Q)` uniformly over covering `(R,Q)`
   -- a finite, separated trade-off (small r => large m_Q, small CV; large r => the reverse).

## Uniformity sweep (script: lrc14_floor_CV_uniformity_sweep_macmini_20260629.py)
Per-`r` requirement (R, Q independent parts): `max_{|R|=13-r} CV(N_R)^2 < cap_r/(1-cap_r)`,
`cap_r = min_{|Q|=r} meas(lonely(Q))`.  Sampling ~500 diverse 14-free R per r (consec, shifted,
AP, {1}+top, and random/large-speed up to ~70):

| r | \|R\| | cap_r | threshold cap/(1-cap) | max CV(N_R)^2 (sampled) | margin |
|---|---|---|---|---|---|
| 2 | 11 | 66/91 = .7253 | 2.6400 | .9010 | **+1.739** |
| 3 | 10 | 55/91 = .6044 | 1.5278 | .7138 | **+0.814** |
| 4 | 9 | 1979/4004 = .4943 | 0.9773 | .4657 | **+0.512** |
| 5 | 8 | 2243/5880 = .3815 | 0.6167 | .3183 | **+0.299** |
| 6 | 7 | 3029/10780 = .2810 | 0.3908 | .2949 | **+0.096** |

The criterion HOLDS for every `r=2..6` over the sampled family.  Margins shrink monotonically;
since `r<=6` (the floor reduction bounds `|Q|<=6`), **r=6 is the binding case** (margin +0.096,
|R|=7).  So the C-S floor criterion is uniform over the sampled bounded-speed family, and the
entire remaining gap localizes to: prove `max_{|R|=7, 14-free} CV(N_R)^2 < cap_6/(1-cap_6) ≈ 0.391`
(a small-R optimization), plus an unbounded-speed discharge for R.

## Spectral form (HYP-3533, mac-mini-S3)
By the Fejer/projection identity `sum_{d=-13}^{13}(14-|d|)e(nd/14)=196.[14|n]`, the CV is EXACTLY
the resonant spectral energy: **`CV(N_R)^2 = (sum_{N!=0}|chat(14N)|^2)/m_R^2`**, so the criterion is
**`sum_{N!=0}|chat(14N)|^2 < m_R^2 . m_Q/(1-m_Q)`** -- the 2-adically-resonant energy of `R-safe`
(its mass on `14Z\{0}`) dominated by `m_R^2 . threshold`.  This localizes the whole `14=2.7`
structure on the `14Z` lattice.  Sub-binomiality (`Var<=14 m_R(1-m_R)`) would reduce the floor to
the ELEMENTARY `1 < m_R+m_Q+13 m_R m_Q` (holds at all caps), BUT it FAILS at even/7-heavy R
(`rho=Var/[14 m_R(1-m_R)]` up to ~2.5 -- the explicit S259 2-adic mechanism), compensated because
high-rho configs have high `m_R` (anti-correlated).  See HYP-3533.

## Open
Uniform closure: prove `CV(N_R)^2 < cap_r/(1-cap_r)` for every `r`, tightest at `r=6` (margin
+0.096); equivalently bound the resonant energy `sum_{N!=0}|chat(14N)|^2 < m_R^2 m_Q/(1-m_Q)` via
the singular series (same object as THM-576/THM-578).  This needs (a) an upper bound on `CV(N_R)^2` over `|R|=13-r` 14-free configs (the sweep
suggests it is bounded well below the threshold, max at large-speed R), and (b) discharge of
unbounded R-speeds (large speeds equidistribute off the 14-grid; expected to lower CV).  Where
Cauchy-Schwarz is too lossy, fall back to HYP-3129's exact-low refinement -- but no sampled row,
including r=6, needed it.

Script: `04-computation/lrc14_floor_CV_sheetcount_bound_macmini_20260629.py`;
output `05-knowledge/results/lrc14_floor_CV_sheetcount_bound_macmini_20260629.out`.
Companion comb-resonance framework: `lrc14_comb_resonance_inclexcl_macmini_20260629.py`.
