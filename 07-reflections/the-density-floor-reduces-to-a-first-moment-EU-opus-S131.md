# The density floor reduces to a first moment: μ_{1/7} ≥ E[uncovered length]

**opus-2026-07-07-S131.** A clean analytic reduction of Route 1's density-floor positivity, found
while working the μ_{1/7} node. It turns the hard *tail-probability* floor into a *first-moment*
bound, which is more amenable to the moment method.

## The reduction

Let `U(x) = ` uncovered length `= Σ_{gaps} (gap − 1/7)_+` (a gap of length `g > 1/7` leaves an
uncovered sub-arc of length `g − 1/7` after the radius-`1/14` arcs cover `1/7` of it). Then
`U(x) ∈ [0,1]` and `μ_{1/7}(E) = meas{x : maxgap > 1/7} = meas{x : U(x) > 0}`.

**Paley–Zygmund.** For `U ≥ 0`, `P[U > 0] ≥ (E[U])² / E[U²]`. Since `U ≤ 1`, `E[U²] ≤ E[U]`, so

> **μ_{1/7}(E) ≥ E[U] = E_x[ Σ_{gaps} (gap − 1/7)_+ ].**

So the density floor `μ_{1/7}(E) > 0` (for all clusters `E`, `k ≤ 13`) follows from
`inf_E E[U](E) > 0` — a lower bound on an **expectation** (first moment) rather than on a tail.

## What E[U] looks like

Computed (grid, exact-in-the-limit), `E[U](E)` for `k = 8..13` at the AP: `0.250, 0.211, 0.183,
0.158, 0.141, 0.127`. Positive and comfortably bounded below. Two structural facts:

- **E[U] is NOT minimized at the AP** (unlike μ_{1/7}). The AP gives `E[U] = 0.127` at `k=13`, but
  "12-AP + outlier" and random spread give slightly *less* (`≈ 0.118`). So the E[U] route is a
  *structure-independent* path to positivity — a different lever than AP-minimality, and possibly a
  cleaner one (no extremal-configuration identification needed, just a uniform floor on the mean).
- **The mean has clean two-term structure.** `E[U] = 1 − E[covered]`, and inclusion–exclusion gives
  `E[covered] = k/7 − C(k,2)·(3/196) + (triples…)`, where the pairwise term `3/196` is
  **structure-independent** (`E_x[overlap of two radius-1/14 arcs at speeds e_i,e_j] = 3/196` for
  every pair, since `(e_i−e_j)x` is equidistributed). The structure enters only at triples and
  higher. The 2-term estimate `1 − k/7 + C(k,2)3/196 ≈ 0.34` (k=13) overshoots the true `0.127` — the
  triple+ terms are large — so a rigorous uniform lower bound on `E[U]` still needs to control the
  higher-order overlaps (the same wall the naive union bound hits).

## Status

- **Rigorous:** the reduction `μ_{1/7}(E) ≥ E[U]` (Paley–Zygmund) is a clean, correct inequality.
- **Empirical:** `E[U] ≈ 0.12–0.25 > 0`, bounded below (`≈ 0.118` at `k=13`) over structured +
  random clusters.
- **Open (the remaining lemma, now cleaner):** `inf_E E[U](E) > 0` uniformly for `k ≤ 13`. This is a
  first-moment bound; the obstruction is controlling the triple+ inclusion–exclusion overlaps
  (structure-dependent). It is a genuinely more tractable target than either the tail `μ_{1/7}` or the
  extremal AP-minimality — a mean, with two-thirds of its terms structure-free.

So Route 1's density-floor positivity now has **two** independent routes, both open at one uniform
lower-bound step: (i) AP-minimality of the tail (exact floor `477/1078`, hard extremal problem), and
(ii) `inf E[U] > 0` (first moment, the moment method). The second is the one to push analytically.
