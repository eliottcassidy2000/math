---
source: opus-2026-07-11-S260
status: PARTIAL + correction of S259's rigor claim. Attempting to prove the Erdős–Turán discrepancy bound for
  the coprime core gives a clean structure (exact Fourier identity + the independent model coreCover ≈
  1−(6/7)^core, CONFIRMED, margin (6/7)^core) BUT the naive Erdős–Turán / L² bounds are ~100–700× too weak,
  because the good set G' has large total variation (N ≈ 341 boundary points) and small measure (~0.1). So the
  naive discrepancy does NOT prove Σε < 1/7. The refined path is MOLLIFICATION of G' (Beurling–Selberg / Fejér
  = LRCFourierCompletion) to cut the total variation, then Erdős–Turán — feasibility = a total-variation vs
  margin trade-off (open). Corrects S259's "rigor via Erdős–Turán".
tags:
  - lrc14
  - covering-min
  - anti-concentration
  - erdos-turan
  - discrepancy
  - mollification
  - correction
---

# The naive Erdős–Turán is too weak; the good set needs mollification

**opus-2026-07-11-S260.** Owner: prove the Erdős–Turán discrepancy bound for the coprime core (S259's rigor
step). Doing so gives clean structure and a confirmed quantitative model — but shows the *naive* Erdős–Turán is
far too weak, correcting S259's optimism.

## The exact structure (clean)

For a core runner `v` and the good set `G'`, the exact Fourier identity is

> `|D_v ∩ G'| − (1/7)|G'| = Σ_{h≠0} b_h ĝ(−hv)`,  `b_h = sin(πh/7)/(πh)`, `ĝ = ` Fourier coeffs of `1_{G'}`.

The **Markov reduction** is then clean: `coreCover ≤ E[W_core | G'] = Σ_core density(D_v in G') = 6/7 + Σε`,
so `coreCover < 1` follows from **`Σε < 1/7`** (where `ε_v = density − 1/7 = (Σ_{h≠0} b_h ĝ(−hv))/|G'|`).

## The independent model — confirmed (the positive)

Under independence (`D_v` decorrelated from `G'` by coprimality), inclusion–exclusion gives
`coreCover ≈ 1 − (6/7)^{core} < 1`, with margin `(6/7)^{core}`. **Verified against the data:**

| `|core|` | mean coreCover | `1−(6/7)^{core}` |
|---|---|---|
| 1 | 0.148 | 0.143 |
| 2 | 0.282 | 0.265 |
| 3 | 0.398 | 0.370 |
| 4 | 0.483 | 0.460 |
| 5 | 0.579 | 0.537 |

The match is tight (actual slightly above, a small positive correlation). So the core **is** nearly independent
in `G'`, `coreCover ≈ 1 − (6/7)^{core} < 1` with a substantial margin `(6/7)^{core} ≈ 0.4–0.85`. This is the
right quantitative picture, and the actual `|ε_v|` is small (mean `0.02`, max `0.086`).

## The naive Erdős–Turán is ~100–700× too weak (the correction)

The naive bound `|ε_v| ≤ N/(6v|G'|)` (from `|ĝ(m)| ≤ N/(2π|m|)` and `Σ|b_h|/|h| ≤ π/6`) is useless here:
`N ≈ 341` boundary points of `G'` (mean; max `532`), `|G'| ≈ 0.1`, so for `v = 41` the bound is
`341/(6·41·0.1) ≈ 14` — versus the actual `ε ≈ 0.02`. **~700× too weak.** The L²/Cauchy–Schwarz variants are
also too weak (they scale like `1/√|G'|` or `N/√v`, all `≫ 1/7`). The culprit is the **large total variation of
`1_{G'}`** (hundreds of intervals, from the many non-core runners) combined with **small `|G'|`**. So the naive
discrepancy machinery does **not** prove `Σε < 1/7`, even though the equidistribution is empirically robust.
S259's "rigor via Erdős–Turán" was too optimistic: the crude bound fails; the truth needs the **cancellation**
among the `ĝ(−hv)` (which coprimality makes small *on average*, not individually).

## The refined path: mollify the good set

The obstruction is `1_{G'}`'s total variation. The fix is **mollification** (Beurling–Selberg / Fejér
smoothing — exactly the project's **LRCFourierCompletion** completion identity `|C_w − b²/q|`): replace `1_{G'}`
by a smoothed minorant `g̃ ≤ 1_{G'}` with **rapidly-decaying** Fourier coefficients, so `Σ_h b_h ĝ̃(−hv)`
converges fast and is small. The cost is the mollification error `|G' \ {g̃ = 1}| ~ N·δ` (boundary smeared by
width `δ`); optimizing `δ` (Beurling–Selberg) trades total variation against decay. Whether the optimized error
beats the margin `(6/7)^{core}` is the open quantitative question — but it is the *correct* analytic route, and
it connects to the fleet's existing Fourier-completion machinery (tasks B.1–B.3).

## Net (honest)

Proving the discrepancy bound: the **structure is clean** (Fourier identity + Markov reduction + the confirmed
independent model `coreCover ≈ 1 − (6/7)^{core} < 1`, margin `(6/7)^{core}`), but the **naive Erdős–Turán is
~700× too weak** (large total variation of `G'`), so it does not close `Σε < 1/7`. This corrects S259's rigor
claim. The equidistribution is real (empirically + the independent model), and the honest remaining step is a
**mollified (Beurling–Selberg / Fejér) discrepancy estimate** — the LRCFourierCompletion direction — whose
success hinges on the total-variation-vs-margin trade-off. So the crux is now a *specific analytic inequality*
(mollified discrepancy of the coprime core against `G'` `< (6/7)^{core}`), not a vague "anti-concentration."

→ opus-S259 (mechanism — rigor claim corrected here), LRCFourierCompletion / tasks B.1–B.3 (the mollification
machinery), opus-S258 (the tools that fail), opus-S255 (runner-1 / near-AP still separate), s558o. Files:
`lrc14_six_core_equidistribution_union_bound_opus_S259.py` (S259 base; this session is analysis — inline
computations: independent-model match, naive-bound weakness, boundary-count `N ≈ 341`).
