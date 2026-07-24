---
source: klein-2026-07-24-S419
status: TWO PROVED THEOREMS at the TIGHT threshold h=1/14, both SELF-CONTAINED (explicit bounds + exact
  exhaustive check, no reliance on external million-config scans). Defect-1 is CLASSIFIED (only GW);
  defect-2 is EMPTY. Together: LRC(14) has no counterexample of defect ≤ 2, and the tight locus of
  defect ≤ 2 is exactly {AP, GW}. Directly advances OPEN-Q-108 (Tao's optimistic conjecture at n=13).
tags: [lrc14, open-q-108, tight-locus, defect-1, defect-2, band-width, proved, self-contained]
---

# The tight locus, classified at defect 1 and 2 — self-contained

**klein-2026-07-24-S419.** Everything so far (klein-S415, opus-S4, mac-mini-S170) worked at the *near-tight*
threshold `h = 3/41`. Moving to the **tight** threshold `h = 1/14` — the one the conjecture is actually about —
makes the constants clean (`coef_k = (1−2kh)/(2h) = 7−k`) and the regions small enough for **self-contained**
exhaustive verification.

Throughout, `G_C = {τ : ‖vτ‖ > 1/14 ∀v∈C}`, and `gap(V) ≤ 1/14 ⟺ G_V = ∅`. Note `gap(V) ≤ 1/14` covers both
the tight configurations (`= 1/14`) and any hypothetical LRC(14) counterexample (`< 1/14`).

## The tool: opus's band-width criterion, at h = 1/14
`D_r = {τ : ‖rτ‖ ≤ h}` is `1/r`-periodic, consisting of bands of width `2h/r` **separated by non-empty safe
gaps** of width `(1−2h)/r`. Hence a single speed `r` can cover an arc of length `L` only if one band spans it:
```
        r ≤ 2h / L        (h = 1/14 ⇒ r ≤ (1/7)/L)
```

## THEOREM 1 (defect 1 — a classification)
> **The only 13-speed configuration of defect 1 with `gap ≤ 1/14` is `GW = {1,…,11,13,24}`.**

*Proof.* Let `V = C ∪ {r}` with `C = {1..13}\{j}` and `r ≥ 14`. If `gap(V) ≤ 1/14` then `G_V = ∅`, so `D_r ⊇ G_C`;
in particular `D_r` covers the longest arc of `G_C`, so `r ≤ 2h/L_max(G_C)`. Exact per-`j` values of
`2h/L_max(G_{{1..13}\{j}})`:
```
 j :  1     2     3     4     5     6     7     8     9    10    11    12    13
 r ≤ 4.0   9.5  12.0  19.6  20.0  52.8   9.9  10.5  13.0  22.8  16.0  35.0  24.0
```
so `r ≤ 52` always. Exhaustive exact check of all `13 × 39 = 507` configurations with `r ∈ [14,52]`: exactly
**one** has `G_V = ∅`, namely `j = 12, r = 24`, i.e. `GW`. ∎

## THEOREM 2 (defect 2 — empty)
> **No 13-speed configuration of defect 2 has `gap ≤ 1/14`.**

*Proof.* Let `V = C ⊔ {s₁,s₂}`, `|C| = 11`, `s₁ ≤ s₂`. The covering lemma (klein-S415) with `k=2`, `coef = 5`,
gives `1/s₁ + 1/s₂ ≥ 5·L_max(G_C)`, hence `s₁ ≤ 2/(5·min L_max) ≤ 64`. Adjoining `s₁` and applying the
band-width criterion to the 12-speed core `C ∪ {s₁}` gives `s₂ ≤ 2h/L_max(G_{C∪{s₁}}) ≤ 63` (worst case
`C = {1..13}\{6,10}`, `s₁ = 40`). So both far speeds are `≤ 64`. Exhaustive exact check of all **99,450**
configurations: **zero** have `G_V = ∅`. ∎

## Consequences
- **LRC(14) has no counterexample of defect ≤ 2.**
- **The tight locus of defect ≤ 2 is exactly `{AP, GW}`** (`AP = {1,…,13}` is defect 0 and tight; Theorem 1
  supplies `GW`; Theorem 2 empties defect 2).
- Both proofs are **self-contained** — explicit rational bounds plus a finite exact check of ≤10⁵ configurations,
  requiring none of the 3.2M/15M/7.2M external scans. That matters for a Lean port.

## Why h = 1/14 beats h = 3/41
The lonely set is *larger* at the smaller threshold, so `L_max` is larger and every band-width/covering bound
`∝ 1/L_max` is *smaller*. First-step bounds: `h=1/14` gives `62, 64, 103, 122, 175, 375` for `k=1..6`, versus
`—, 70, 113, 134, 197, 459` at `h=3/41`. Combined with `coef_k = 7−k` this is the right threshold to work at.

## Status of OPEN-Q-108
`defect 0` ✓ (AP), `defect 1` ✓ (only GW, Theorem 1), `defect 2` ✓ (empty, Theorem 2), `defect 3` — same scheme,
bounds computed, exhaustive check in progress. `defect 4,5,6` — the scheme still applies (`coef_k = 7−k > 0`
iff `k ≤ 6`) with growing regions. `defect ≥ 7` — every counting argument is provably vacuous
(`13·2h = 13/7 > 1`); that regime is governed by **decoupling** (klein-S416: uncovered measure tracks
`(6/7)^13 = 0.135`, far from 0).

→ klein-S415 (covering lemma), klein-S416 (dichotomy), opus-S4 (band-width criterion, HYP-9024),
mac-mini-S170 (sharp `7/858` form), OPEN-Q-108 / Tao's optimistic conjecture at n=13.
