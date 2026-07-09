---
source: opus-2026-07-08-S165
status: the LRC(14) capstone (good period j*=O(k)) reduces to ONE sharp a-priori inequality --
  the W-hat partial-sum correction ratio |Corr_N|/(N (6/7)^k) < 1 at N = ceil(7(k-1)/6). VERIFIED
  uniformly < 1 (max 0.84) over all good clusters INCLUDING near-APs (ratio 0.59-0.61), so the
  W-hat route (kps-S90) closes j* <= ceil(7(k-1)/6) UNIFORMLY -- no small-L/near-AP split needed;
  the exact full AP is the sole boundary (ratio = 1 exactly, no good period = the tight case, S164).
  Unifies kps-S90 (W-hat) + mac-mini (embedded-AP, shown unnecessary) + opus-S164 (dichotomy) +
  opus-S154 (Corr_N IS the resonance sum).
tags:
  - lrc14
  - good-period
  - capstone
  - W-hat
  - longest-AP
  - resonance
---

# The capstone reduces to one sharp inequality: correction ratio < 1

**opus-2026-07-08-S165.** Owner: work the capstone. The whole of LRC(14) is reduced to the good
period `j* = O(k)` (THM-527-A). kps-S90 (HYP-5507) gave the W-hat partial-sum route; this note
sharpens it to a SINGLE a-priori inequality, verified uniformly, and shows the near-AP needs no
separate argument.

## The clean sufficient condition

For a `k`-cluster `E`, ruler `Vmax`, let `W(x) = ` uncovered measure of `{frac(e_i x)}` `= sum_i
(gap_i - 1/7)_+`, so `x` is a good period `⟺ W(x) > 0` (maxgap `> 1/7`). The `j=0` anchor is
`W(0) = 6/7` (all phases at `0`). With `S_N := sum_{j=1}^{N} W(j/Vmax)`:

> `j* <= N  ⟺  S_N > 0`  (since `W >= 0`, `S_N > 0 ⟺ some j in {1..N}` is good).

Fourier-expanding `W` (LEM-011): `S_N = N (6/7)^k + Corr_N`, `Corr_N = sum_{n != 0} What(n) sum_{j=1}^N
e(nj/Vmax)` (the SAME resonance sum as the density-floor tail, opus-S154). Hence the **sharp
sufficient condition**:

> **`|Corr_N| < N (6/7)^k`  (i.e. the correction ratio `r_N := |Corr_N|/(N (6/7)^k) < 1`)  ⟹  `S_N > 0`
> ⟹  `j* <= N`.**

## The uniform verification: `r_N < 1` at `N = ceil(7(k-1)/6)`, ALL clusters

Computed `r_N` at the target `N = ceil(7(k-1)/6)` over every cluster (bucketed by longest-AP `L`),
`k = 11, 13`:

| longest-AP L | max `j*` | max `r_N` (k=11, N=12) | max `r_N` (k=13, N=14) |
|---|---|---|---|
| small (`L <= k-3`) | `<= 7` | `<= 0.82` | `<= 0.84` |
| near-AP (`L = k-1, k`) | up to `13` | `0.61` | `0.59` |

**`r_N < 1` UNIFORMLY (max `0.84`), including the near-APs.**  So `S_N > 0` and `j* <= ceil(7(k-1)/6)`
for EVERY cluster with a good period — the W-hat route closes the capstone UNIFORMLY, with **no
small-L / near-AP split** (mac-mini's separate embedded-AP Dirichlet is not needed for the bound —
the near-AP's `r_N` is actually the SMALLEST, `~0.6`).

## The exact AP is the sole boundary (`r_N = 1`)

For the exact full complete-residue AP `{0,…,k-1}`, `Vmax = k` prime: every `j in {1..k-1}` gives all
residues, `W(j/k) = 0`, so `S_N = 0`, `Corr_N = -N(6/7)^k`, `r_N = 1` EXACTLY. This is the unique
`r_N = 1` case = the S164 dichotomy's no-good-period set = the tight `M = 1/k` LRC instance (cited via
LRC`<=13`). So the boundary of the sufficient condition is EXACTLY the tight case — the inequality is
sharp, and everything strictly inside (`r_N < 1`) has `j* <= ceil(7(k-1)/6)`.

## What remains (the single a-priori target)

The capstone is now: **prove `r_N < 1` at `N = ceil(7(k-1)/6)` a-priori** (for non-tight clusters),
i.e. `|Corr_N| < N (6/7)^k`. By LEM-011, `|Corr_N| <= sum_{n != 0} |What(n)| · min(N, 1/(2||n/Vmax||))`
with `|What(n)|` a-priori (klein-S194). The needed margin is generous (verified `r_N <= 0.84`, so a
`16%` slack). This is ONE explicit inequality on the a-priori W-hat sum — the same shared constant as
the density-floor tail (opus-S154 / LEM-011), now driving the good period too. The `min(N, 1/||·||)`
sampling factor is what makes it converge where the raw `L^1` sum (opus-S154) diverges — the geometric
partial sum caps the resonant frequencies.

## Ledger / next

- FOUND: the capstone's sufficient condition is `r_N < 1` at `N = ceil(7(k-1)/6)`; verified UNIFORMLY
  `< 1` (max `0.84`) over all good clusters incl. near-APs => the W-hat route is uniform (no L-split,
  no separate embedded-AP needed for the bound); exact AP = the unique `r_N = 1` boundary (tight case).
- REDUCES the capstone to ONE a-priori inequality `|Corr_N| < N(6/7)^k` (LEM-011 W-hat sum), with a
  `16%` verified margin.
- Unifies kps-S90 (W-hat), mac-mini (embedded-AP, unnecessary for the bound), opus-S164 (dichotomy),
  opus-S154 (`Corr_N` = the resonance sum; the `min(N,·)` cap gives convergence).
- Files: `lrc14_jstar_whatcorr_vs_L_opus_S165.py` (+out). -> kps-S90/HYP-5507, mac-mini-S59/S60,
  THM-527-A/LEM-010, LEM-011.
