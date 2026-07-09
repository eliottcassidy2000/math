---
id: LEM-012
title: The gap-splitting good-period bound — if a k-element co-offset set E contains an arithmetic progression of length L ≥ k−5, then a good period exists at j ≤ ⌈7(L−1)/(L−k+6)⌉ = O(k), by an ELEMENTARY argument (Dirichlet clustering of the sub-AP + a pigeonhole gap-split): cluster the L-term AP into a circular arc of span < (L−k+6)/7 so its complement is a single gap of length > (k−L+1)/7, which the remaining m=k−L ≤ 5 points can split into at most m+1 pieces, the largest > 1/7. This elementarily closes the near-AP branch (longest-AP ≥ k−5) of the general j*=O(k) capstone — including mac-mini's exact-AP case (L=k, bound = ⌈7(k−1)/6⌉) — leaving only the deeply-dissociated case (longest-AP ≤ k−6) to kps-S91; the two ranges tile all L with no gap
status: PROVED (elementary: Dirichlet's simultaneous-approximation theorem in ONE dimension + a pigeonhole on the split gap; no Weyl / equidistribution / resonance sum). VERIFIED: the Dirichlet-cluster j leaves a >1/7 gap in 100% of constructed long-AP clusters (L=k..k−5, m=0..5, k=11,12,13, Vmax=91/200/400), and j* ≤ Q in all sampled hard clusters with L ≥ k−5. Closes THM-527-A's finite check for the near-AP branch WITHOUT the 𝒲̂-resonance bound
source: klein-2026-07-09-S196
depends_on:
  - THM-527   # good period j ⟺ maxgap{frac(e_i j/Vmax)} > 1/7; the finite-Vmax glue THM-527-A
  - LEM-010   # the deterministic good-period lemma; this bounds its j* for the near-AP branch
related:
  - THM-664   # the Weyl route (the analytic alternative this sidesteps for L ≥ k−5)
  - LEM-011   # the 𝒲̂ route (needed only for the dissociated L ≤ k−6 complement, kps-S91)
external: Dirichlet's approximation theorem; the three-distance / pigeonhole.
---

# LEM-012 — The gap-splitting good-period bound (elementary, near-AP branch)

## Statement

Let `E = {e_0=0, e_1, …, e_{k−1}} ⊂ ℤ` (`k ≤ 13`) be a co-offset set, `V > max E` the ruler. A
**good period** is `j ∈ {1,…,V−1}` with `maxgap{ e·j mod V : e∈E } > V/7` (THM-527: `⟹ M(S) ≥ 1/14`).
Suppose `E` contains an **arithmetic progression of length `L`** (common difference `d`), i.e.
`{a, a+d, …, a+(L−1)d} ⊆ E`, with

> **`L ≥ k − 5`.**   Then a good period exists at some `1 ≤ j ≤ Q := ⌈7(L−1)/(L−k+6)⌉ = O(k)`.

(`Q ≤ 49` for all `k ≤ 13`; for `L = k`, `Q = ⌈7(k−1)/6⌉` — mac-mini-S59's exact-AP bound, recovered.)

## Proof (elementary)

Write `m := k − L ≤ 5` (the number of points of `E` off the AP) and `S := (L−k+6)/7 = 1 − (m+1)/7 > 0`
(the hypothesis `L ≥ k−5` is exactly `S > 0`).

**Step 1 — cluster the sub-AP (Dirichlet).** By Dirichlet's approximation theorem, for the integer
`Q = ⌈(L−1)/S⌉ = ⌈7(L−1)/(L−k+6)⌉` there is `j ∈ {1,…,Q}` with `‖j d / V‖ < 1/Q ≤ S/(L−1)` (this is
UNCONDITIONAL in `d, V`). Put `s := (jd \bmod V)` in `(−V/2, V/2]`, so `|s| < V·S/(L−1)`. At this
dilation the AP maps to `{(a+id)j \bmod V} = {aj + i·s \bmod V : i=0,…,L−1}` — an AP of step `s`
spanning `|(L−1)s| < V·S < V`, so it lies in a **single circular arc `A` of length `< S·V`** (no wrap).

**Step 2 — the complement is one big gap.** Since the `L` AP-points sit inside the arc `A` of length
`< S·V`, the complementary arc `Γ := (circle) ∖ A` is a single gap of length `> (1−S)·V = (m+1)V/7`,
containing NO AP-point.

**Step 3 — the `m` stray points cannot fill it (pigeonhole).** The remaining `m = k−L` points of `E`,
at this dilation, sit somewhere on the circle. Only those falling inside `Γ` subdivide it; there are at
most `m` of them, so they cut `Γ` into at most `m+1` sub-arcs whose lengths sum to `|Γ| > (m+1)V/7`.
Hence the **largest sub-arc `> |Γ|/(m+1) > V/7`** — a gap `> V/7`. So `maxgap{e·j mod V} > V/7`: `j` is a
good period, and `j ≤ Q`. ∎

(Coincidences only help: if stray points merge with each other or with AP-points, there are fewer
cuts and the surviving gap is only larger.)

## Verification

Machine-checked (`lrc14_gapsplit_mechanism_klein_S196`): over constructed clusters `E = AP_L ∪ {m
random points}` for every `L ∈ {k,…,k−5}` (`m = 0,…,5`), `k = 11,12,13`, `V = 91/200/400`, hard
(`spread ≥ 6V/7`, `j=1` fails) — the Dirichlet-cluster `j` leaves a `> 1/7` gap in **100%** of cases,
and `j* ≤ Q` always. Also (`lrc14_gapsplit_klein_S196`) over sampled hard clusters, `j* ≤ Q` holds
whenever `L ≥ k−5` (0 failures).

## Consequence — the near-AP branch of `j*=O(k)` is ELEMENTARY

The general `j*=O(k)` capstone (the last math gap of the covering case, THM-527-A / LEM-010) splits on
the **longest AP** `L(E)`:

- **`L ≥ k−5` (near-AP):** LEM-012 — `j* ≤ ⌈7(L−1)/(L−k+6)⌉`, **elementary** (Dirichlet + pigeonhole,
  no equidistribution). This is the *hard* branch — the near-AP is where `j*` is largest (`≈ k`), and
  where mac-mini had only the exact AP and opus needed the `r_N` (𝒲̂) verification. Now it is a
  one-page elementary proof.
- **`L ≤ k−6` (deeply dissociated):** kps-S91 — `j* ≤ 3` (verified; a-priori via the LEM-011
  𝒲̂-smallness, which holds because a short longest-AP means few resonances).

The two ranges **tile `{2,…,k}` with no gap** (`k−5` and `k−6` are consecutive). So `j*=O(k)` holds for
every cluster, closing THM-527-A's finite check to `V ≤ O(k)` (trivial) — **modulo only** the a-priori
`𝒲̂`-smallness of the dissociated branch (or an elementary argument for it). The near-AP branch — the
one everyone expected to need the analysis — is now elementary.

## Remaining (the dissociated branch)

`L ≤ k−6` (only reachable for `k ≥ 8`; `L ≤ 7` at `k=13` down to `L=2` Sidon at `k=8`). Empirically
`j* ≤ 3` (kps-S91). Two honest routes: (a) the 𝒲̂ partial-sum bound (LEM-011, the same object as the
density-floor tail) — `dissociated ⟹ few small-`n` resonances `n·E≡0` ⟹ the correction is small ⟹
`E_grid[W] ≈ (6/7)^k > 0` at small `N`); (b) an elementary gap-argument for short-longest-AP hard sets
(open). Either closes `j*=O(k)` and hence LRC(14)'s covering case outright.
