---
id: THM-761
title: Multi-exception sheet covering bound — a scaled core plus up to six exceptional runners preserves the 1/14 threshold whenever the exceptions' bad-sheet budget is under c; every gcd stratum admitted; closes pull-card A1 and the r=1 case of A2
status: PROVED (elementary; LRC(<=13) enters only through the core margin M(P) >= 1/(14-r))
source: opus-2026-07-14-S299
depends_on:
  - THM-760   # the r=1 coprime case (with the stronger 1/2-1/(2c) clearance)
  - LRC(<=13) # settled external input accepted by project policy
related: [THM-531, THM-668, THM-724, THM-755, THM-757, THM-758, HYP-6780, HYP-6785, HYP-6820, HYP-6825]
verification: 04-computation/lrc14_multi_exception_sheet_bound_opus_S299.py
  (+ 05-knowledge/results/lrc14_multi_exception_sheet_bound_opus_S299.out)
---

# THM-761 — Multi-exception sheet covering bound

## Statement

Let `P` be a nonempty finite set of distinct positive integer speeds, let `c >= 2`,
and let `W = {w_1, ..., w_r}` (`r >= 1`) be distinct positive integers none of which
is divisible by `c` (so `cP` and `W` are disjoint and `V = cP ∪ W` has `|P| + r`
distinct speeds). Write `g_a = gcd(w_a, c)` — a proper divisor of `c` — and fix a
clearance target `delta` with `0 < delta <= 1/4`.

Let `t0` be any time with core margin `min_{p in P} ||p t0|| >= delta_core`, and
define the **sheets**

`t_k = (t0 + k)/c`,  `k = 0, 1, ..., c-1`.

**(i) Core exactness.** For every sheet, `min_{v in cP} ||v t_k|| = min_{p in P} ||p t0||`.
The core margin is preserved EXACTLY on every sheet — no measure, no continuity, no
raw-height dependence.

**(ii) Bad-sheet count.** Call sheet `k` *bad for* `w_a` if `||w_a t_k|| < delta`. Then

`#bad(w_a) <= g_a * (floor(2*delta*c/g_a) + 1) <= 2*delta*c + g_a`.

**(iii) Union bound / free sheet.** If

`sum_{a=1}^r g_a * (floor(2*delta*c/g_a) + 1) <= c - 1`,

then some sheet `t_k` is bad for no exception, and at that sheet

`min_{v in V} ||v t_k|| >= min(delta_core, delta)`.

Taking `t0` a maximin point of the core: **`M(cP ∪ W) >= min(M(P), delta)`**.

## The LRC(14) corollary

Take `|P| = 13 - r` (thirteen speeds total) and `delta = 1/14`, so the bad arc has
length `1/7` and `#bad(w_a) <= g_a*(floor(c/(7*g_a)) + 1) <= c/7 + g_a`. Settled
LRC(14-r) gives `M(P) >= 1/(14-r) > 1/14` for `1 <= r <= 6`. Hence:

> **If `sum_a g_a*(floor(c/(7*g_a)) + 1) <= c - 1`, then `M(cP ∪ W) >= 1/14`.**
> In particular (all exceptions coprime to `c`, `g_a = 1`):
> **any `r <= 6` exceptional runners dodge at every scale `c >= 43`**, and the exact
> per-`(r, c)` criterion `r*(floor(c/7) + 1) <= c - 1` is decidable by inspection.
> For `r = 1` the criterion holds for EVERY `c >= 2` and EVERY gcd stratum
> `g_1 | c`, `g_1 < c` — extending THM-760 to non-coprime exceptions.

Exact coprime failure sets (where the criterion does not fire; computed exactly in
the companion script, PART 2):

| r | criterion fails exactly for c in |
|---|---|
| 1 | (never — and every proper gcd stratum also closes, c in [2,300] checked) |
| 2 | {2} |
| 3 | {2, 3} |
| 4 | {2, 3, 4, 7, 8} |
| 5 | {2, 3, 4, 5, 7, 8, 9, 10, 14, 15} |
| 6 | {2,...,12, 14,...,18, 21, 22, 23, 24, 28, 29, 30, 35, 36, 42} |

The windows are non-monotone because the per-exception `+1` overhead interacts with
the `floor(c/7)` jumps at multiples of 7 (note r=6 already closes at c = 13, 19, 20,
25, 26, 27, 31–34, 37–41). `c >= 43` is uniform for every `r <= 6`; `c = 42` is the
last uniform failure.

**Sharp `r` threshold.** At `delta = 1/14` each exception burns an open phase arc of
length `2*delta = 1/7`; seven exceptions can burn total measure `7 * (1/7) = 1` — the
whole sheet cycle. So `r <= 6` is structurally sharp for the plain union bound: the
`r = 7` tight case is a TILING of the sheet cycle `Z_c` by rotated arc-preimages
(the unnumbered opus-S299 follow-on: the discrete 7-clock, one level down;
compare HYP-6825's independent node/fibre tiling discipline). For general LRC(n) the same
argument closes `r <= ceil(n/2) - 1` exceptions.

## Proof

**(i)** For `p in P` and any integer `k`:
`||(cp) t_k|| = ||p(t0 + k)|| = ||p t0 + pk|| = ||p t0||`, since `pk` is an integer.

**(ii)** Write `g = gcd(w, c)`, `w = g w'`, `c = g c'` with `gcd(w', c') = 1`. The
exception's phase at sheet `k` is

`w t_k = (w t0)/c + (w k)/c  (mod 1)`.

As `k` runs over `Z_c`, the residue `wk mod c = g*(w'k mod c')` runs over the
subgroup `g*Z_{c'}` with every value attained exactly `g` times (multiplication by
`w'` permutes `Z_{c'}`). So the phase multiset is the translated arithmetic grid
`{theta_0 + j/c' mod 1 : j = 0, ..., c'-1}`, `theta_0 = (w t0)/c mod 1`, each grid
point carrying multiplicity `g`.

The bad condition `||phi|| < delta` is an open arc of length `2*delta < 1`. A grid
of spacing `1/c'` meets an open arc of length `L` in at most `floor(L*c') + 1`
points: if `m` grid points lie in the open arc, the outermost two differ by
`(m-1)/c' < L` strictly, so `m - 1 < L*c'`, i.e. `m <= floor(L*c') + 1` (the
integer-`L*c'` boundary case gives `m <= L*c'` outright, which is smaller). With
`L = 2*delta` and multiplicity `g`: `#bad <= g*(floor(2*delta*c/g) + 1)`.

**(iii)** If the bad counts sum to at most `c - 1`, some `k* in Z_c` lies in no bad
set. At `t_{k*}`: every `w_a` has `||w_a t_{k*}|| >= delta` (the closed complement of
the strict bad condition — exactly the closed inequality LRC(14) needs, consistent
with the AP/GW boundary atoms), and the core sits at its exact `t0` margin by (i).
Choosing `t0` with `min_p ||p t0|| = M(P)` (attained: continuous function on the
circle) gives `M(V) >= min(M(P), delta)`. ∎

**LRC(14) corollary.** `|P| = 13 - r` distinct speeds form an instance of the settled
`(14-r)`-runner conjecture (`14 - r <= 13` for `r >= 1`), so `M(P) >= 1/(14-r) >= 1/8
> 1/14` for `r <= 6`, and `min(M(P), 1/14) = 1/14`. The threshold table is a finite
computation on `r*(floor(c/7)+1) <= c-1`; for `c = 7m + s >= 43`, `r <= 6`:
`6(m+1) <= 7m + s - 1` iff `m >= 7 - s`, which holds for all `c >= 43` (worst case
`s = 1, m = 6`), with `c = 42` (`s = 0, m = 6`) the last uniform failure. ∎

## Why this matters (the scale frontier)

- **It is the multi-exception theorem pull-card A1 asked for** (codex atlas
  2026-07-14), with the coupling guardrail honored: the union bound runs over the
  SHARED sheet index, using no independence between exceptions.
- **It extends THM-760 along both named axes**: several exceptions (`r <= 6`) and
  arbitrary gcd strata (the `sum g_a` budget); at `r = 1` EVERY gcd stratum closes at
  EVERY scale `c >= 2`, so card A2's single-exception case is closed entirely.
- **It is scale-native where bounded-q witness banks are scale-blind**: codex-S3's
  refutation family `26*{1..12} ∪ {339}` (first rational witness at `q = 27 > 25`,
  killing the uniform q<=25 finish) is exactly a `c = 26, r = 1` sheet instance —
  THM-760/761 close it analytically with no witness search. The good-period
  denominator grows with `c`; the sheet count IS `c`. On scale rays, sheets are the
  correct clock.
- **It converts HYP-6780's obstruction into a dichotomy**: `v*(cP) = c v*(P)` says
  the capped-envelope band inflates linearly in `c`; THM-761 says large `c` dies by
  sheets. The two regimes meet at an explicit constant (`c = 43` for the uniform
  coprime form): **above it sheets win; below it the band inflation is bounded by
  the same constant** — the unnumbered opus-S299 scale-uniformization route.

## Residual (named honestly)

1. **`r >= 7`** (core `|P| <= 6`): the union bound cannot fire (`7 * 1/7 = 1`); the
   tight configurations are cyclic tilings of `Z_c`. The wall is REAL: the script
   found `c = 7`, core `{1..6}`, `W = {12, 38, 72, 96, 151, 169, 188}` with ALL
   seven sheets bad at the core optimum `t0 = 1/7` — yet that family still has
   `M = 1/7 >= 1/14` (lonely by a non-sheet route). So the wall is a wall of the
   METHOD, not (so far) a wall of the conjecture. Core margin is large
   (`M(P) >= 1/8`), so there is room: sweeping `t0` over the core's whole safe set
   (not one optimum), a sharper count exploiting the arcs' EQUAL length `1/7`, or
   one round of the capped envelope on the exceptions should decide it. NOT claimed
   here.
2. **Small-`c` failures** of the per-`(r, c)` criterion (all confined to `c <= 42`):
   finite in `(r, c)` but not in `P`; these route to the band protocol whose
   inflation is bounded (`v*(cP) <= 42 v*(P)`), or to the exact per-family check.
3. **Large `sum g_a`** (deep gcd entanglement): descend `c -> c/g` recursively
   (cards A2/C10); the recursion terminates but its bookkeeping is not written here.
4. Families with **no scale structure at all** (no `c >= 2` divides `>= 7` elements):
   outside this theorem's hypotheses; they are the low-`r_P` cores where the capped
   envelope is already strong (the unnumbered opus-S299 follow-on records the
   precise complementarity claim to prove).

## Assumption challenge and tournament analysis

Vertices: witness sheets `Z_c` (following THM-760), against runner/arc/pair-ruler
alternatives — runners fail (guardrail 15: they lose the shared-`k` coupling), and
pair-sum rulers are the wrong grain (the sheet grid is a SINGLE ruler `q = c`
carrying all `r` exceptions simultaneously). Pairwise observable: sheet `k` beats
sheet `k'` if its exception-clearance vector lexicographically dominates. The free
sheet is a SOURCE of this tournament; the union bound is the statement that the
source exists once the bad in-degree budget is under `c`. What the sheet quotient
destroys — endpoint owners, pair-sum labels — is exactly what the residual cases
(r >= 7 tilings, gcd descent) need restored; that is where the blocker complex
(HYP-6785) takes over.
