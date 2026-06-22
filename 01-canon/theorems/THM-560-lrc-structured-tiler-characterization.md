---
id: THM-560
title: LRC(14) structured-tiler characterization -- the difference-closed exact tilers are exactly the dilated intervals d*{1..13}, all tight
status: PROVED (elementary)
author: kind-pasteur-2026-06-22-S31n
depends_on:
  - HYP-+2888   # mac-mini: exact-coverage extremal is scaling-invariant (this resolves its structured case)
  - OPEN-Q-108  # tight-locus finiteness (this isolates the residual SPORADIC case)
verified: computational (lrc_tiling_rigidity_kps.py, lrc_tight_vs_counterexample_kps.py) + proof below
---

# THM-560: the difference-closed exact tilers of LRC(14) are exactly d*{1..13}

## Setup
For a set `S` of 13 distinct positive integers (the runner speeds of LRC(14), relative frame),
define the **lonely constant** `M(S) = max_{t in [0,1)} min_{s in S} ||s t||` and the **safe set**
`Safe(S) = { t in [0,1) : ||s t|| >= 1/14 for all s in S }`. LRC(14) for `S` <=> `Safe(S)` nonempty.
Call `S` an **exact tiler** when `meas(Safe(S)) = 0` (the danger arcs `U_s={||st||<1/14}` cover
`[0,1)` up to measure zero). Call `S` **difference-closed** when `|s - s'| in S` for every
`s != s'` in `S`.

mac-mini (HYP-+2888) reframed the LRC(14) wide-bound crux as: *which `S` are exact tilers?* and
conjectured "only `d*{1..13}`." THM-560 resolves this **for the difference-closed sets** (and shows
the conjecture is INCOMPLETE -- a sporadic tiler exists).

## Statement
**(a)** A 13-element set `S` of positive integers is difference-closed **iff** `S = d*{1,...,13}`
for some integer `d >= 1` (a dilated interval).

**(b)** Every dilated interval `d*{1,...,13}` is an exact tiler, and is **tight**:
`M(d*{1..13}) = 1/14`, with `Safe = { j/(14 d) : gcd(j,14)=1 }` -- exactly 6 points, measure 0.

**(c)** Hence among difference-closed 13-sets, *exactly* the AP-dilates tile, all at `M = 1/14`.
The remaining exact tilers are **non-difference-closed (sporadic)**; one exists -- the
Goddyn-Wong-type set `S_GW = {1,...,11, 13, 24}` (verified `M = 1/14` exactly, witness `t = 5/14`,
`meas(Safe)=0`). Their classification is the residual hard core (OPEN-Q-108).

## Proof
**(a, <=)** `d*{1..13}` has nonzero differences `d*{1,...,12} ⊆ d*{1,...,13}`. Difference-closed.

**(a, =>)** Let `S` be difference-closed, `d = min(S)`. For `s in S` with `s > d`, the difference
`s - d > 0` lies in `S`. Iterating `s -> s-d -> s-2d -> ...` stays in `S` and `>= d` (minimality),
so the chain terminates at `d`; thus `s = (k+1)d` is a multiple of `d`. So `S = d*S'` with `S'`
difference-closed and `min(S') = 1`. For `s' in S'`, `s'-1 in S'` (difference with `1`), so
`s', s'-1, ..., 1` all lie in `S'`; hence `S' ⊇ {1,...,max(S')}`. With `|S'| = 13`, `max(S') = 13`
and `S' = {1,...,13}`. Therefore `S = d*{1,...,13}`. ∎

**(b)** Fix `S = d*{1..13}` and `t in Safe(S)`. Put `u = (d t) mod 1`. Consider the 14 points
`P = { k u mod 1 : k = 0,1,...,13 }`. Any two have circle distance
`|| (i-j) u || = || |i-j| u ||` with `|i-j| in {1,...,13}`. Now `||m u|| = ||m d t|| = ||(md) t||`
and `md in d*{1..13} = S` for `m in {1..13}`, so by `t in Safe(S)` every pairwise distance is
`>= 1/14`. Fourteen points on a circumference-1 circle, pairwise `>= 1/14`: the 14 cyclic gaps sum
to `1` and each (a distance between adjacent points) is `>= 1/14`, forcing **all gaps `= 1/14`**.
So `P` is the equally-spaced set `{0, 1/14, ..., 13/14}`, which forces `u = j/14` with
`gcd(j,14)=1` (else `{ku}` misses residues). Hence `d t ≡ j/14`, i.e. `t in {j/(14d) : ...}` --
finitely many points, so `meas(Safe) = 0`. At such `t`, `min_s ||st|| = min_k ||kj/14|| = 1/14`
(attained at the `k` with `kj ≡ ±1 mod 14`), and no `t` gives `min > 1/14` (that would be an open
safe set, contradicting `meas(Safe)=0`); therefore `M = 1/14`. ∎

**(c)** Immediate from (a),(b); `S_GW` checked computationally (not difference-closed: e.g. `24-1=23
∉ S_GW`; tight with `M = 1/14`). ∎

## Significance
- **Resolves the structured half of mac-mini's crux**: the difference-closed (equivalently
  "interval-like") exact tilers are *completely* pinned down -- exactly the AP-dilates, all tight,
  by a one-paragraph pairwise-distance/equal-spacing argument. No additive energy, no convexity, no
  measure-LP -- just difference-closure forcing equal spacing.
- **Corrects HYP-+2888**: "only `d*{1..13}` tile" is FALSE -- `S_GW` tiles too. The tight locus is
  at least `{ AP-dilates } ∪ { GW-type }`, matching the known LRC tight locus `{AP, Goddyn-Wong}`.
- **Isolates the genuine hard core**: OPEN-Q-108 = "classify the SPORADIC (non-difference-closed)
  exact tilers." The structured case is no longer the obstruction; the sporadic finiteness is.
- Connects to the project frame: difference-closure = the speed set being its own difference set =
  the AP being the unique *self-similar* (Sidon-opposite) speed structure; the equal-spacing rigidity
  is the 1-D shadow of the regular/transitive extremality that runs through the whole tournament story.

-> OPEN-Q-108, HYP-+2888 (mac-mini), HYP-2885, `lrc-coverage-transcends-the-h-level-...md`.
