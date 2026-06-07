---
id: THM-427
title: The signed pairwise LRC floor is a max-cut LRC — lower bounds 1/(2·r_min) (unconditional) and 1/(r_min+1) (LRC); positive for every fixed n, but below 1/n for n ≥ 6
status: PROVED (lower bounds elementary + LRC-conditional; computationally verified n≤7)
source: monad-explorer-2026-06-06-S3
depends_on:
  - THM-426    # sign patterns = cuts of K_{n-1}; Gstar(S)=max over cuts of M(relative-speed set)
related:
  - HYP-2293   # REFUTED: Gstar can be < 1/n (S2 example 3/19 — now seen to be a small-B artifact)
  - HYP-2294   # forward: the n-asymptotic of inf_S Gstar(S) (this theorem brackets it)
  - THM-420    # opus-S700 k-clock/shell-partner: the SAME LRC loneliness bound, one level down
  - THM-369    # n∤v ⇒ M_obs ≥ 1/n (the observer floor this theorem is compared against)
---

# THM-427 — the signed pairwise floor is a max-cut LRC

**Convention** (repo canon): `n` runners total, observer speed `0`, movers
`S = {v_1 < ⋯ < v_r}` distinct positive integers, `r = n−1`, `gcd(S)=1`, `‖x‖` = distance to
nearest integer.

By **THM-426** the signed pairwise gap depends on the sign pattern `ε∈{±1}^r` only through the
bipartition `(A,B)` of the movers, and
```
   Gstar(S) = max over bipartitions (A,B) of  M(W(A,B)),
   W(A,B)   = { |v_i−v_j| : i,j same side } ∪ { v_i+v_j : i,j across the cut },
   M(W)     = max_t min_{w∈W} ‖w t‖     (the LRC loneliness of the relative-speed set;
              only the DISTINCT values of W matter).
```
So `Gstar` is a **max over cuts of a Lonely-Runner loneliness** — an LRC quantity wrapped in a cut
optimisation. This theorem bounds it below by LRC itself.

Let `k(W)` = number of distinct (nonzero) values in `W`, and
```
   r_min(S) = min over bipartitions of k(W(A,B))     (the cut-minimised distinct-rel-speed count).
```
All elements of every `W` are positive (`v_i+v_j>0`; `|v_i−v_j|>0` since speeds distinct), and
`r_min(S) ≤ C(r,2) = C(n−1,2)` (any single cut already has at most `C(r,2)` pairs).

## Two loneliness lemmas

> **Lemma A (measure / union bound — UNCONDITIONAL).** For any finite set `W` of `k` distinct
> positive integers, `M(W) ≥ 1/(2k)`.
>
> *Proof.* For each speed `w`, `{t∈[0,1): ‖wt‖<δ}` is a union of `w` intervals of length `2δ/w`,
> total measure `2δ` (for `0<δ≤½`). The union over the `k` speeds has measure `< 2δk`. If
> `2δk < 1` the union is a proper subset of `[0,1)`, so there is `t*` with `‖w t*‖ ≥ δ` for every
> `w`, giving `min_w ‖w t*‖ ≥ δ`. Hence `M(W) = sup_t min_w ‖wt‖ ≥ δ` for every `δ < 1/(2k)`, so
> `M(W) ≥ 1/(2k)`. ∎

> **Lemma B (LRC bound — CONDITIONAL).** If the Lonely Runner Conjecture holds for `k` runners
> (`M(W) ≥ 1/(k+1)` for any `k` distinct positive integer speeds), then `M(W) ≥ 1/(k+1)`. This is
> **unconditional for `k ≤ 6`** (Barajas–Serra 2008, the 7-runner case). It is open for `k ≥ 7`.

## The theorem

> **THM-427.** For every `S`:
> 1. **(unconditional)** `Gstar(S) ≥ 1/(2·r_min(S)) ≥ 1/(2·C(n−1,2)) = 1/((n−1)(n−2))`.
> 2. **(LRC)** `Gstar(S) ≥ 1/(r_min(S)+1)`; unconditional when `r_min(S) ≤ 6`, conditional on LRC
>    otherwise.
>
> *Proof.* Pick the bipartition achieving `k(W*) = r_min`. Then `Gstar(S) ≥ M(W*)`; apply Lemma A
> (resp. Lemma B) with `k = r_min`. The absolute bound uses `r_min ≤ C(n−1,2)`. ∎

> **Corollary 1 (fixed-`n` positivity).** For each fixed `n`, `inf_S Gstar(S) ≥ 1/((n−1)(n−2)) > 0`.
> The signed pairwise floor is **bounded away from 0 for every `n`**, even though (HYP-2293) it
> dips below the *observer* floor `1/n` for `n ≥ 6`.

> **Corollary 2 (when the observer floor survives).** If `r_min(S) ≤ n−1`, then `Gstar(S) ≥ 1/n`
> (unconditional when `n ≤ 7`, since then `r_min ≤ 6`; conditional on LRC otherwise). Contrapositive:
> `Gstar(S) < 1/n ⟹ r_min(S) ≥ n` — **every cut must leave ≥ n distinct relative speeds.**

So the breakdown of the observer floor in the pairwise problem is exactly an obstruction to making
the relative-speed set small: a config drops below `1/n` only when **no** 2-colouring of the movers
keeps the distinct-relative-speed count `≤ n−1`. The mechanism is an irreducible **cluster of small
speeds** (e.g. three near-consecutive movers): a cut can convert a small *difference* into a large
*sum*, but it cannot do so for all pairs of the cluster at once, so a small relative speed and extra
distinct values survive in every colouring (see the reflection).

## Computational verification (exact arithmetic)

`04-computation/signed_lrc_rmin_bound_monad_s3.py` (+ `.out`):
- Lemma A bound (i), Lemma B bound (ii, for `r_min ≤ 7`), and Corollary 2 hold with **0 failures**
  over all `gcd=1` sets at `n = 3..8` (within the speed bounds tested).
- **Every** below-`1/n` set has `r_min ≥ n` (Corollary 2 contrapositive verified, 0 exceptions).
- The observer floor `Gstar(S) ≥ 1/n` holds for **all** `gcd=1` `S` at `n ≤ 5` (n=5 robust to
  `B ≤ 24`, 10061 sets, `inf = 1/5` exactly); the first violation is at `n = 6`, and the floor also
  breaks at `n = 7, 8` (consistent with HYP-2293 / THM-426).

`04-computation/signed_lrc_inf_highB_monad_s3.py` shows the inf is approached **from above** as the
speed bound `B` grows:
```
   n=6:  inf = 2/13, 3/20, 4/27, 6/41   for B = 10,13,16,18   (n·inf = .923,.900,.889,.878)
   n=7:  inf = 2/15, 1/8,  5/42         for B = 9,11,13       (n·inf = .933,.875,.833)
   n=8:  inf = 1/9,  2/19               for B = 9,11          (n·inf = .889,.842)
```
So the S2 value `3/19` was a small-`B` artifact. The true `inf_S Gstar(S)` for fixed `n` is positive
(Corollary 1), `n·inf_S Gstar` is `< 1` and decreasing for `n ≥ 6`, and its `n→∞` behaviour (does it
stay `Θ(1)` or decay toward the `1/((n−1)(n−2))` regime?) is the open question HYP-2294.
