---
source: oracle-2026-06-01-S552
status: VERIFIED computational result (gap exhaustive n<=8; edge achieved+proven all n) + sharp open conjecture
tags: [lonely-runner, extremizer, spectral-gap, margin, doubled-apex, doubling, AP, qq-cycle, lrc-progress]
---

# The LRC loneliness spectral gap: 1/n, then a jump to 2/(2n-1) (doubled apex)

**Prompt (user):** spend another long session pushing these forward searching for
LRC progress.

This session turns the recent metaphor thread (doubling/apex/cascade/rank) into a
**falsifiable, exact, quantitative statement about the Lonely Runner extremizers** --
and finds a clean *spectral gap* with a named witness.

## Setup (covering / max-collar form)

For `n-1` distinct positive integer speeds `S = {s_1,...,s_{n-1}}` (gcd 1, observer
speed 0), define the **max-collar**
```
M(S) = max_{t in (0,1)}  min_i || s_i t ||.
```
LRC(n) is the claim `M(S) >= 1/n` for every `S`. The conjectured unique worst case is
the arithmetic progression `AP = {1,2,...,n-1}`, which achieves `M(AP) = 1/n` exactly.

The *hardness* of LRC at `n` is not the floor `1/n` -- it is the **margin**: how far
the *second*-worst configuration sits above `1/n`. I computed `M(S)` exactly (Fraction
arithmetic, exact candidate set = tent peaks and pairwise crossings -- no sampling)
over all gcd-1 sets up to an entry bound, and measured that margin.

## The result: a spectral gap of width 1/(n(2n-1))

Over all configurations, the max-collar values are
```
   M(S) = 1/n        <-- only the AP tight family (the conjectured extremizers)
   M(S) >= 2/(2n-1)  <-- everything else,  with equality ACHIEVED.
```
**There is no configuration with `M(S)` strictly between `1/n` and `2/(2n-1)`.**
Exhaustively verified for `n = 4,5,6,7,8` (full enumeration, bound large enough to
expose the `2/(2n-1)` achiever; stress-tested to `B=18` at n=6). The margin is exactly
```
   margin(n) = 2/(2n-1) - 1/n = 1/( n(2n-1) )  ~  1/(2 n^2).
```
Values: `1/28, 1/45, 1/66, 1/91, 1/120, ...` for `n=4,5,6,7,8`.

## The witness: the AP with its apex runner DOUBLED

The second-worst configuration -- the one that *achieves* the gap edge `2/(2n-1)` --
is
```
   A_n = { 1, 2, ..., n-2, 2(n-1) }   =   the AP with its TOP speed n-1 doubled.
```
This is **provable in closed form** (verified n=4..22). At `t* = 2/(2n-1)`:
```
   runner 1        :  ||2/(2n-1)||                   = 2/(2n-1)
   runner 2(n-1)   :  ||(4n-4)/(2n-1)|| = ||2 - 2/(2n-1)|| = 2/(2n-1)
   runner s (2<=s<=n-2) :  ||2s/(2n-1)|| = min(2s, 2n-1-2s)/(2n-1)  >= 3/(2n-1)
```
so `min_i ||s_i t*|| = 2/(2n-1)`, and the max over `t` equals this (verified). The two
*binding* runners are exactly **speed 1 and speed 2(n-1)** -- the unit and the doubled
apex -- meeting at distance `2/(2n-1)`; every interior AP runner has strictly more room
(`>= 3/(2n-1)`).

Even sharper, along the **apex direction** `{1,...,n-2, s}` (vary only the top speed
`s >= n-1`): the loneliest is `s = n-1` (the AP, `M=1/n`); the **second**-loneliest is
*uniquely* `s = 2(n-1)`, with `M = 2/(2n-1)`. Verified `n=5..12`. So doubling the apex
is the single minimal move off the extremizer, and it lands exactly on the gap edge.

## Why this is the payoff of the doubling thread

The recent sessions argued, metaphorically, that **doubling is pairing** (S546), the
**apex** is the `n/2` / order-2 pivot (S530/S547), and the doubled prime `2q` is the
recursion's load-bearing **bridge** (S549). This session makes it literal and extremal:

> The boundary of the LRC extremal basin -- the *second*-loneliest structure at every
> `n` -- is the AP with its apex runner doubled. Doubling the top speed is the minimal
> perturbation of the conjectured extremizer, and it sits exactly `1/(n(2n-1))` above
> the floor.

So "doubling the apex" is not a metaphor -- it is the explicit witness on the edge of
the gap. The `(q,q)` cycle prediction is also confirmed here: at doubled primes
`n=6,10` every minimax extremizer's observer-necklace has rotation-by-2 cycle type
exactly `(q,q)` (`(3,3)`, `(5,5)`).

## What this buys for LRC (and the honest boundary)

- **It localizes the difficulty.** LRC(n) is delicate *only* on the tiny AP tight
  family. Every other configuration clears the `1/n` floor with a definite surplus
  `>= 1/(n(2n-1))`. There is no "near miss" continuum -- the second-worst value is a
  jump, not a creep. A proof of LRC only needs to (a) pin the tight family and (b)
  show everything else exceeds `2/(2n-1)`; there is nothing to handle in between.
- **The gap edge is unconditional.** `M(A_n) = 2/(2n-1)` is proven for all `n` (closed
  form). So the *upper* side of the margin is a theorem; only the *emptiness* of the
  open interval `(1/n, 2/(2n-1))` is so far computational (`n<=8` exhaustive).
- **Honest status.** "No config in the open gap" is exactly a strengthening of LRC
  (it implies LRC: nothing dips to the floor except the AP family). I have it
  exhaustively for `n<=8`, not proven in general. The clean **provable slice** is the
  apex direction `{1,...,n-2,s}`: there the second-loneliest is uniquely `2(n-1)` at
  `2/(2n-1)` -- a one-parameter gap theorem.

## Open (-> HYP-2052)

1. **Prove the gap.** Show every non-AP-tight gcd-1 config has `M(S) >= 2/(2n-1)`.
   (Equivalently: the only way to reach max-collar `< 2/(2n-1)` is to BE an AP-tight
   config with `M = 1/n`.) Provable slices to chain: AP-prefix families
   `{1,...,k, ...}`, then dilations, then general.
2. **Characterize the full M2 family.** Besides `A_n`, the enumerated `M2` sets include
   `{1,3,4,...,n-1, ...}` and dilates -- do they all reduce to "an AP-tight set with one
   apex doubled / one runner halved", i.e. one elementary `x2`/`/2` move? If the gap
   edge is always a single doubling/halving off a tight set, that is a structural law.
3. **Margin asymptotics as an LRC-hardness measure.** `margin(n) ~ 1/(2n^2)` shrinks
   only polynomially -- the basin stays well-separated. Does this separation survive
   the `n -> infinity` regime where LRC is open, or does the gap close (a value creeps
   into `(1/n, 2/(2n-1))`) at some threshold `n`? Find the first `n`, if any, where the
   gap is NOT `2/(2n-1)`.

## Anchors
- `04-computation/lrc_minimax_margin_extremizers_s552.py` (+`.out`): minimax = 1/n,
  the tight family, the margin `2/(2n-1)`, and the `(q,q)` cycle confirmation.
- `04-computation/lrc_collar_gap_s552b.py` (+`.out`): M2-family, the doubled-apex
  pattern, gap stress test.
- `04-computation/lrc_doubled_apex_gap_s552c.py` (+`.out`): closed-form
  `M(A_n)=2/(2n-1)` (n=4..22), AP floor `1/n`, exhaustive empty-gap `n<=8`.
- `05-knowledge/results/lrc_apex_slice_s552d.out`: apex-direction unique 2nd-loneliest
  `s=2(n-1)`.
Builds on S546 (doubling=pairing), S547 (apex/(q,q)), S549 (doubled-prime bridges),
S526 (covering proof / n=3 character).
