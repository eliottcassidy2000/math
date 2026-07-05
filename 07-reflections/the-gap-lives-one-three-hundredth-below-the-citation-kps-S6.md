# The gap lives one three-hundredth below the citation

**kind-pasteur-2026-07-05-S6 (HYP-4112)**

## The observation

The spectral-gap crux of LRC(14) asks: can a primitive 12-set have
`M(W) ∈ (1/13, 2/25)`?  The window's ceiling is `2/25 = 0.08`.  The LRC(≤13)
citation, applied one peel level down, gives every **11**-runner family a time
with margin `1/12 = 0.0833…`.

`1/12 > 2/25`.  By exactly `1/300`.

That strict inequality is a proof engine. In a gap violator, no time has margin
`2/25` for all twelve runners. But the citation hands each 11-runner complement
a `1/12`-margin witness `t'`, and Lipschitz keeps the complement's margin above
`2/25` on a window of half-width `(1/300)/B` around `t'`. So the excluded runner
must cover that entire window with its own `2/25`-teeth — and a comb of speed
`w` has teeth only `(4/25)/w` wide. One-tooth containment forces

    w ≤ (2/25)·B / (1/300) = 24·B.

**Every runner of a gap violator is at most 24× the largest of the others**
(`peel_height_bound`, now kernel-pure). In particular `w_max ≤ 24·w_2nd`: gap
violators are 24-compressed (`gap_compressed_24`). No 11-level rigidity, no
tower recursion, no census — just the citation, a Lipschitz window, and a comb
that is too sparse to hide a whole interval.

## Why this transcends the session

1. **The gap is trapped between two floors of the same conjecture.** The
   spectral gap `(1/13, 2/25)` sits below the 11-runner citation floor `1/12`.
   The conjecture's own settled cases press down on its first open case from
   one level up — the tower is not just a proof strategy, it is a vice. The
   margin of pressure, `1/300`, is the difference between "the settled cases
   help" and "they don't": had the 12-runner second value been `1/12` or more
   (instead of `2/25`), the peel would give nothing. The numerology `25 > 24`
   (i.e. `2·12 < 25`) is load-bearing.

2. **Compression is contagious downward.** The same peel at deeper levels
   (remove `l` runners; citation floor `1/(13−l)`; the `l` removed runners must
   jointly cover the window, fee-ledger style) gives a whole ladder of
   constraints — the `l = 1` rung is the sharpest and is now formal. klein's
   multi-killer fee machinery (S134/S135) is exactly the tool for the deeper
   rungs; run at target `2/25` instead of `1/14`, it would bound every order
   statistic against the next. The gap violator's shape is squeezed from the
   top (compression), the bottom (spread, ratio > 11.5), and the middle
   (big binding pair ≥ 38): a creature that must be everywhere large and
   nowhere dominant.

3. **The crux is now a bounded-geometry problem.** Before today the gap
   violator's height structure was unconstrained. Now: `w_max ≤ 24·w_2nd`,
   two elements sum ≥ 38, ratio > 11.5 overall, all moduli ≤ 12 covered
   inside the set, all moduli ≤ 25 pinned. What is NOT yet bounded is the
   overall scale (dilation-like freedom for non-AP shapes remains). The
   remaining question has the flavor: can a 24-compressed, spread, everywhere-
   pinned 12-set cover the circle at radius 2/25 more efficiently than its own
   arithmetic allows? That is a finite-geometry question per scale — the census
   machinery's home turf.

Related: [[the-gap-violator-wears-handcuffs-at-every-prime-kps-S3]],
[[the-loose-witness-is-one-decidable-inequality-kps-S2]].
