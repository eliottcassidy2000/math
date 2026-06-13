# Tiling Sequences and Cut Geometry

**Session:** opus-2026-05-27-S2

## The Good-Cuts Lattice

Every tiling τ has a good-cut set G(τ) ⊆ {1,...,n-1} (cuts with ≥1 upward tile). The distribution of |G(τ)| across all 2^m tilings follows a remarkable pattern:

|G|=0: 1 tiling (transitive)
|G|=1: IMPOSSIBLE
|G|=2: n-2 tilings
|G|=3: 5(n-3) tilings
...
|G|=n-1: SC(n) tilings

The impossibility of |G|=1 is the cleanest result: every tile spans ≥2 consecutive cuts (because tiles are non-consecutive pairs, so x-y ≥ 2 cuts). This means any upward tile "opens" at least 2 cuts simultaneously. A tournament can be "barely non-transitive" (open exactly 2 consecutive cuts) but cannot "barely non-transitive" via a single cut.

This is the analog of: **you can't have a half-bridge**. A bridge either connects two components completely (cuts all the way through) or not at all (doesn't help). Here, every upward tile "bridges" across at least 2 cuts.

## Sub-SC Counts as Coefficients

The coefficient 5 in "5(n-3) tilings with exactly 3 good cuts" is exactly SC(4) — the number of strongly connected tilings on 4 vertices. Similarly, 50 appears in the d=4 formula as SC(5). The pattern:

exactly-d-good ≈ (n-d) × SC(d+1) + (correction)

interprets beautifully: to have exactly d consecutive good cuts {k,...,k+d-1}, the free tiles within the interval must form an SC tiling of the sub-problem (a (d+1)-vertex tiling model). The "correction" comes from non-consecutive good sets, which appear first at d=4.

The sub-tiling counts SC(d+1) = 1, 5, 50, 903, 30773, ... (our new SC tiling sequence!) are the building blocks for the exact d-good-cut formulas.

## The Non-SC Asymptotics: P(non-SC) → 8/2^n

The non-SC tiling probability approaches 8/2^n = 2^{3-n}:
- n=7: 1995/32768 ≈ 1/16.4 vs 8/128 = 1/16
- n=11: 137B/35T ≈ 1/256 vs 8/2048 = 1/256

This means: in a "random" tiling, the probability of being non-SC halves every time n increases by 1. Or equivalently: the probability of SC approaches 1, with the deficit shrinking as 8/2^n.

The constant 8 = 2^3 comes from the two boundary cuts (k=1 and k=n-1, which have the smallest f-values = n-2) each contributing 2^{m-(n-2)} = 2^{m-n+2} to the non-SC count. Together: 2·2^{m-n+2} = 2^{m-n+3} = 8·2^{m-n} = 8·(2^m / 2^n).

## The f(S) Formula: Intersection Extremes

The most beautiful result: f(S) depends on intersections |tiles_k₁ ∩ tiles_k₂ ∩ ...| which reduce to min(S)·(n-max(S)). Only the EXTREME cuts in S matter for deep intersections! The intermediate cuts in S are "invisible" in the intersection formula.

This resonates with the staircase geometry: the "reach" of cut k (how many tiles cross it) = k(n-k)-1, which is symmetric, peaking at k = n/2. The deepest cut (most constraining) is the central one. But for intersections of cuts, only the top and bottom of the staircase matter — the minimum and maximum of S.

This is the staircase analog of the "extreme value principle": for crossing probability in a random walk, only the range (max-min) matters, not the intermediate positions.

## Two New OEIS Sequences

The SC tiling sequence 1, 5, 50, 903, 30773, 2032504, 264271477, ... and the non-SC sequence 1, 3, 14, 121, 1995, 64648, 4163979, ... are both new (not in OEIS as of 2026-05-28). They are computationally accessible via the IE formula (fast, no brute force needed) and would be natural additions to OEIS's tournament section.

The SC sequence in particular is the "strongly connected base-path-fixed tournaments" — a natural variant of A054946 (labeled SC = all base paths free) and A051337 (unlabeled SC).

## The H=7 Gap at n=5

H=7 is impossible at n=5 because the OCF constraint 1+2α₁+4α₂=7 ⟹ α₁+2α₂=3 has no realizable solution in a 5-vertex tournament. This is not obvious from tournament axioms alone — it requires the cycle structure to be globally consistent. The OCF turns a mysterious combinatorial gap into an arithmetic impossibility.

This suggests looking for similar gaps at n=6,7,8: what H values are impossible? At n=6 the achievable set seems much denser, but there may still be gaps. Worth computing.
