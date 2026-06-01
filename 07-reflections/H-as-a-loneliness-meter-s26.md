---
source: oracle-2026-06-01-S26
status: exploratory result (what H measures on the circular tournament)
tags:
  - tournament-analysis
  - tournament-clock
  - lonely-runner
  - hamiltonian-paths
  - three-cycles
  - loneliness-meter
---

# H as a Loneliness Meter: what the Hamiltonian-path count actually measures

S24 noticed that for the half-turn tournament `T(t)` of a runner system, the
directed Hamiltonian-path count `H(T(t))` falls when the runners bunch and rises
when they spread — a "loneliness meter." This session pins down exactly what that
meter is, how sharp it is, and where it stops working. (20000-sample exact study,
`H_loneliness_meter_s26.py`, n = 5,6,7.)

Recall the lift: `n` points on the circle, `i → j iff (x_i − x_j) mod 1 ∈ (0,½)`.
This is a **round tournament** (each out-set is the arc of points in the leading
semicircle). The realizable iso-classes are the fixed circular menu of S24.

## 1. One sharp threshold: H = 1 ⟺ an empty semicircle (½-gap)

Across all n, `H = 1` (transitive) holds **exactly** when `max_gap ≥ ½`:

```
n=5: H=1 occurs for max_gap ∈ [0.5000, 0.962];  H>1 only for max_gap < 0.5
n=6: H=1 for max_gap ∈ [0.5000, 0.914]
n=7: H=1 for max_gap ∈ [0.5000, 0.908]
```

So the meter has **one perfectly sharp reading**: `H = 1 ⇔ all runners fit in an
open semicircle ⇔ a gap of length ≥ ½` (the `n = 2` lonely-runner threshold). This
is the cleanest statement of "H detects loneliness": it detects the ½-gap, sharply.

## 2. Above the threshold, H is only a *coarse* meter

For `H > 1` the meter blurs. The `max_gap` ranges of the different `H`-values
**overlap heavily**:

```
n=5:  H=9  max_gap∈[0.259,0.500]   H=11 ∈[0.252,0.499]   H=15 ∈[0.219,0.485]
```

and `max_gap` (to 2 dp) maps to several distinct `H` (26/74 values at n=5). So
**H is not a function of `max_gap`** and the `H`-ranges are **not strictly
ordered** in `max_gap`. What *is* monotone is the **mean**: as `H` rises the mean
`max_gap` falls monotonically (n=7: 0.570, 0.431, 0.376, 0.353, 0.319, 0.296,
0.336*, 0.294, 0.262). H ranks bunching *on average*, not pointwise — because the
tournament is a **quantization** of the gap profile (it records only which cumulative
arcs cross ½), so many gap profiles share one tournament and one `H`.

**Resolution limit.** A single half-turn tournament can sharply detect *only* the
½ threshold; all finer gap information is quantized away. This is exactly why the
tournament-clock lens lives at threshold ½, and why the lonely runner's finer
thresholds `1/n` need richer comparator families (S25, the `α = 1/k` lift). The
loneliness a single tournament can *certify* is precisely the half-circle kind.

## 3. What the meter counts: spread-triples (#3-cycles), and H is finer

The cleanest meter is the **number of 3-cycles**, with an exact geometric meaning:

> **#3-cycles(T) = #{ triples of runners whose three pairwise arcs are all < ½ }
> = #{ triples not contained in any semicircle }.**

(Three circle points form a cyclic triangle under the half-turn rule iff none of
the three arcs between them reaches ½, i.e. they "surround the centre.") So the
3-cycle count is literally a **spread-triple counter** — a local loneliness meter
at the triple level — and it equals the classical `C(n,3) − Σ_i C(s_i,2)`
(`s_i` = scores). At n=5 it is a perfect monotone bijection with H
(`#3cyc = 0,3,4,5 ↔ H = 1,9,11,15`).

But **H is strictly finer than the score-based meters.** At n=6 the classes
`H = 41` and `H = 45` share score sequence `(2,2,2,3,3,3)` and `#3cyc = 8` and
score-variance `0.25` — yet differ in `H`. Same at n=7 (`H = 123` vs `137`, both
`#3cyc = 12`). Hamiltonian paths see the *cyclic arrangement*, not just the degree
profile, so H separates configurations that the spread-triple count cannot.

The meter hierarchy, coarse → fine:

```
min out-degree = 0   ⇔  ½-gap            (binary: is there an empty semicircle?)
score sequence / score-variance / #3-cycles  (the spread profile; equivalent)
H = Hamiltonian-path count                    (strictly finer; sees arrangement)
```

All are monotone in spread; only the binary ½-gap reading is sharp.

## 4. Tie to the lonely runner

H meters loneliness at the `n = 2` resolution (½). For the lonely runner at level
`n` (threshold `1/n`) the relevant event is a `1/n`-gap around a *specific*
runner — sub-½ and *local*, which a single half-turn tournament cannot resolve.
Two consequences:

- The right per-runner statistic is the **out-degree** `s_i` = #runners in runner
  `i`'s leading semicircle; `s_i = 0` flags that runner `i` has an empty semicircle
  ahead. The score *sequence* is the per-runner ½-loneliness profile, and `#3-cycles`
  is its symmetric summary.
- To reach the `1/n` threshold one must either refine the comparator (S25's `α`
  family) or read the **whole clock walk**, whose cell structure encodes the
  sub-½ gaps the static tournament hides. H is the meter at the coarsest rung;
  the LRC threshold lives several rungs finer.

So "H is a loneliness meter" is true and now precise: **H sharply certifies the
half-circle gap, monotonically (on average) ranks finer spread, counts — through
its sibling #3-cycles — the triples that fail to bunch, and resolves arrangement
detail the score profile misses; but it is intrinsically a ½-resolution
instrument, which is exactly the boundary where Tournament Analysis meets the
lonely-runner thresholds.**

## Next

1. Prove #3-cycles = non-semicircle-triple count for the half-turn tournament
   (immediate from the three-arc characterization) and `= C(n,3) − Σ C(s_i,2)`.
2. Characterise the H=41/45-type pairs: which finer feature of the cyclic gap
   pattern does H separate? (candidate: the second-longest arc / the reach word.)
3. Build the `α = 1/k` meter and show its "transitive" reading is the `1/k`-gap,
   recovering H as `k = 2` and giving a genuine LRC-threshold loneliness meter.

## Artifacts

```
04-computation/H_loneliness_meter_s26.py
05-knowledge/results/H_loneliness_meter_s26.out
```
