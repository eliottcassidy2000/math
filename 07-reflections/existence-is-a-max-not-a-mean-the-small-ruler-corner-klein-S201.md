---
source: klein-2026-07-09-S201
status: The good-period leg of the LRC(14) covering case has a structural axis nobody had named — large ruler
  (V ≥ Q+1) vs small ruler (V ≤ Q). The extremal tight AP {0,…,12} at its own ruler V=13 (the M=1/14 instance)
  has NO good period; both proposed a-priori closures of the near-AP branch (opus-S170 smooth grid-MEAN, and
  LEM-012's V>maxE Dirichlet bound) FALSELY claim one there, by the same mechanism: at the resonant ruler V=k the
  grid j/V lands on the equidistribution NULLS of maxgap. The small-ruler corner is DENSITY-FLOOR territory. The
  recurring lesson (3rd time now): good-period EXISTENCE is a MAX statement, never a MEAN/COUNT.
tags: [lrc14, good-period, covering-case, existence-is-a-max, density-floor, resonance, meta-pattern]
---

# Existence is a MAX, not a mean — the small-ruler corner of the covering case

**klein-2026-07-09-S201.** Owner: completely close the covering case. In the same route-breaking mindset
that killed route (c) and the arc-count bound, I stress-tested the good-period leg's *own* proposed
closures on the one family everyone kept quoting as the win — the tight AP `{1..13}`. Both closures break,
and the break reveals a structural axis of the covering case that had not been named.

## The three-times lesson

The covering case reduces to: for every hard co-offset cluster `E`, either a **good ruler period** `j`
exists (`maxgap{e·j mod V} > V/7`) or the **density floor** holds (`μ_{1/7}(E) ≥ bar_k`). The good-period
leg has been attacked three ways, and the same error recurred each time — **substituting an AVERAGE for the
MAX that existence actually is:**

1. **Arc-count (MISTAKE-127, S192):** "a good period exists because `#arcs < ρ*·Vmax`." The arc *count* is
   a discrepancy error term; 7-structure fragments `Good_E` into many thin slivers, spiking `#arcs` without
   lowering the best maxgap. Vacuous on the near-AP extremal.
2. **`c < D3` / arc-count-again (MISTAKE-128, S199–S200):** `#arcs/spread < D3` — false at all spreads for
   7-structured difference sets; the *count* is fooled by the `1/7` resonance while the best maxgap is huge.
3. **Smooth grid-MEAN (this session, refuting opus-S170):** "a good period exists because the ruler-grid
   mean `E_j[maxgap] > 1/7`, and `E_j ≈ E_x[maxgap] > 1/7` up to a `0.006` discrepancy." The *mean* is fooled
   by the resonant ruler.

Each time: the averaged/counted invariant clears the bar on generic clusters, and the EXTREMAL/RESONANT
cluster silently violates the hidden gap between the average and the max. Existence is `max_j maxgap(j/V) >
1/7`. You cannot certify a max with a mean unless you control the gap — and at the exact clusters that
matter, that gap is maximal.

## The sharpest instance: the extremal AP has NO good period

`E = {0,1,…,12}` (`= {1..13}` shifted) is the extremal LRC(14) instance — `M(S) = 1/14` exactly — at ruler
`V = 13`. For **every** `j ∈ {1,…,12}`, `{i·j mod 13 : i=0..12}` is a permutation of all 13 residues, so
`maxgap = 1/13 = 0.077 < 1/7`. **No good period exists.** Yet:

- `E_x[maxgap] = 0.211 > 1/7` (the continuum mean — for a *generic* real `x` the 13 points do leave a
  `~0.2` gap). opus read this as "the tight AP closes via the mean."
- `E_grid[maxgap] = 1/13 = 0.077` (the ruler-grid mean). The grid `j/13` lands **exactly on the
  equidistribution nulls** of `maxgap(x)` — its minima. `disc = |E_grid − E_x| = 0.134`, not `≤ 0.006`.

The `α>1` Fourier-decay argument (opus's engine) fails not because `maxgap` isn't smooth, but because the
resonance `nV = 13, 26, …` hits the **head** of `maxgap`'s spectrum: for consecutive velocities `maxgap(x)`
is strongly `1/13`-periodic, so its `m=13` coefficient is dominant, not tail. The grid samples the function
in anti-phase with itself. This is the covering-case twin of the Mertens caution (S198): a smooth surrogate
tames the L¹ divergence but **not** the resonance when the sampling period *is* the resonance.

## The un-named axis: large ruler vs small ruler

The good-period leg secretly lives on `V ≥ Q+1`, `Q = ⌈7(L−1)/(L−k+6)⌉ = O(k)`:

- **`V ≥ Q+1` (large ruler):** Dirichlet's clustering `j ∈ {1,…,Q} ⊂ {1,…,V−1}` is a valid non-trivial
  period; the AP clusters, the complement is one big gap, the `≤5` strays can't fill it. LEM-012 (now with
  the corrected hypothesis) + LEM-013 close this. **MAX-based** — they exhibit the good `j`.
- **`V ≤ Q` (small ruler):** the AP wraps and equidistributes on the coarse grid; Dirichlet's `j` collapses
  to `j ≡ 0 (mod V)`, the excluded trivial period; often **no good period exists**. This is the
  extremal/dense corner — and it is **DENSITY-FLOOR territory** (`μ_good({0..12}) = 0.44 ≥ bar_13`).

So the covering-case dichotomy is really a `2×2`: {near-AP, dissociated} × {large ruler, small ruler}. The
good-period leg owns three cells; the fourth (near-AP ∩ small-ruler = the extremal tight AP) is the density
floor's. Both prior closures tried to make the good-period leg own all four, and the fourth cell is exactly
where they had no right to.

## Why this matters for "completely closing"

Nothing is lost — the covering case still closes (good-period `V ≥ Q+1` ∪ density-floor `V ≤ Q`). But the
formalization target changes: do **not** try to prove `exists_good_of_smooth_mean`'s hypotheses for the tight
AP (they are false); do **not** state LEM-012 with `V > max E` (Lean would let you, and the theorem would be
wrong). The honest capstone routes the small-ruler resonant corner to the density floor **explicitly**. The
math keeps insisting: to prove something exists, show the witness (the best `j`) or change legs — never
average.
