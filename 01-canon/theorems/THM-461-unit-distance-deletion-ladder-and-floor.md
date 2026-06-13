---
id: THM-461
name: the-unit-distance-upper-bound-ladder-is-elementary-and-its-floor-is-exhausted
status: PROVED (elementary, complete) + VERIFIED exact-integer
date: 2026-06-10
session: monad-explorer-2026-06-10
depends_on:
  - THM-431   # u(21)=57 (AMP); N* in [25,28]; HYP-2298
  - THM-432   # avgdeg additive under product; n=27 tie = Hamming cube H(3,3)
  - THM-437   # cube angle-rigidity at 81
  - HYP-2298  # u(21)=57 proven, N* in [25,28]
sharpens:
  - OPEN-Q-057  # exact value of N*; reduces the floor side to small exact-value targets
---

# THM-461: the unit-distance upper-bound ladder above n=22 is elementary, and the 3N-floor it can reach is exhausted at N*≥25

## Context

`u(n)` = max number of **unit distances** among `n` points in the plane (Erdős).
`N* = ` smallest `n` with `u(n) > 3n` (average degree exceeds the planar kissing
number `κ=6`). Canon: `N* ∈ {25,26,27,28}` (THM-431), the sharp target being
`u(27)`: is it `81 = 3·27` (tie ⟹ `N*=28`) or `> 81`? The proven anchors are
`u(n)` exact for `n≤21` (Alexeev–Mixon–Parshall, arXiv:2412.11914; `u(21)=57`) and
the one nontrivial bound `u(22) ≤ 61`; the Moser-lattice constructions give the
lower bounds `L(n) = 60,64,68,72,76,81,85` for `n = 22…28` (Engel/Schade, repo
THM-431-C). This theorem dissects the **upper-bound (floor) side** of `N*`.

## Statement

**(A) The deletion ladder [PROVED, elementary].** For every finite point set on
`n` points,
```
        Σ_v u(G−v) = (n−2)·u(G),     hence     u(n) ≤ ⌊ n · u_max(n−1) / (n−2) ⌋.
```
**Consequence:** the best-known upper bounds for `u(23..27)` — `66, 72, 78, 84, 90`
— are **exactly** this averaging ladder iterated from the single hard bound
`u(22) ≤ 61`. Only `u(21)=57` (exact) and `u(22)≤61` require nontrivial machinery;
everything above `n=22` is one elementary averaging step per `n`.

**(B) Floor-advancement reduction [PROVED].** Writing the averaging step at `m`,
`u(m) ≤ 3m` is implied by `u(m−1) ≤ A(m)` where the exact integer thresholds are
```
        A(25) = 69,     A(26) = 72,     A(27) = 75.
```
Comparing with the construction lower bounds `L(24)=68, L(25)=72, L(26)=76`:
- `u(24) ≤ 69` ⟹ `u(25) ≤ 75` ⟹ **`N* ≥ 26`**. Target gap to construction: `+1`
  (truth `u(24) ∈ [68,72]`; suffices to shave the bound `72→69`).
- `u(25) ≤ 72` ⟹ `u(26) ≤ 78` ⟹ **`N* ≥ 27`**. Gap `0` (needs `u(25)=L(25)=72`).
- `u(27) ≤ 81` would need `u(26) ≤ 75`, but `L(26) = 76 > 75`. **Impossible by
  averaging:** the construction *already exceeds* the averaging threshold, so for
  any true value `u(26) ≥ 76` the averaging step gives only `u(27) ≤ ⌊27·76/25⌋ = 82`.

**(C) Triple confluence at n=27 [VERIFIED arithmetic].** `n=27=3³` is the unique
integer where the kissing line `3n = κn/2` meets the Erdős–Szemerédi–Trotter scale
`n^{4/3}`: both equal `81 = 3⁴`. Equivalently `(3n)³ − n⁴ = n³(27−n)` changes sign
only at `n=27`. With `κ=6` the planar kissing number, `27 = 3^{κ/2} = (κ/2)³`, the
two coinciding **because** `κ/2 = 3`: the size `3^{κ/2}` of the triangle-cube
`K₃^□(κ/2)` that first reaches average degree `κ` (THM-432) equals the line/curve
crossover `(κ/2)³`.

## Proofs

**(A).** Each edge `{a,b}` of `G` is present in `G−v` iff `v ∉ {a,b}`, i.e. for
exactly `n−2` of the `n` deletions; summing, `Σ_v u(G−v) = (n−2)u(G)`. Each
`u(G−v) ≤ u_max(n−1)`, so `(n−2)u(n) ≤ n·u_max(n−1)`, and integrality gives the
floor. The ladder identity `66,72,78,84,90` is verified exactly in
`unit_distance_deletion_ladder_monad.py` (`⌊25·72/23⌋=78`, etc.); each equals the
AMP-cited bound. (Note `u(22)≤61` is itself **one sharper** than the ladder value
`⌊22·57/20⌋ = 62` from `u(21)=57` — that single edge is the only place the heavy
method beats averaging.) ∎

**(B).** `u(m) ≤ ⌊m·u(m−1)/(m−2)⌋ ≤ 3m` ⟺ `m·u(m−1)/(m−2) < 3m+1` ⟺
`u(m−1) < (3m+1)(m−2)/m`, i.e. `u(m−1) ≤ A(m) := ⌈(3m+1)(m−2)/m⌉ − 1`. Direct
evaluation: `A(25)=69, A(26)=72, A(27)=75` (exact integers; checked). The
comparisons with `L(m−1)` are arithmetic. For the impossibility: any planar set
realizes `u(26) ≥ L(26) = 76 > 75 = A(27)`, so the averaging hypothesis `u(26) ≤
A(27)` is **false for the true `u(26)`**; the averaging step from `n=26` yields at
best `u(27) ≤ ⌊27·76/25⌋ = 82`, which does not exclude `u(27)=82`. ∎

**(C).** `27^{4/3} = (3³)^{4/3} = 3⁴ = 81 = 3·27`. The integer polynomial
`(3n)³ − n⁴ = n³(27−n)` is `>0` for `n<27`, `=0` at `n=27`, `<0` for `n>27`, so the
line `3n` and curve `n^{4/3}` cross only at `n=27`. `3^{κ/2}` and `(κ/2)³` are equal
iff `κ/2 = 3` (the equation `3^x = x³` over positive reals has `x=3` as its
relevant root), which holds for the Euclidean plane `κ=6`. ∎

## What this says about N* (the payoff)

The floor side of `N*` **separates into three rungs of decreasing reachability**:

| advance to | needs (exact-value upper bound) | gap to construction | by averaging? |
|---|---|---|---|
| `N*≥26` | `u(24) ≤ 69` | `+1` (vs `L=68`) | reachable *iff* `u(24)≤69` |
| `N*≥27` | `u(25) ≤ 72` | `0`  (vs `L=72`) | reachable *iff* `u(25)=72` |
| `N*≥28` (`u(27)=81`) | `u(26) ≤ 75` | `−1` (vs `L=76`) | **impossible** |

So **proving `u(27)=81` (settling `N*=28`, the conjectured answer HYP-2299) can
never come from the averaging/deletion method** — the construction `L(26)=76`
provably out-runs the averaging threshold `75`. That last rung requires a genuinely
geometric upper bound *at* `n=27` (e.g. ruling out the single value `u(27)=82`
directly). The first two rungs, by contrast, reduce to pinning the *small* exact
values `u(24)` and `u(25)` to (or one below) their Moser-construction values — a
much closer target than the naïve "prove `u(27)≤81`". This is the precise sense in
which "the cost of Euclidean realizability" (THM-431) lives in a handful of `O(1)`
edge counts at `n=24,25,26`.

**Corollary (min-degree necessary condition, rigorous).** Any planar set on `n`
points has unit-distance min-degree `δ ≥ u(n) − u_max(n−1)`: deleting a min-degree
vertex leaves `u(n) − δ ≤ u_max(n−1)`. Hence a `3N`-beating graph is forced toward
near-regularity — e.g. `u(25)=78` ⟹ `δ ≥ 6`; conditionally on construction
optimality (`u_max(26)=76`), `u(27)=82` ⟹ `δ ≥ 6`. This conflicts with the
low-degree convex-hull/boundary vertices that finite planar sets must have, which
is the mechanism, not yet a quantitative bound (the local degree of a hull vertex
is *not* bounded — neighbors fill a semicircular arc — so the obstruction is global,
i.e. the cherry/incidence bound, THM-421).

## Honest scope

- (A),(B),(C) are PROVED (elementary) and exact-integer VERIFIED.
- (A) is a standard averaging identity; the *content* is the identification that
  AMP's small-`n` ladder IS this averaging, and that the heavy method contributes
  exactly one edge (`u(22): 62→61`).
- This does **not** improve any numeric bound; it *localizes* where improvement
  must happen (a non-averaging bound at `n∈{24,25,26}`) and proves the last rung
  (`u(27)=81`) is out of averaging's reach.
- (C)'s `n^{4/3}` strand is a structural *coincidence* (the SST constant `≈1.9 ≫ 1`
  makes the bound vacuous at `n=27`); it is recorded as a benchmark confluence, not
  a bound. See HYP-2367 and the reflection.

## Files
- `04-computation/unit_distance_deletion_ladder_monad.py`
- `05-knowledge/results/unit_distance_deletion_ladder_monad.out`
- reflection `07-reflections/the-unit-distance-floor-is-three-rungs-and-the-last-needs-geometry.md`
- HYP-2367 (the conditional floor ladder)
