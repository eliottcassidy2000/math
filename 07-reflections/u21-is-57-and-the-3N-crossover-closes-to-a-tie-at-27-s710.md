---
source: monad-explorer-2026-06-07-S710 (deep-research; dispatched seed = prove/pin u(21))
status: REFLECTION on THM-431 (which sharpens THM-421). The dispatched n=21 campaign
  is SETTLED by the literature: u(21)=57 is PROVEN (Alexeev-Mixon-Parshall 2024). The
  live mathematics moves one notch over: combining AMP's proven small-N table with
  THM-421 collapses the "first N to beat 3N" interval from [17,32] to [25,28], and the
  best constructions close the deficit to 3N to EXACTLY a tie at n=27 before breaking
  through at n=28.
tags: [unit-distance, erdos, u(21), AMP-2412.11914, THM-421, THM-431, 3N-floor,
  kissing-number, eisenstein, moser-lattice, cartesian-product, wheel, tie-at-27]
---

# u(21) = 57, and the 3N crossover closes to a tie exactly at n = 27

## The seed is already a theorem (cite, don't re-derive)

The campaign "prove/pin `u(21)`" is **closed by the literature**: Alexeev–Mixon–
Parshall (arXiv:2412.11914, Dec 2024) prove the Erdos unit-distance maximum
*exactly* for every `n ≤ 21` — matching the Schade lower bounds with a new
computer-assisted upper bound (forbidden-`F`-free enumeration + a bespoke
unit-distance embedder that dodges cylindrical algebraic decomposition) — and they
fully enumerate the extremal graphs. **`u(21) = 57`**, settling the prior gap
`57 ≤ u(21) ≤ 68`. The honest research move is to *cite* this and ask what it
unlocks elsewhere in the project — not to rebuild a 30-page computer search.

What it unlocks is two things: a clean structure for the extremal graph, and a
sharp tightening of one of our own theorems.

## The extremal graph is a Cartesian product — and the repo had the wrong split

AMP's extremal `n=21` graph is the **generic-angle Minkowski sum of a unit triangle
and a unit wheel**, i.e. the graph Cartesian product

```
   u(21) = e(K₃ □ W₇) = e(T)·n(W) + n(T)·e(W) = 3·7 + 3·12 = 21 + 36 = 57.
```

`W₇` = hub + `C₆` is exactly the **Eisenstein rosette/flower** (6 spokes + 6 rim =
12 edges), the κ=6 sixfold star that runs through the whole unit-distance story.
The S630 reflection had guessed `57 = 20 + 37` (Hamiltonian spine + centered-hex
bulk); the *proven* split is the product split `57 = 21 + 36`. Worth fixing in the
mental model: at the optimum, `n=21` is not "a spine plus a hex shell," it is
"three wheels stacked over a triangle" (or seven triangles stacked over a wheel) —
a **product**, not a core-plus-halo. The product view also explains why `21 = 3·7`
is not a coincidence with forbidden `H=21`: it is `n(K₃)·n(W₇)`, pure combinatorics
of the two factors, nothing to do with Hamiltonian-path counts.

## The real payoff: the 3N floor collapses from [17,32] to [25,28]

THM-421 framed `N*` = the smallest `N` with `u(N) > 3N` (the kissing threshold
`3N = (κ/2)N`) and proved `N* ∈ [17,32]`: a combinatorial common-neighbour floor
of 17, and a `√7` Eisenstein construction at 32. AMP's table pinches both ends:

```
 n    21  22  23  24 | 25  26  27 | 28  29  30
u≤    57  61  66  72 | 78  84  90 | 96 103 110
u≥    57  60  64  68 | 72  76  81 | 85  89  93
3n    63  66  69  72 | 75  78  81 | 84  87  90
```

- **Floor `N* ≥ 25`:** `u(n) ≤ 3n` is *proven* for every `n ≤ 24` (exact for
  `n≤21`; `u(22)≤61<66`, `u(23)≤66<69`, `u(24)≤72=3·24`). No ≤24-point set beats 3N.
- **Ceiling `N* ≤ 28`:** the *realizable* `u(28) ≥ 85 > 84 = 3·28` (Engel's Moser
  lattice). 28 points beat 3N.

So `N* ∈ [25, 28]` — four candidates, down from sixteen.

## The shape of the crossover: deficit closes to a TIE at 27

Track the best-known construction's deficit `u≥(n) − 3n`:

```
 n        22  23  24  25  26  27  28
 u≥ − 3n  −6  −5  −4  −3  −2   0  +1
```

The deficit climbs `−6,−5,−4,−3,−2` one step at a time, then **jumps −2 → 0**
(skipping −1) to a **clean tie at `n = 27`**, and breaks positive at 28. Two things
are striking:

1. **The tie lands exactly at `27 = 3³`.** The "3" of `3N` is `κ/2`, the half-rosette
   of the Eisenstein sixfold symmetry — the same "3" as `χ=3`, the third-roots, and
   the `3` in the LRC worry-set (cf. the-common-neighbour-bound reflection, S630).
   The realizable density reaches *exactly* average-degree-6 (`u = 3·27`) at the cube
   of that 3, then exceeds it. Disciplined reading (per S630): this is suggestive,
   not a license — but `27` being the precise tie point, bracketed by `N* ∈ [25,28]`,
   is the kind of "too-clean" landing the project flags as structure worth probing.

2. **The `√7` family is the wrong lattice for `N*`.** Our own THM-421 ceiling came
   from `√7` Eisenstein disks (first beat at `n≈32`, or `n≈39` for the plain disk —
   re-confirmed by exact `Q=7` recount this session, `U=119>117` at 39). But Engel's
   **Moser lattice** beats 3N already at 28. The `√7` rosette maximises *interior*
   degree (12), but pays too much at the boundary for small `N`; the Moser lattice —
   built from a non-lattice unit-vector set with extra coincidences — wins the
   small-`N` frontier. Lesson: *the asymptotically-densest layer is not the
   first-to-cross layer.* For the moderate-`N` extremal question, boundary economy
   beats interior density, and that selects a different arithmetic (Moser, not `√7`).

## Where to push (HYP-2298)

`n = 27` is the sharp target. The known construction *ties* `3·27 = 81`; is
`u(27) = 81` or `> 81`? An exact-integer Moser-lattice / product construction giving
`82` at `n=27` would prove `N* ≤ 27`; an upper bound `u(27) ≤ 81` would (with `25,26`)
prove `N* = 28`. Either way the crossover index becomes a *named constant* of the
plane — the smallest point count at which Euclidean unit distances overtake the
penny/kissing rate. THM-421's two floors then read cleanly: design-theoretic floor
**17** (Shrikhande-tight), Euclidean-realizable crossover **`N* ∈ [25,28]`**, and the
gap between them is exactly the price of putting the design in the plane.

## The transferable shape

> When an extremal count has a hard asymptotic regime (here: density → 6 in the `√7`
> rosette) and you want the *first crossing* of a threshold, do not optimise the
> asymptotic family. The first crossing is a **boundary-dominated** event: the family
> that minimises perimeter deficit at moderate size wins, even if a different family
> wins at infinity. Two different arithmetics (`√7` Eisenstein vs Moser lattice) own
> the two regimes.
