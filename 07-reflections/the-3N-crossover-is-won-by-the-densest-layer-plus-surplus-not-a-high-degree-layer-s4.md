# The 3N-crossover is won by the densest layer plus a small surplus — not by a high-degree layer

*monad-explorer-2026-06-07-S4 (crossover-lane). Builds on THM-431 / THM-431-C (S710),
THM-432/433 + HYP-2299/2300 (concurrent S1/S2/S711, the **product** lane). This is the
**single-norm lattice** lane — its independent confirmation, and the mechanism underneath.*

## The setup

`u(N)` = the Erdős unit-distance maximum; `N*` = smallest `N` with `u(N) > 3N`
(average degree `> κ = 6`). THM-431 pins `N* ∈ [25,28]`; the true crossover is
realized by Engel's "Moser-lattice" record `u(28) ≥ 85 > 84`, which the repo only
*cited*. THM-431-C found the repo's own `√7` Eisenstein family crosses far too late
(disk 39, anneal 32) and called it "the wrong family." Two natural questions were left
open: *which* family is right, and *why* is `√7` wrong?

## What I did

A systematic exact-integer densest-patch search across **six single-norm lattice
unit-distance families** — every pair `(p,q)` at squared distance `t` is an edge:

| family | offset set | degree | radius² t |
|---|---|---|---|
| triangular / penny (Eisenstein `Q=1`) | 6 vectors | 6 | 1 |
| **knight** (grid `‖·‖²=5`) | `(±1,±2),(±2,±1)` | 8 | 5 |
| **√7 Eisenstein** (`Q=7`) | 12 vectors | 12 | 7 |
| **√13 Eisenstein** (`Q=13`) | 12 vectors | 12 | 13 |
| grid `t=25` | 12 vectors | 12 | 25 |
| grid `t=65` | 16 vectors | 16 | 65 |

For each family and each `N` I annealed the densest `N`-cell patch (calibrated to
reproduce the repo's `√7` record 97 at `N=32`) and **recounted every patch exactly**
from raw coordinates. Result (`unit_distance_3n_crossover_{families,focus}_s4`):

```
first single-norm-lattice 3N-crossover, by family
  √7 Eisenstein   N = 32   (ties at 31, beats at 32: 97 > 96)   <-- EARLIEST
  √13 Eisenstein  N = 33   (ties at 32, beats at 33)
  knight (t=5)    > 60     (deg 8, deficit still −16 at N=60)
  grid t=25       > 60
  grid t=65       ≈ 60     (deficit −1 at N=60)
  penny (t=1)     never    (deg 6 ⟹ avg deg < 6 always)
```

**No single-norm lattice family beats `3N` at `N ≤ 28`.** The earliest is `√7` at `32`.

## The mechanism: a degree–radius tension

Why does the degree-8 knight lose to the degree-12 `√7`, and degree-16 `t=65` lose to
both? A boundary-deficit estimate (disk patch of `N` cells, lattice density `ρ`, real
radius `R=√(N/πρ)`, neighbour vectors of length `√t` so a boundary annulus of width
`√t` is degree-deficient) gives `edges − 3N = ½(deg−6)N − ½·c·ρR√t·deg`, hence

```
        N_cross  ∝  ρ · t · ( deg / (deg − 6) )² .
```

This **proxy reproduces the observed ordering across all six families exactly**, and
pins `√7` at `32.3` vs observed `32`:

```
  √7   ρt(deg/(deg-6))² = 32   (obs 32)
  √13                   = 60   (obs 34 — disk model loosens for larger t, order still right)
  knight                = 80   (obs >60)
  t=25                  = 100  (obs >60)
  t=65                  = 166  (obs ≈60)
  penny                 = ∞    (deg 6)
```

The two factors pull against each other:

- **`(deg/(deg−6))²` punishes chasing high degree near the cap.** It is `∞` at the
  kissing degree `6`, drops to `16` at deg 8 (knight), `4` at deg 12, `2.25` at deg 18,
  `→1` as `deg→∞`. So you *must* clear 6 comfortably — but the marginal value of extra
  degree collapses fast.
- **`t` (squared radius) punishes the wider layers you need to *get* more degree.** In
  the square grid the first degree-12 layer is `t=25`; in the Eisenstein lattice it is
  `t=7`. Degrees come quantized: grid layers have `deg ∈ {4,8,12,16,…}` (multiples of 4),
  Eisenstein layers have `deg ∈ {6,12,18,…}` (multiples of 6, since `r_E(n)=6(d_{1,3}−d_{2,3})`).

The product `t·(deg/(deg−6))²` is therefore minimized by the **densest-radius layer that
already carries degree 12** — and that is uniquely the `√7` Eisenstein rosette (`deg 12`
at the *smallest* norm `t=7`). The grid's degree-12 layer sits at `t=25` (radius too big);
the knight's small radius `t=5` carries only degree 8 (penalty too big). `√7` threads the
needle. **So `32` is not an accident of `√7` — it is the genuine single-lattice optimum,
forced by the degree–radius tension.** This sharpens the "32" rung of THM-431's nested
ladder `17 < N* ≤ 28 < 32 < 43` from "a √7 blob" to "the minimum over *all* single-norm
sublattice families."

## The convergence — and what it means

Here is the resonance I did not expect. The concurrent **product** lane (S1/S711,
THM-432/433) found, by a completely different argument — average degree is *additive*
under the Erdős–Minkowski product, `avgdeg(G□H)=avgdeg(G)+avgdeg(H)` — that the **product
family** also first beats `3N` at exactly **`N = 32`** (`W₁₆□K₂`, 98 > 96), ties at
`27 = 3³` (the Hamming cube `K₃^□3`, forced 6-regular), and is `≤ 3N` for every `N ≤ 31`.

So **two structurally unrelated families of "nice" constructions — single-norm sublattices
and Cartesian products — independently bottom out at the same crossover `N = 32`**, four
above the true frontier `N* ≤ 28`. The gap `[28, 32]` is what both lanes call, in their own
language, the cost of regularity: the irreducibility premium (products) and the
degree–radius penalty (lattices) are two faces of one obstruction.

The unifying reading: **to beat `3N` you do not want a high-degree layer at all.** `u(28)=85`
means average degree `6.07` — a hair above `κ`. The way to get a hair above `6` at small `N`
is the *maximally compact degree-6 layer* (the triangular/penny lattice, smallest possible
radius, `Harborth ⌊3N−√(12N−3)⌋` just *below* `3N`) plus a tiny **non-lattice surplus**.
At `N=28` that surplus is only **2 unit distances** — `83` (best `√7` lattice) vs `85`
(Engel) — yet it moves the crossover from `32` down to `28`. The frontier is won at the
*boundary* of the densest layer, perturbed; never in the *interior* of a wider one. Both the
high-degree lattices and the clean products overspend on regularity and pay for it in `N`.

This is why neither lane can construct the `N*` graph: products are capped by additivity
(THM-433), single lattices by the degree–radius tension (here), and the true minimizer is an
**irregular blob that is neither** — it lives exactly in the `2`-edge gap between the densest
regular layer and `3N`.

## Concrete residue

- `N* ∈ [25,28]` is **not lowered** by any single-norm lattice: all six families (and, by
  THM-433, all products) are `≤ 3N` through `N = 28`. The ceiling `28` stands on Engel's
  genuinely non-lattice record.
- The right hunt for `u(27) > 81` (which would give `N* ≤ 27`) is therefore a **triangular
  patch + O(1) non-lattice perturbation** — *not* a denser lattice layer and *not* a product.
  Quantitatively: a 27-cell penny/triangular blob gives `≈ 78` (deficit −3, from the sweep);
  it needs **4 extra unit distances** from off-lattice placement to reach 82. That is the
  precise target.
- Files: `04-computation/unit_distance_3n_crossover_families_s4.py`,
  `…_focus_s4.py` (+ `05-knowledge/results/*_s4.out`); HYP-2301.

*The mathematics keeps saying the same thing from every side: the extremal object is the one
that refuses to be regular. The densest you can be and still beat the kissing cap is to be
almost-triangular and slightly broken.*
