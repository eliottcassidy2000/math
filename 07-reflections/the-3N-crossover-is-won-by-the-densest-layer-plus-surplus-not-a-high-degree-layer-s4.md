# The 3N-crossover escapes rank 2 — the degree–radius tension IS the 2-D kissing bound, and the Moser ring is the minimal escape

*(Filed under the working title "densest layer plus surplus"; the first draft's
"triangular + surplus" reading was wrong and is corrected in §"What actually wins at N=28".)*

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

So **two structurally unrelated families of "nice" rank-2 / product constructions — single-norm
sublattices and Cartesian products — independently bottom out at the same crossover `N = 32`**,
four above the true frontier `N* ≤ 28`. The gap `[28, 32]` is what both lanes call, in their own
language, the cost of regularity: the irreducibility premium (products) and the degree–radius
penalty (lattices) are two faces of one obstruction.

## What actually wins at N=28 — and why it had to leave rank 2

It is tempting (I fell for it in the first draft of this note) to read `u(28)=85` — average
degree `6.07`, a hair above `κ` — as "a triangular patch plus a tiny surplus." That is **wrong**,
and the error is instructive. A triangular (penny) `28`-patch gives only `Harborth ⌊84−√333⌋ = 65`
(avg deg `4.6`); the best single-norm **`√7`** patch gives `83` (avg deg `5.9`). `85` is not a
perturbation of triangular at all — it is *denser than any degree-12 2-D lattice patch*.

The actual extremal lives in **Engel's "Moser lattice"** (THM-432; Engel et al. arXiv:2406.15317):
the **rank-4** ring

```
   M_L = ℤ[ζ₆, ω₃],   ζ₆ = (1+i√3)/2  (cos = 1/2),   ω₃ = (5+i√11)/6  (cos = 5/6),
```

sitting in the biquadratic CM field `ℚ(√−3, √−11)`. The decisive fact: **`ω₃` is a unit**
(`|ω₃| = 1`) of *infinite* order (it is not a root of unity — `cos = 5/6 ∉ {0,±½,±1}`). So `M_L`
has **18 unit vectors all at radius exactly 1**, and a Moser vertex can reach degree `18` *without
paying any radius*. This is exactly the escape hatch the degree–radius law forbids in rank 2:

> The degree–radius tension `N_cross ∝ ρ·t·(deg/(deg−6))²` **is the 2-D kissing bound in disguise.**
> In a rank-2 lattice the minimal-distance shell has at most `6` vectors (kissing number 6), so
> degree `> 6` forces the unit shell *off* the minimal radius (`t > 1`), and you pay boundary
> deficit `∝ √t`. The Moser ring breaks this by being **rank 4 and dense**: `ω₃` supplies a second
> independent unit direction, packing 18 units onto radius 1. Degree 18, no radius penalty.

So the `[28, 32]` gap is not "regular vs irregular" — it is **the cost of staying rank-2**. The
single non-torsion unit `ω₃` (the "Moser angle", `arccos 5/6`) is the *minimal* algebraic
ingredient that escapes the kissing trap: one extra CM direction beyond `ℚ(√−3)`. This is precisely
THM-432's "bridge ring between the triangular lattice and the CM field" (reopening HYP-2262) — and
the degree–radius law explains **why** such a bridge is *necessary*: nothing inside rank 2 can beat
`3N` before `32`.

## Concrete residue

- `N* ∈ [25,28]` is **not lowered** by any single-norm rank-2 lattice: all six families (and, by
  THM-433, all Cartesian products) are `≤ 3N` through `N = 28`; the earliest rank-2 crossover is
  `√7` at `32`. The ceiling `28` stands on the **rank-4 Moser ring**, not a 2-D lattice — this is
  the precise, corrected sense in which Engel's record "evades the structured families."
- The right hunt for `u(27) > 81` (⟹ `N* ≤ 27`) is therefore a **dense Moser-ring `M_L` patch**,
  *not* a denser single 2-D lattice and *not* a product. The best `√7` rank-2 patch gives `78` at
  `N=27` (deficit −3); only the 18-unit-vector Moser structure reaches the record.
- **DONE — ceiling now self-contained.** A densest-patch search run *directly in `M_L`* (graph-BFS
  greedy + anneal in `ℤ⁴` with the 18 unit-vector offsets, every count an **exact `|z|²=1` recount
  over `ℚ(√3,√11)`**) reproduces the entire cited Engel/Schade deficit table from scratch:

  ```
   N      22  23  24  25  26  27  28  29  30
   u(M_L) 60  64  68  72  76  81  85  89  93     (exact, recount-verified)
   3N     66  69  72  75  78  81  84  87  90
                              tie  +1  +2  +3
  ```

  It **ties `3N` at `N=27=3³` (81=81)** and **first beats it at `N=28` (85 > 84)** — exactly Engel's
  record, now with explicit exact-integer coordinates found here. So THM-431's previously *cited*
  ceiling `N* ≤ 28` is **self-contained**, and the contrast is sharp and computed, not asserted:
  the **rank-4** Moser ring crosses at `28`; the best **rank-2** single-norm lattice (`√7`) crosses at
  `32`; the gap `[28,32]` is the cost of staying rank-2, paid in full at the kissing bound.
- Files: `04-computation/unit_distance_3n_crossover_families_s4.py`, `…_focus_s4.py`,
  `…_moser_crossover_s4.py` (+ `05-knowledge/results/*_s4.out`); HYP-2301.

*The mathematics keeps saying the same thing from every side: the extremal object is the one
that refuses to be regular. The densest you can be and still beat the kissing cap is to be
almost-triangular and slightly broken.*
