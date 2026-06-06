---
id: HYP-2264
status: CONFIRMED (computational; global claim CONJECTURE)
source: monad-explorer-2026-06-06-S703
related: [THM-411, HYP-2262, HYP-2257, T755, T756]
---

# HYP-2264 — Triangular lattice is the finite-crossover minimizer

## Statement

Among all 2D lattices (positive-definite binary quadratic forms `Q`), the
**triangular lattice** (`Q = x^2+xy+y^2`, disc `-3`) minimizes the finite
crossover

```text
   N_cross(Q) = min { N : some N-point patch of Q has more than 3N unit distances },
```

with `N_cross = 43` (for `> 3N`) and `61` (for `> 3.5N`). Equivalently, the
triangular lattice beats Harborth's `3n` floor with the fewest points.

## Mechanism (THM-411 density quantization)

Density `r_Q(D)/2` is a multiple of `w/2` (`w` = # proper automorphs = roots of
unity). Triangular (`w=6`) is forced to density `6` (skips 4,5) but reaches it at
the smallest popular norm `D* = 7`. Boundary-lens law:
`N_cross ~ c^2 D* density^2/(density-3)^2`; triangular minimizes it (model `28`
vs square `80`, disc -15 `100`, disc -12 `112`).

## Evidence

`04-computation/lrc_density_quantization_crossover_s703.py` (exact integer
arithmetic, disk patches `Q(g) <= T`):

| form | disc | density | `D` | `N(U>3N)` | `N(U>3.5N)` |
|------|------|---------|-----|-----------|--------------|
| (1,1,1) triangular | -3 | 6 | 7 | **43** | **61** |
| (1,1,4)            | -15| 5 | 16| 71 | 117 |
| (1,0,3)            | -12| 6 | 28| 79 | 115 |
| (2,1,2)            | -15| 6 | 32| 83 | 115 |
| (1,1,2)            | -7 | 5 | 16| 97 | 175 |
| (1,0,1) square     | -4 | 4 | 5 | 101 | 421 |

Density quantization `w | r_Q(D)` verified for ALL reduced primitive forms to
disc `< -200`.

## Patch-shape robustness (S703b)

`04-computation/lrc_crossover_patch_robustness_s703b.py`: triangular wins under
ALL three patch families (the main caveat below is thereby resolved):

| form | disk `N` | box `N` | euclidean `N` |
|------|----------|---------|----------------|
| triangular | **43** | **49** | **43** |
| square     | 101 | 169 | 101 |
| disc-12    | 79  | 121 | 79  |
| disc-15 (density 5) | 71 | 81 | 71 |

(`disk` = Q-ball, `eucl` = true round Euclidean disk in the embedded plane —
these agree exactly since `Q(i,j)` IS the squared Euclidean length; `box` =
index-box, the genuinely different shape.)

## Honest limits

- CONFIRMED for all competitive forms (smallest `D*`/model) by exact enumeration,
  ROBUST across disk / box / euclidean patches (S703b).
- The **global** minimality (over EVERY lattice) is a well-supported CONJECTURE:
  only small-model forms were exact-checked. No form came within `28` of
  triangular's model and none beat its exact `43` under any patch.
- Corrects S702: square crosses at `N=101` (radius `sqrt5`), not `121`.
