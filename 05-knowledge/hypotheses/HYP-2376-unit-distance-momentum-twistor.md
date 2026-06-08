---
id: HYP-2376
name: unit-distance-momentum-twistor-angle-is-log-of-rotation
status: VERIFIED (twistor linearization, direction decomposition, CM layers, n=22 search) + a unification conjecture
date: 2026-06-08
session: claudebox-2026-06-08-S719
depends_on:
  - THM-445   # LRC momentum twistor = discrete log (S718): the parallel construction
  - HYP-2350  # unit-distance product anatomy (u(21)=K3 [] W7), n=22 obstruction (S712)
  - THM-441   # convolution-correlation-adjoint duality / autocorrelation operator (S714)
  - THM-431   # u(21)=57, AMP sqrt(13) layer
external:
  - "Alon-Bloom-Gowers-... grid-disproof: CM fields, bounded-height modulus-1 elements (S614 memory)"
provisional_id: true
---

# HYP-2376: The unit-distance momentum twistor is the angle map (the log of the rotation group)

The LRC momentum twistor (THM-445) is the discrete log because the LRC dual conformal symmetry is the
MULTIPLICATIVE group `(Z/m)*`. The unit-distance dual conformal symmetry is GLOBAL ROTATION — the
multiplicative group `U(1)` of modulus-1 numbers, since a unit edge IS a difference vector on the unit
circle. So the unit-distance momentum twistor is the **angle map** `e^{i theta} -> theta` (the log of
`U(1)`), and it is the exact analog of the LRC discrete log.

## (1) VERIFIED: the angle twistor linearizes the rotation (dual conformal) symmetry

Global rotation of a configuration by `phi` shifts every unit-edge angle by `phi` (verified exactly on
`K3 [] W7`): in the twistor (angle) coordinate the dual conformal symmetry is the additive translation
`theta -> theta + phi` — exactly as the LRC multiplier becomes `ell -> ell + ell(a)`.

## (2) VERIFIED: the angle-twistor autocorrelation u = sum over directions

Resolve unit distances by edge LINE-direction `phi in [0,pi)`: `u(config) = sum_phi m(phi)`, where
`m(phi)` = #unit edges along direction `phi`, and `m(phi) <= n-1` (a single direction is a union of
unit-vector chains). The proven optimum `K3 [] W7` (`u=57`) has **exactly 6 line-directions** with counts
`[12,12,12,7,7,7]` (sum `=57`): three directions carrying `7 = |W7|` (the `K3` edges replicated over the
7 wheel vertices) and three carrying `12 = e(W7)` (the wheel edges replicated over the 3 triangle
vertices). The product structure is read directly off the twistor: `57 = 3*7 + 3*12`. A triangular
22-blob has only 3 line-directions and `u=49`.

**The maximization principle:** `u` is large iff the unit edges CONCENTRATE into FEW directions, each
near-saturated — i.e. the edge-directions form a ROTATION GROUP (roots of unity). This is the twistor
explanation of why optimal unit-distance configs are CM/cyclotomic (the grid-disproof's modulus-1 supply):
maximizing the angle-autocorrelation peak forces a multiplicative (rotation-group) structure on directions.

## (3) VERIFIED: CM layers; r_Q(D) = degree = rotation-group order = the twistor Z/r_Q(D)

In the Eisenstein layer where the "unit" is the norm-`D` vector, the number of unit-neighbours (degree) is
`r_Q(D)` = #representations of `D` by `a^2+ab+b^2`, and the rotation group has order `r_Q(D)` (its log is
`Z/r_Q(D)` — the twistor). `r_Q(1)=6` (triangular, deg 6), `r_Q(7)=r_Q(13)=...=12` (deg 12 — the AMP
`sqrt(13)` layer, THM-431/S710), `r_Q(49)=18`, `r_Q(91)=24`. The optimum uses `sqrt(13)` (deg 12), not
`sqrt(1)` (deg 6): more directions = denser. The discrete log of the rotation group is the unit-distance
analog of the LRC discrete log.

## (4) BLEEDING EDGE (n=22): a sqrt(13) construction reaches u(22) >= 57, beating the triangular 49

Densest 22-point subset found by greedy+swap in Eisenstein norm-`D` layers: norm-1 (triangular) `49`,
norm-3 `49`, norm-7 `48`, **norm-13 `57`**. So the `sqrt(13)` CM layer gives an explicit `u(22) >= 57`
(equal to `u(21)`, far above the triangular `49`), confirming the twistor thesis (the deg-12 rotation group
beats the deg-6 one). The record `u(22) >= 60` (Engel) and the open `61` exceed any single lattice patch
because the true optimum is a **non-lattice rigid** graph (S710: lattice DISK crosses `3N` only at `N=43`,
true optimum at `N=28`). Twistor reading: reaching `60/61` needs MORE saturated directions than a lattice
patch supplies at `n=22` — a rigid multi-direction (CM-rotation) NON-PRODUCT graph (S712: `22=2*11` caps
products at `57`).

## The unification conjecture (LRC twistor = unit-distance twistor)

Both twistors are the logarithm of the multiplicative dual conformal symmetry:
| | LRC (THM-445) | unit distance (here) |
|---|---|---|
| dual conformal group | `(Z/m)*` (multiplier) | `U(1)` (rotation) |
| twistor (log) | discrete log -> `Z/phi(m)` | angle -> `R/2piZ` (CM: `Z/r_Q(D)`) |
| special conformal | inversion `v->v^{-1}` -> negation | reflection `theta -> -theta` |
| extremal/critical structure | half-system (tiles `Z/phi`) | edges fill a rotation orbit (CM module) |
| the answer | dodge = additive covering | `u` = angle-autocorrelation peak |

**Conjecture:** the LRC transversal core (half-system tiling of `Z/phi(m)`) and the extremal unit-distance
configs (edge-directions a rotation orbit / CM module) are the SAME phenomenon — a difference set filling a
coset structure of the rotation group — read in the two twistors. The grid-disproof's CM modulus-1 supply
is the unit-distance avatar of the LRC cyclotomic shell tower (THM-439).

## Next
- exact-arithmetic `sqrt(13)` rigid (non-lattice) 22-graph search toward 60 (the AMP/Engel regime);
- prove `m(phi) <= n-1` gives a usable upper bound when combined across the (few) CM directions;
- the unification conjecture: match half-systems in `Z/phi(2n-1)` to rotation-orbit unit-distance configs.
