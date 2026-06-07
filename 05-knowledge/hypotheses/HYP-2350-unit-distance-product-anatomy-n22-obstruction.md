---
id: HYP-2350
name: unit-distance-product-anatomy-n22-obstruction
status: MIXED — Prop. P1 PROVED (elementary); the n=22 resolution stays open
date: 2026-06-07
session: claudebox-2026-06-07-S712
depends_on:
  - THM-431   # u(21)=57 (AMP24); extremal = K3 [] W7
  - THM-412   # density quantization of the unit-distance layer (Eisenstein)
external:
  - "Alexeev-Mixon-Parshall, arXiv:2412.11914 (2024): u(n) exact for n<=21, extremals enumerated."
  - "Engel et al., arXiv:2406.15317: u(22) >= 60 (Moser-lattice / diverse beam search)."
provisional_id: true   # THM/HYP counter contested (THM-421 and THM-431 each already collided); renumber at PR
---

# HYP-2350: The certificate anatomy of u(n), n<=21, is the Minkowski-product structure — and why n=22 is the first open case

## The anatomy (understanding n=1..21)

The proven small-n unit-distance optima are **Cartesian products (generic-angle Minkowski sums)** of a
few small unit-distance "atoms": K2 (edge), K3 (unit triangle), the rhombus R4, and the wheel/flower
W7 = hub + C6 (the Eisenstein rosette, 6 spokes + 6 rim, 12 edges). For a product `G [] H`, vertices
multiply (`|V|=|VG||VH|`), edges follow `e(G[]H) = e(G)|VH| + |VG|e(H)`, and **degrees ADD**
(`deg(g,h)=deg_G g + deg_H h`). The flagship is

```
   u(21) = e(K3 [] W7) = 3*7 + 3*12 = 21 + 36 = 57    (THM-431, AMP24)
```

with degree sequence `{2+6 = 8  (x3 hubs), 2+3 = 5  (x18 rim)}` = the **eighteen-5's-three-8's** of
S711. So S711's measured anatomy IS the product: three degree-8 vertices = the three (K3-vertex, W7-hub)
points, each the center of a unit-circle of 8 neighbors. K_{2,3}-free and omega=3 follow from the
factors. VERIFIED by reconstructing K3 [] W7 from scratch (s712 part B): 21 vtx, 57 edges, 5x18+8x3,
0 K_{2,3}, matching both S711 and the stored AMP "core 1".

## P0 (product lower bound, classical). 
`u(ab) >= u(a)*b + a*u(b)`, realized by `A [] B` at a generic relative angle (no edge collisions).
Tabulating `u_prod(n) = max over factorizations` vs the true `u(n)` (s712 part A) gives the
**product-optimal composites n<=21 = {6, 8, 9, 12, 21}** (gap 0), the rest having a small positive gap.
The optimum is a clean 2-factor product exactly at these n; 21 = 3*7 is the largest, the last proven
value, and product-optimal.

## P1 (the n=22 product obstruction — PROVED, elementary).
**If `u(22) >= 58`, then no extremal 22-point set is a Minkowski sum `A (+) B` of two sets each of size
`>= 2`.**
*Proof.* A 22-point Minkowski sum with `|A|=a, |B|=b, ab=22` forces `{a,b} = {2,11}` (22 = 2*11 is the
only factorization into factors `>= 2`). To have exactly 22 points the sum must be collision-free
(generic), so its unit-distance count is exactly `e(A)*11 + 2*e(B) <= u(2)*11 + 2*u(11) = 1*11 + 2*23 =
57 < 58`. A non-generic angle creates point-collisions, dropping below 22 points. Hence any 22-point set
with `>= 58` unit distances is not such a Minkowski sum. ∎

This is the clean structural reason 22 is the first awkward case: it is the first frontier value whose
**only** factorization (2*11) forces the weak K2-doubling against an indecomposable prime 11, so the
optimum (`>= 60 > 57`) is necessarily **non-product** — it must come from an additive (lattice / Moser-
ring) construction, not a multiplicative one. 21 (multiplicative, solved) -> 22 (additive, open) is a
genuine arithmetic phase change, not just "the computation got harder."

## P2 (single-vertex extension of the maximum 21-graph — VERIFIED negative).
Adding one point at unit distance from `k` existing points requires `k` of them to lie on a common
unit circle. Searching this over the entire 1-parameter K3 [] W7 family (2000 angles, near-degenerate
angles discarded — see Caution) and on the stored AMP "core 1": the **maximum single-vertex extension
degree is 2** (`57 + 2 = 59`), never 3 or 4. So:
- the record `u(22) >= 60` is **NOT a single-vertex extension of the maximum 21-graph** — it is a fresh,
  non-product configuration (consistent with Engel's Moser-lattice 60);
- by the deletion argument (S614): a `61`-edge UDG on 22 vtx minus a min-degree-(4 or 5) vertex leaves a
  `56` or `57`-edge 21-core. P2 shows the `57`-core route (extend the maximum graph by a degree-4
  vertex) is **closed**. Therefore **if `u(22) = 61`, its 21-core is a `56`-edge graph extended by a
  degree-5 vertex** — the search should target the 56-edge UDGs, not the maximum ones.

## Caution (a methodological catch worth recording).
The first (buggy) sweep reported a "degree-4 extension => u(22) >= 61" at `theta ~ 0.0003`. This was a
**float artifact**: as `theta -> 0` the product degenerates toward an aligned triangular-lattice patch
with near-colliding points, and a fixed `1e-7` tolerance counts phantom unit distances (the tell: it
also reported impossible "60-edge 21-cores", exceeding the proven max 57). Guarding on minimum pairwise
separation `>= 0.05` removes the entire degenerate regime; all 1808 surviving angles give exactly 57
core edges and max extension degree 2. Near-degenerate angles must never be trusted in float UDG search.

## Next tests
- run P2's extension search on the OTHER AMP-enumerated 57-edge extremals (not just K3 [] W7);
- obtain the 56-edge 21-UDGs and search degree-5 extensions (the only remaining route to 61);
- the LRC transfer (S710 3N reflection): is the LRC tight family likewise a suboptimal *multiplicative*
  (AP/shell) rung beaten by an *additive* (V*-type) configuration — the same product-vs-additive split?
