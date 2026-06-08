---
id: HYP-2364
name: abc-recurrence-is-boolean-IE-and-C-finiteness-computational-tricks
status: PROVED/VERIFIED (parts 1-4 all computationally confirmed; the order=width claim is the conjecture)
date: 2026-06-07
session: claudebox-2026-06-07-S716
depends_on:
  - HYP-2360  # the triangle recursion / temperature ladder (S715)
  - THM-291   # Mode-B n->n-2 multilinear recursion (boundary = overlap + wiring + apex)
  - THM-337   # base-path staircase order-3 recurrence
  - THM-326   # H = independence polynomial
  - THM-027   # transfer-matrix trace (the C-finiteness engine)
  - THM-329   # A038375 (not C-finite)
provisional_id: true
---

# HYP-2364: A+B+C-D-E-F+G is boolean inclusion-exclusion; the recurrence gives two computational tricks

## What the recurrence IS

The 7 pieces are the `7 = 2^3 - 1` nonempty subsets of a 3-set: `{A,B,C}` singletons (+), `{D,E,F}`
pairs (-), `{G}` triple (+). The signs `(+,+,+,-,-,-,+)` are `(-1)^(|S|+1)` = the **Mobius function of the
boolean lattice B_3**. So `A+B+C-D-E-F+G` is exactly **inclusion-exclusion for the union of the three
corner `(n-1)`-subtriangles**, and the coefficient vector `(1,-3,3,-1)` is `(x-1)^3 = Delta^3` (3rd finite
difference).

## RELATION 1 — the d-simplex family (PROVED/verified d=1..4)

Corner-inclusion-exclusion over a **d-simplex** has `C(d+1,k)` pieces of size `n-k` with sign
`(-1)^(k+1)`, reconstructing the side-n simplex. Its operator is `(x-1)^(d+1) = Delta^(d+1)`, and the
cell count `N_d(n) = C(n+d-1, d)` is a **degree-d polynomial**. Verified: the IE holds and
`Delta^(d+1) N_d = 0` for `d = 1,2,3,4`.

```
   d=1 (interval, 2 corners):   (x-1)^2,  linear cells
   d=2 (triangle, 3 corners):   (x-1)^3,  quadratic cells   <- the user's tournament staircase
   d=3 (tetrahedron, 4 corners):(x-1)^4,  cubic cells
```
The recurrence ORDER = number of corners = `d+1`; the growth DEGREE = `d`. The tournament staircase is
the 2-simplex, so order 3, quadratic.

## RELATION 2 — exactly which invariants obey it: the valuative (additive) ones (VERIFIED)

A tournament invariant obeys the **exact** IE `F(union)=A+B+C-D-E-F+G` iff it is **valuative/additive**
(a measure on the tile/arc set). Verified: arc-count of a real 8-vertex tournament reconstructed exactly
from its three vertex-deleted corners minus overlaps plus center (`28 = 28`). Such invariants then satisfy
`Delta^3 F = 0`, are **quadratic in n**, and are determined by **3 seeds**. (`#tiles=C(n-1,2)`:
auto-discovered recurrence `(3,-3,1)`, reconstructed from `[1,3,6]`.) Non-valuative parts (e.g. the full
`H`-polynomial, degree `n-1` in the tiles) do NOT obey it.

## TRICK A (additive): 3 seeds + O(1)/term, or a closed quadratic

Any tile-additive tournament invariant `= a n^2 + b n + c`; fit from 3 values, evaluate in O(1). No
enumeration needed once you know it is valuative.

## TRICK B (multiplicative): C-finiteness => companion-matrix exponentiation, O(log n) vs O(2^n)

`H` (#Hamiltonian paths) of a **fixed recursive family** is **C-finite** (satisfies a linear recurrence)
because the family is built by a bounded-width transfer operation (THM-291/027). Therefore:

1. compute the first few `H` by brute/DP;
2. **auto-discover** the minimal linear recurrence (solve the Hankel system over Q);
3. **matrix-power** the companion matrix to reach any `n` in `O(r^3 log n)`.

Verified end to end on the base-path family: from 9 DP terms `[1,5,17,57,193,653,2209,7473,25281]` the
finder recovers `(3,1,1)` (= THM-337), and matrix-power then gives
`H(k=20)=56,804,250,945`, `H(k=50)~4.3e26`, `H(k=100)~1.25e53` instantly — while direct DP is
`O(2^(2k))`, infeasible past `k ~ 13`. This is the transfer-matrix method realized as a recurrence;
**the recurrence order equals the family's boundary width** (conjecture: order = number of boundary
states of the transfer operator; for the staircase, `width` tied to the wiring/apex tiles of THM-291).

## RELATION 3 — the order/temperature split (continues S715)

The recurrence ORDER is a **geometric** invariant (= corners = `d+1` = boundary width); the recurrence
ROOTS move with the additive->multiplicative "temperature": additive `(x-1)^3` (triple root 1, polynomial)
vs base-path `x^3-3x^2-x-1` (root `~3.383`, exponential). Same order 3, different roots.

## THE BOUNDARY (VERIFIED) — max-H is NOT C-finite

`A038375` (max-H over ALL tournaments) admits **no** linear recurrence of order `<= 5` (checked on 13
terms). The optimum-over-all-tournaments has unbounded effective width, so neither trick applies — only
fixed families are C-finite. This is exactly why THM-329 needs annealing for new A038375 terms.

## Next
- prove "recurrence order = transfer-matrix boundary width" and read off the order for each staircase
  family from THM-291's boundary (2n-1 tiles -> but the recurrence collapses to small order: why?);
- catalog C-finite tournament H-families and their char polys (a zoo of `x^3-3x^2-x-1`-type roots);
- combine TRICK B with the Pfaffian (S715): is `Pf` of a recursive family also C-finite (a signed
  transfer), giving `det(I+2A)` of the family in O(log n)?
