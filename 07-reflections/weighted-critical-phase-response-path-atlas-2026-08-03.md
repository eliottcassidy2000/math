# Weighted critical phase and the response-path geodesic backbone

**Date:** 2026-08-03
**Status:** promoted as
[THM-3277](../01-canon/theorems/THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas.md),
`PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`. The first hostile audit confirmed every finite graph
census but caught two load-bearing statement defects: the optimized objects
are vertex-simple paths, not arbitrary walks, and the primary phase gauge is
the full-Jacobian generator `[16-17]`, not the secondary singleton generator
`[7-17]`. Both were repaired below, and a fresh independent replay passed
without further correction.

## Inheritance pass and the reframe

The closest proved mechanism is `THM-3269`: its exact analytic edge weight
selects an ordered tree pair, the weak tree's center row 17, and the unique
incident primitive full-Jacobian direction `[16-17]`. `THM-3273` supplies the
intrinsic order-twelve quotient. The canonical hostile example is their
six-class vertex image: twelve rows are not a `C_12` torsor. The corrected
near miss is `THM-3260`'s rank-eleven bridge: a dimension match plus a
bispanning chart is not a physical cyclic action.

The live concept board was:

1. critical divisor difference rather than vertex label;
2. clutch product rather than unweighted hop length;
3. weak tree versus strong cotree;
4. symmetric exchange with a preserved target;
5. the quotient loss from endpoint pairs to one phase difference.

The bold question was whether the missing cyclic data should be sought as a
vertex torsor at all.  A weaker but physically better-typed target is an
oriented response path whose endpoint divisor difference is prescribed.  On
this target, the analytic and critical data lock together unexpectedly well.

## Audited dependency boundary

The declared companion first passed normal, optimized, and stored replay
byte-for-byte.  The new companion then rebuilt the result by a different
route:

- it loaded only the pre-main exact coefficient definitions of `THM-3238`;
- reconstructed all 22 integral response rows on the eleven reset-link
  states;
- solved each one-edge trapping inequality directly, without importing
  `THM-3254`'s interval or cover functions;
- recomputed every overlap strength and its maximizing witness pairs;
- checked reversal and a nontrivial exact positive row rescaling;
- enumerated all `C(22,11)` subsets, obtaining all 74,748 spanning trees,
  and only then filtered the 9,920 ordered complementary-tree charts;
- recovered the unique minimizer, its complement, exact second gap and
  broken intrinsic involution;
- independently recovered rigidity, center 17, radius four, diameter eight,
  the four leaves, both incidence determinants and the integral transition
  matrix.

The independent strength-bank and transition digests are respectively

```text
c7b8441f54bd893f33f1767db913da0d40f392f20d74f6d177782a0490dc911a
93f9ae1643ad9606f15826d8f1d3809fd0b04ce1d74c11d92d354edbdd737d4f
```

and the second-to-first chart ratio is again

```text
111771229679330735594490655 / 106891942612503954628585203.
```

No assertion or scope repair was found in the polarization, root, full
critical ordering, and abstract augmentation bridge used here. A later
six-edge unimodular sharpening of THM-3269 is a separate delta and is not a
dependency of the path minimization below.

## Gauge ledger: full-Jacobian primary, singleton generator secondary

Let

```text
J_12 = Jac(G_0) / 12 Jac(G_0),
```

the intrinsic cyclic order-twelve quotient from `THM-3273`.  The unique
weighted weak tree `T_*` has unique center

```text
r = 17.
```

Use the rooted Abel--Jacobi map

```text
A_r(v) = [v-r] in J_12.
```

THM-3269's unique incident primitive full-Jacobian direction fixes the primary
gauge

```text
g_full=[16-17] -> 1 in J_12.                            (A)
```

An independent reduced-Laplacian and adjugate reconstruction also shows that
the generator-valued classes in `A_r(V)` have precisely two fibres:

```text
one order-twelve class has fibre size 1,
one order-twelve class has fibre size 3.
```

Therefore the singleton class is a second intrinsic generator of `J_12`. Its
unique vertex is row 7, and in the primary gauge

```text
[7-17]=7[16-17].                                        (B)
```

Both characterizations are intrinsic, but they are not the same gauge.
Equation `(A)` is primary because it is the image of THM-3269's canonical
full-critical generator. Equation `(B)` is retained as a useful secondary
normalization and a hostile check against silent multiplication by a unit.

The resulting canonical internal response-phase chart is

```text
row:     2   3   7  10  11  13  16  17  18  19  21  22
phase:  10  10   7  10   4   4   1   0   1   4   8   1.
```

Its fibre multiplicities are still

```text
0:1, 1:3, 4:3, 7:1, 8:1, 10:3.
```

The twelve vertices still emphatically fail to form a torsor. The path target
below is only their relative divisor difference in the primary gauge.

## Candidate 1: every phase has a canonical weighted vertex-simple path

For an oriented core edge `u -> v`, define

```text
d(u,v) = phase(v)-phase(u) mod 12.
```

The oriented one-edge targets are exactly

```text
(Z/12) \ {4,8}.
```

So phases `4` and `8` have unweighted response distance two; every other
phase has distance one.  The weak-tree edges alone carry

```text
{1,2,3,5,6,7,9,10,11},
```

whereas the strong-cotree edges alone carry

```text
{0,1,3,5,6,7,9,11}.
```

Neither edge set is phase-complete.

For a nonempty oriented **vertex-simple** path `P`, use the same
multiplicative objective as `THM-3269`:

```text
K(P) = product_(e in P) kappa_e.
```

There are exactly 21,226 oriented nonempty vertex-simple paths. The minimum
for each phase is unique, with the unavoidable reversal pair for the
self-inverse phases zero and six:

| phase | minimizing oriented path | exact product |
|---:|---|---|
| 0 | `3-11-10` and reverse | `25025006903228968482171/9358007889540925731463` |
| 1 | `7-22-21` | `4117478020950840745301900523488000/181067356219656228631993554957801` |
| 2 | `21-2` | `2556604002388690672/627821797117010301` |
| 3 | `18-19` | `4566023170631/1510274189347` |
| 4 | `19-2-21` | `5283860475339823645447139037140308586/736904568303910466989977698132844195` |
| 5 | `21-22` | `325385963940378812800/15279959903050007919` |
| 6 | `7-22` and reverse | `12654135326210/11849988963879` |
| 7 | `22-21` | same as phase 5 |
| 8 | `21-2-19` | same as phase 4 |
| 9 | `19-18` | same as phase 3 |
| 10 | `2-21` | same as phase 2 |
| 11 | `21-22-7` | same as phase 1 |

For phase zero, the empty path is excluded. Arbitrary walks are deliberately
not admitted: the backtrack `7-22-7` has phase zero and cost

```text
(12654135326210/11849988963879)^2,
```

strictly below the listed simple-path minimum. This is the canonical hostile
which forced the repaired object definition. Two non-obvious weighted detours
remain: phases zero and `+/-1` prefer weak-tree detours to an available direct
edge. Phases `+/-4` require two hops even before weighting.

Most strikingly, every one of these global phase-target geodesics lies in
`T_*`.  Their seven-edge union is the forest

```text
3--11--10

18--19--2--21--22--7

and isolated vertices 13,16,17.
```

Call this union the **phase-geodesic backbone**.  The strong cotree carries
none of the twelve global phase minima: the weak/strong polarization is
`12/12` versus `0/12` on this target.

This is not an ordinary shortest-path-tree statement.  The weak tree carries
only 48 of the 132 endpoint-specific weighted geodesics and is merely the
60th-lightest of all 74,748 spanning trees.  It becomes globally geodesic
only after endpoint pairs are quotiented by their preserved critical phase.
That quotient, rather than a hidden metric-tree theorem, is the mechanism.

## Candidate 2: a target-preserving symmetric-exchange subatlas

A tree carries all twelve global phase geodesics exactly when it contains the
seven-edge backbone.  Exact exhaustion gives

```text
205 spanning trees contain the backbone,
36 ordered complementary-tree charts contain it on their first side.
```

Join two of the 36 charts by a symmetric exchange that retains the backbone.
The induced exchange graph has

```text
36 vertices, 120 edges,
minimum/maximum degree 4/10,
one connected component,
diameter 5.
```

The weighted chart `T_*` has eight such one-exchange neighbors.  Thus phase
target preservation does not freeze the whole polarization: it leaves a
connected 36-chart subatlas.  The exact clutch product then selects `T_*`
uniquely inside that subatlas.

This gives symmetric exchange a genuine preserved predicate, absent from the
unqualified exchange graph of `THM-3260`: retain one globally minimizing path
for every critical phase target.

## Map contract and loss ledger

The connection is

```text
source: nonempty oriented vertex-simple paths in the weighted core
target: J_12, canonically coordinated by [16-17] -> 1
map:    P=(v_0,...,v_m) |-> sum d(v_(i-1),v_i)
       = [v_m-v_0]
preserved predicate: prescribed critical phase and exact minimum clutch product
needed sidecars: THM-3254 response bank, THM-3269 weights/root,
                 THM-3273 critical quotient
cheapest decisive test: exhaustive exact simple-path and symmetric-exchange scan
```

The destroyed information is load-bearing:

| forgotten coordinate | exact evidence | consequence |
|---|---|---|
| endpoint identities | only `48/132` pair geodesics lie in `T_*` | phase geodesicity cannot be lifted to all row pairs |
| intermediate path and edge multiset | many paths telescope to one divisor difference | the phase label is not a path reconstruction |
| `6229`-primary critical coordinate | target is only the Hall-`{2,3}` quotient | no full-Jacobian or vertex-separation claim |
| six vertex collisions | phase fibres are `1,1,1,3,3,3` | no `V(G_0) <-> C_12` bijection |
| maximizing trap states and interval endpoints | `K(P)` retains only products of maxima | path cost is a ranking functional, not a composition law for traps |
| response signs and positive cones | graph/path data use only exact scalar strengths | no positive current or Gaussian moment follows |
| physical owner/Singer phase | no map from rows to `THM-3255` words was constructed | the internal phase is not an LRC phase observable |

In particular, saying that two clutch overlaps multiply along a graph path
does **not** prove that the corresponding two state traps can be realized
simultaneously.  The product is deliberately the chart-selection functional
inherited from `THM-3269`, nothing stronger.

## What this changes and the next hostile tests

The useful reframe is that an internal phase need not enumerate vertices to
organize relative response operations. The weighted graph has a canonical
root and full-critical generator from THM-3269, plus a canonical minimum
vertex-simple response path for every relative phase. THM-3269 now supplies
the abstract cyclic vertex ordering needed by the augmentation bridge; what
remains missing for LRC is a physical owner/ancestry intertwiner.

The next exact tests should be:

1. attach a literal same-ancestry `(q,R)` record to an oriented backbone path
   and test whether its Singer-norm phase equals the critical path target;
2. retain the full `C_74748` divisor coordinate and determine exactly which of
   the six mod-12 collisions it separates, without pretending the
   `6229`-primary residue is physical;
3. replace the product objective by bottleneck and lexicographic endpoint-bank
   objectives and test whether the same seven-edge backbone survives;
4. classify the contracted five-component multigraph obtained from the
   backbone and explain the 36-chart connected subatlas structurally rather
   than only by exhaustion.

## Exact companion

Run

```text
python3 04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py
python3 -O 04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/gmc_weighted_critical_phase_path_atlas_scout_20260803.out
```

The current hashes are

```text
script  0beca08f9214bd6befeafdc0ccc7648be33dd8c1c2d2f12560d2e56708f2cfcf
output  a55cf2bb5c130d06aaae1c84dc70d2e78e570d82253c8fdb624c4f9d86ede2e0
```

Both normal and optimized output agree byte-for-byte with the stored output.
The companion uses exact integers and rationals only, has no randomness or
floating point, and raises on every failed certificate.

The repaired package is promoted as THM-3277.
