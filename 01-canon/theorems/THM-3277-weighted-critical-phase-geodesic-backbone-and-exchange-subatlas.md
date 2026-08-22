---
id: THM-3277
title: "Weighted critical-phase geodesic backbone and exchange subatlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the canonically
  weighted THM-3269 response core, every target in the genuine critical C12
  quotient has an exact minimum nonempty oriented vertex-simple path. All
  twelve target minima lie in the selected weak tree and their union is a
  seven-edge forest. Exactly 205 spanning trees and 36 ordered bispanning
  charts contain this backbone; the induced target-preserving symmetric-
  exchange graph is connected with diameter five, and the exact clutch weight
  again selects one chart. The theorem is about simple paths and an internal
  quotient, not arbitrary walks, response composition, an owner phase, or LRC.
source: root/2026-08-03
audit: >
  The exact companion independently rebuilds the response rows, trap
  strengths, all 74,748 spanning trees, 9,920 ordered bispanning charts, full
  critical coordinate and C12 phase chart; enumerates all 21,226 oriented
  nonempty vertex-simple paths; and exhausts the backbone carriers and
  symmetric exchanges. The first hostile audit reproduced every finite
  census but rejected an arbitrary-walk reduction and caught a generator-
  gauge mismatch. The package was repaired to the vertex-simple path object
  and THM-3269 full-generator gauge, with the cheap backtrack retained as a
  hostile. A fresh independent re-audit matched every path, product, tie,
  carrier and scope statement. Normal, optimized and stored outputs agree.
depends_on:
  - THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary
  - THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization
  - THM-3273-critical-group-c12-quotient-and-relative-c7-equivariance-boundary
related:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-3274-norm-fibre-constrained-phase-transfer-and-refinement-invoice
script: 04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py
output: 05-knowledge/results/gmc_weighted_critical_phase_path_atlas_scout_20260803.out
script_sha256: 0beca08f9214bd6befeafdc0ccc7648be33dd8c1c2d2f12560d2e56708f2cfcf
output_sha256: a55cf2bb5c130d06aaae1c84dc70d2e78e570d82253c8fdb624c4f9d86ede2e0
hash_basis: LF-normalized bytes
---

# THM-3277 -- weighted critical-phase geodesic backbone and exchange subatlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Canonical phase gauge and admissible paths

Let `G_0` be the 12-vertex, 22-edge response core of
[THM-3260](THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary.md).
[THM-3269](THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization.md)
assigns every edge `e` a positive rational clutch strength `kappa_e>1` and
selects the unique minimum ordered complementary-tree chart

```text
(T_*,T_*^c),             root r=17.                     (1)
```

Its unique incident primitive full-critical direction fixes the genuine
order-twelve gauge

```text
g=[16-17] -> 1 in J_12.                                 (2)
```

For reference, the secondary singleton generator class is
`[7-17]=7g`. It is not the primary gauge.

The rooted vertex phases in gauge `(2)` are

```text
row:     2   3   7  10  11  13  16  17  18  19  21  22
phase:  10  10   7  10   4   4   1   0   1   4   8   1. (3)
```

They have fibres of sizes `1,3,3,1,1,3`; the vertices are not a `C_12`
torsor.

An admissible object in this theorem is a **nonempty oriented vertex-simple
path**

```text
P=(v_0,...,v_m),          v_i all distinct.             (4)
```

Its target and weight are

```text
d(P)=phase(v_m)-phase(v_0) mod 12,
K(P)=product_(e in P) kappa_e.                          (5)
```

Arbitrary walks are not included. This restriction is necessary: the
backtrack `7-22-7` has target zero and weight

```text
(12654135326210/11849988963879)^2,                      (6)
```

strictly below the minimum simple-path weight for target zero. Thus no
walk-to-simple-path reduction is claimed.

## 2. Complete weighted phase-path atlas

The core has exactly 21,226 nonempty oriented vertex-simple paths. The one-
edge targets are all residues except 4 and 8. Exact exhaustion gives the
following global minima:

| target | minimizing oriented path(s) | exact product |
|---:|---|---|
| 0 | `3-11-10`, `10-11-3` | `25025006903228968482171/9358007889540925731463` |
| 1 | `7-22-21` | `4117478020950840745301900523488000/181067356219656228631993554957801` |
| 2 | `21-2` | `2556604002388690672/627821797117010301` |
| 3 | `18-19` | `4566023170631/1510274189347` |
| 4 | `19-2-21` | `5283860475339823645447139037140308586/736904568303910466989977698132844195` |
| 5 | `21-22` | `325385963940378812800/15279959903050007919` |
| 6 | `7-22`, `22-7` | `12654135326210/11849988963879` |
| 7 | `22-21` | same as target 5 |
| 8 | `21-2-19` | same as target 4 |
| 9 | `19-18` | same as target 3 |
| 10 | `2-21` | same as target 2 |
| 11 | `21-22-7` | same as target 1 |

The only ties are the unavoidable reversal pairs at the self-inverse targets
zero and six. The empty path is excluded at target zero.

Every one of these twelve target minima lies in `T_*`. Their seven-edge union
is

```text
E_back={
 (2,19),(2,21),(3,11),(7,22),(10,11),(18,19),(21,22)
}.                                                       (7)
```

It is the forest

```text
3--11--10

18--19--2--21--22--7

and isolated vertices 13,16,17.                         (8)
```

Call `(7)` the **critical-phase geodesic backbone**. The weak/strong
polarization is exact on this quotient target:

```text
T_* carries 12/12 global phase minima,
T_*^c carries 0/12.                                     (9)
```

## 3. Why the quotient matters

The backbone property is not an ordinary all-pairs shortest-path-tree
theorem. Among the 132 ordered pairs of distinct vertices, `T_*` contains
only 48 endpoint-specific minimum product paths. Among all 74,748 spanning
trees, it has total clutch-product rank 60 rather than rank one.

Thus the exact mechanism is

```text
endpoint pairs -> critical phase difference -> one minimum per target.   (10)
```

The quotient in `(10)` destroys endpoint identities and allows a tree which
is not globally all-pairs geodesic to carry every target-class minimum.

## 4. Target-preserving symmetric-exchange subatlas

A spanning tree carries all twelve phase minima if and only if it contains
`E_back`. Exact exhaustion, independently checked by a contracted matrix-tree
determinant, gives

```text
205 spanning trees contain E_back,
36 ordered complementary-tree charts contain E_back on their first side. (11)
```

Join two of the 36 charts when one symmetric exchange transforms the first
into the second while retaining the backbone. The induced exchange graph has

```text
vertices=36, edges=120,
degree range=4,...,10,
one connected component,
diameter=5.                                              (12)
```

The selected chart `T_*` has degree eight in `(12)`. The phase predicate does
not freeze the chart; it leaves a connected 36-chart subatlas. The exact
clutch product then selects `T_*` uniquely within that subatlas.

This gives symmetric exchange a genuine preserved predicate:

```text
retain one globally minimizing vertex-simple path for every critical phase.
                                                                  (13)
```

## 5. Map contract and loss ledger

The typed connection is

```text
source: nonempty oriented vertex-simple paths in the weighted core;
target: J_12, coordinated by [16-17] -> 1;
map:    P |-> sum edge increments = [v_m-v_0];
preserved: target phase and exact minimum clutch product;
sidecars: THM-3269 weights/root/generator and THM-3273 quotient.  (14)
```

It forgets endpoint identities, intermediate vertices, the edge multiset,
the `6229`-primary critical coordinate, maximizing trap witnesses, response
signs and every physical owner label. Six phase-fibre collisions remain.

The product `K(P)` is a ranking functional. It does not say that the trap
witnesses on consecutive edges can occur simultaneously, and it is not a
response-composition law. The theorem supplies no vertex torsor, positive
current, Gaussian-moment functional, Singer owner phase, physical ancestry
map or `LRC(14)` decrement.

## 6. Reproduction

Run

```text
python3 04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py
python3 -O 04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py
```

and compare LF-normalized bytes with the declared output. The companion uses
exact integer and rational arithmetic only, has no assertion node, floating
literal, randomness or graph-library dependency, and freezes the complete
path-bank digest.

QED.
