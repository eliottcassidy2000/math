---
id: HYP-3038
title: LRC14 q=23 drop/add Haar square repair
status: EVIDENCE / local residual-pair grid and zipper target; not a proof
source: codex-2026-06-26-S201
tangent: T1119
script: 04-computation/lrc14_q23_drop_add_haar_square_codex_s201.py
result: 05-knowledge/results/lrc14_q23_drop_add_haar_square_codex_s201.out
related:
  - HYP-3037
  - HYP-3036
  - HYP-3032
  - HYP-3031
  - HYP-3027
  - HYP-3026
  - HYP-3024
  - HYP-3023
  - HYP-2991
  - HYP-2989
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3038: LRC14 q=23 Drop/Add Haar Square

## Claim

The HYP-3032 squarefree `q=23` petal/covering residual is not an analytic
scalar obstruction.  It is the diagonal of a two-coordinate drop/add square.
The first repair tooth is the mixed Haar/fixed-margin coordinate that keeps
the diagonal doubling match; the second repair tooth is endpoint-owner strip
data that separates petal from covering inside the `q=23` diagonal.

In short:

```text
raw analytic q=23 shadow
-> exact_M_zeta_grid / diagonal-doubling match
-> endpoint-owner strip
-> labelled petal or covering route
```

## Computation

Script:

```text
04-computation/lrc14_q23_drop_add_haar_square_codex_s201.py
```

Stored output:

```text
05-knowledge/results/lrc14_q23_drop_add_haar_square_codex_s201.out
```

The square uses AP13 as base.  Rows are `(base - drop_pair) union add_pair`.

```text
AA: drop(10,13), add(20,26)  PETAL-DIAGONAL
AB: drop(10,13), add(16,24)  CROSS-Q-WITNESS
BA: drop(8,12),  add(20,26)  CROSS-Q-WITNESS
BB: drop(8,12),  add(16,24)  COVER-DIAGONAL
```

The opposite diagonal corners are the HYP-3032 residual pair.

## Readout

The diagonal-doubling corners have the same exact maximin scale:

```text
AA petal: M=2/23 at 7/23, safe=59/6370
BB cover: M=2/23 at 3/23, safe=263/25872
```

The cross corners open into easier direct witnesses:

```text
AB cross: M=1/10 at 1/10
BA cross: M=1/8 at 1/8
```

The local Haar/fixed-margin zeta is nonzero for exact `M`:

```text
zeta(M) = 2/23 - 1/10 - 1/8 + 2/23 = -47/920
```

Other nonzero mixed coordinates include safe measure, bar count, longest bar,
midpoint slack, boundary count, and zero-sum active-pair count.  Magnitude
height and magnitude delta have zeta `0`, so raw magnitude misses this local
interaction.

## Quotient Lesson

Local quotient stress:

```text
status_only                         mixes all four open rows
diagonal_doubling_match             mixes AA with BB
exact_M                             mixes AA with BB
exact_M_plus_safe_body              splits all rows
exact_M_plus_endpoint_owner_strip   splits all rows
diagonal_plus_endpoint_owner_strip  splits all rows
```

The endpoint-owner strip distinction is visible even though both diagonal rows
have `B18Z6` in the coarse endpoint-current word:

```text
AA external owners: 12:26x6,6:20x4
BB external owners: 2:16x6
```

Thus `B18Z6` is a lossy endpoint scalar; the owner strip must retain which
external speed/residue owns the boundary facets.

## Tournament Analysis

Vertices are local proof teeth, not runners:

```text
raw_analytic_q23_shadow
drop_add_row_column_shadow
diagonal_doubling_match
exact_M_zeta_grid
safe_component_body
endpoint_owner_strip
labelled_packet_route
```

Pairwise observable:

```text
status_retention,
route_retention,
zeta_retention,
endpoint_owner_retention,
topology_retention,
family_transfer,
low_cost
```

Switch/gauge:

Orient by majority of retained proof coordinates, with tie Hamiltonian path

```text
labelled_packet_route
> endpoint_owner_strip
> safe_component_body
> exact_M_zeta_grid
> diagonal_doubling_match
> drop_add_row_column_shadow
> raw_analytic_q23_shadow
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

## Proof Target

Candidate local lemma:

```text
For any double-pair drop/add square, diagonal doubling either
  (a) opens the off-diagonal corners to direct q-witnesses,
  (b) descends through a family q-diagonal,
  (c) or exposes endpoint-owner strip data that routes petal, covering,
      K33, or F7/THM-572 debt.
```

For the tested `q=23` square, the observed repair class is:

```text
nested_refinement_to_q23_diagonal_then_owner_strip
```

## Assumption Challenge

Alternate vertices considered: runners, dropped pairs, added pairs, gaps,
endpoints, fixed sections, residues, Fourier modes, Haar tiles, cover arcs,
matroid circuits, and proof obligations.

The chosen vertices are proof teeth because the preserved LRC predicate is
open/boundary status plus theorem-route schedulability after quotienting.
This destroys raw runner names, row/column drop/add margins, exact `q=23`
without owner identity, and analytic clocks unless the mixed coordinate and
endpoint-owner strip are retained.
