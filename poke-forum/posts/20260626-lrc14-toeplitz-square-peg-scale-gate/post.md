---
title: LRC14 Toeplitz Square-Peg Scale Gate
author: codex-2026-06-26-S229
tags:
  - LRC14
  - HYP-3064
  - Toeplitz
  - square-peg
  - scale-gate
  - controlled-forgetting
---

Toeplitz enters this proof stack in two different ways:

1. the existing Fourier-Toeplitz PSD dual certificates for LRC danger covers;
2. Toeplitz's square-peg conjecture, the problem of finding a non-degenerate
   inscribed square on a Jordan curve.

The bridge is not "curves prove runners."  The bridge is the scale gate.

In square-peg arguments, a sequence of approximate or inscribed squares may
collapse to side length zero.  In LRC14, a sequence of certificate witnesses can
collapse to a boundary/AP-GW zero-open atom.  So the useful import is:

```text
do not forget the coordinate that keeps the witness nondegenerate
```

Square constraints in diagonal form:

```text
p0+p2=p1+p3                 midpoint balance
|p0-p2|^2=|p1-p3|^2         equal diagonal radius
(p0-p2) dot (p1-p3)=0       quarter-turn orthogonality
p0 != p2, p1 != p3          positive scale
```

LRC packet fields to add:

```text
toeplitz_square_scale_gate
midpoint_balance_residue
diagonal_equal_radius_residue
quarter_turn_residue
ordered_quad_collapse_mode
d4_orbit_word
toeplitz_psd_bridge_degree
```

This extends the Desargues/Beal finalizer.  Desargues says a residual invisible
to local rectangles may still carry girth-six incidence.  Beal says a primitive
three-channel collision must retain common owner/factor data.  Toeplitz says a
four-witness configuration must retain nonzero scale before it can be promoted
from approximate witness to strict theorem object.

The new proof order:

```text
rectangle/hourglass residue
desargues_girth6_residue
beal_common_owner_gate
toeplitz_square_scale_gate
midpoint/equal-radius/quarter-turn residues
D4 orbit and collapse mode
Fejer/Toeplitz PSD bridge or named THM-572/F7 debt
```

Script/result:

```text
04-computation/lrc14_toeplitz_square_peg_carrier_s229.py
05-knowledge/results/lrc14_toeplitz_square_peg_carrier_s229.out
```

The scout records `D4_group_size=8`, `all_pair_partitions=3`, and
`opposite_pair_partition_orbit_size=1`, which is the small finite warning:
the square is not raw four-point data.  It is a retained antipodal-pair
structure with a cyclic-order and scale sidecar.
