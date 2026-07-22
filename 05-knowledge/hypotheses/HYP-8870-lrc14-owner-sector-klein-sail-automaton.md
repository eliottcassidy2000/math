---
id: HYP-8870
title: "LRC14 owner-sector Klein-sail automaton"
status: >
  OPEN completion program downstream of THM-2055. Replace raw two-anchor
  matrices and round-disk scans by signed hull vertices, rational normal cones,
  and primitive lattice points under the tangent-disk roof. The target is a
  finite Farey/Stern--Brocot automaton whose accepting states carry exact
  pair-sum or owner-labelled Euler certificates. No such uniform automaton is
  proved here.
source: codex-2026-07-21-LRC-normal-fan-sail
related:
  - THM-2052
  - THM-2053
  - THM-2054
  - THM-2055
  - HYP-2108
  - HYP-2647
  - HYP-2896
  - HYP-2986
  - HYP-8841
  - HYP-8846
  - MISTAKE-225
---

# HYP-8870 -- owner sectors and ordinary continued fractions

THM-2055 changes the finite object. For each THM-2052 two-anchor star:

1. compute the saturated column configuration and the signed hull `K_U`;
2. delete nonvertex columns from the determinant sidecar only;
3. split parameter space by the rational normal fan of `K_U`;
4. intersect each owner cone with its single tangent disk and the positivity
   cone of the speed row;
5. enumerate primitive slopes by the ordinary Stern--Brocot tree/Klein sail;
6. label every node by exact pair-sum margin, relative-Fejer resonance packet,
   and HYP-2986 open-tope/boundary-cocircuit state;
7. accept only with a strict phase or owner-labelled Euler endpoint.

The key hoped-for compression is **interval acceptance**: adjacent Farey nodes
with the same active hull owner and the same signed phase-height wall word
should admit one symbolic certificate for the whole mediant interval. A failed
interval must split at an explicit event:

```text
hull-owner tie,
positivity or collision wall,
pair-sum ruler change,
relative-Fejer resonance,
endpoint-owner exchange.
```

All event equations are rational linear or quadratic in the two parameters,
so every fixed star has a finite exact cell structure. The missing uniform
theorem is to bound the number of event words over the entire bounded
two-anchor coefficient atlas without generating every matrix.

## Assumption challenge

The vertices are not runners, arcs, primes, or form classes. They are
**signed hull-owner sectors**. This quotient preserves the determinant gate
and rational parameter address; it destroys non-hull runner constraints, which
must re-enter through the phase-height sidecar before any LRC conclusion.

## Tournament analysis

Proof-carrier vertices:

```text
owner_sector_sail
exact_pair_sum_fan
relative_Fejer_cell
endpoint_owner_cocircuit
raw_tangent_disk_scan
Heegner_form_class
raw_relation_matrix
```

Pairwise observable: `(gate exactness, LRC-predicate retention, symbolic cell
coverage, sidecar debt, cost)`. The switch prefers a carrier only if it keeps
the specified primitive slope and either certifies a whole cell or emits an
exact split event. The tie path starts

```text
owner_sector_sail
> exact_pair_sum_fan
> relative_Fejer_cell
> endpoint_owner_cocircuit
> raw_tangent_disk_scan.
```

`Heegner_form_class` and `raw_relation_matrix` rank last because neither
preserves the determinant owner or the weak phase-height predicate.
