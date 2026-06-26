---
id: HYP-3062
title: Roth-Minkowski Diophantine lattice fence
status: SYNTHESIS / proof-interface sidecar; not a proof
source: codex-2026-06-26-S226
tangent: T1144
related:
  - HYP-3061
  - HYP-3058
  - HYP-3054
  - HYP-3009
  - HYP-3008
  - HYP-2998
  - HYP-2982
  - HYP-2963
  - HYP-2764
  - HYP-2614
  - HYP-2613
  - HYP-2612
  - HYP-2608
  - THM-538
  - OPEN-Q-108
---

# HYP-3062: Roth-Minkowski Diophantine Lattice Fence

## Claim

Roth's theorem and Minkowski's theorem should enter the LRC14 proof surface as a
paired controlled-forgetting sidecar, not as two scalar estimates.

Minkowski supplies an existence/volume gate: once the relation lattice, ambient
rank, covolume, and convex body are named, a short nonzero lattice relation may
be forced.  Roth supplies the algebraic near-miss fence: once the algebraic
target, field degree, height scale, approximation exponent, and finite
exception list are named, an infinite family of too-good rational approximants
cannot remain hidden behind a quotient.

The proof carrier is therefore not "volume is large" or "approximation exponent
is too good."  It is the combined ledger:

```text
roth_minkowski_sidecar:
  relation_lattice
  ambient_rank
  covolume
  convex_body
  successive_minima_profile
  short_vector_certificate
  algebraic_target
  field_degree
  height_bound
  approximation_exponent
  epsilon_margin
  exceptional_approximants
  low_height_wall_ledger
  deleted_anti_cosets
  residue_signed_tail
  proof_exit
```

## Controlled-Forgetting Rule

A quotient using geometry-of-numbers or Diophantine-approximation pressure is
legal only after the quotient retains, reconstructs, annihilates, descends, or
explicitly names debt for the following coordinates:

```text
lattice side:
  relation lattice
  ambient rank
  covolume
  convex body
  successive minima
  short-vector certificate

height side:
  algebraic target
  field degree
  height bound
  approximation exponent
  epsilon margin
  finite exceptional approximants

LRC packet side:
  exact M/Farey address
  low-height wall class
  deleted anti-coset
  residue signed tail
  route/status handoff
```

If any of these fields is forgotten without a legal exit, the estimate remains a
scout, not a theorem carrier.

## LRC14 Transfer

The old support-six Minkowski frontier in HYP-2612, HYP-2613, and HYP-2614
already asks for a geometry-of-numbers finish: delete the dangerous wall
classes, control the relation lattice tail, and certify the remaining packet.
HYP-3062 adds the missing Diophantine-approximation guardrail.  If the tail
argument produces rational approximants to an algebraic target, Roth says that
too-good approximants are finite after the target and height scale are fixed.
The finite set cannot be forgotten; it becomes `exceptional_approximants` or
`low_height_wall_ledger`.

This gives a three-stage proof carrier for support-six and related packets:

1. Finite wall deletion: classify low-height walls, anti-cosets, and residue
   signed tails before passing to a tail estimate.
2. Minkowski tail: name the relation lattice, covolume, convex body, and
   successive minima profile; force or exclude short relations with an explicit
   certificate.
3. Roth fence: for algebraic near misses, name the algebraic target, height
   bound, approximation exponent, epsilon margin, and exceptional approximants;
   prove that any infinite leakage would violate the Roth fence or descend to a
   named exact packet family.

The same sidecar links to:

- HYP-2608a wide-spread bound: use the lattice/height fields to distinguish
  genuine spread from quotient-blind approximation.
- HYP-2764 zonotope and geometry-of-numbers carriers: attach covolume,
  convex-body, and successive-minima payloads before projecting.
- HYP-3061 geometry-regime audit: use `geometry_regime_signature` for axis and
  curvature language, then use HYP-3062 when the axis becomes a lattice/height
  estimate.
- HYP-3058 hyperbolic reciprocal packets: do not treat hyperbolic margin as an
  exponent proof unless the algebraic target and height sidecar are named.
- HYP-3008/HYP-3009 automata and Fermat-Catalan packets: power-lift guards and
  automatic language shadows need `height_bound` and
  `exceptional_approximants` before claiming infinite-family exclusion.
- HYP-2998/HYP-3003 Farey/additive-basis lanes: `p/q`, `p+q`, `p*q`, and
  power-stress lanes can feed approximation data only with exact Farey address
  and low-height wall class retained.
- HYP-2982 analytic sieve lanes: large-sieve or smoothing estimates should
  state whether their exceptional set is analytic, lattice, algebraic-height,
  or packet-topological.
- HYP-2963 packet ledger and OPEN-Q-108: add the Roth-Minkowski fields to any
  row that invokes "Minkowski count", Diophantine approximation, height,
  algebraicity, or low-height exceptional walls.

## Tournament Analysis

Vertices are proof carriers and sidecar columns, not runners:

```text
labelled_packet_sheaf
low_height_wall_ledger
relation_lattice_covolume
minkowski_successive_minima_gate
roth_algebraic_height_fence
residue_signed_tail
hyperbolic_reciprocal_signature
automatic_gap_language
raw_volume_or_exponent_scalar
```

Pairwise observable:

```text
retained exact M
relation lattice
covolume
successive minima
algebraic target
field degree and height
exceptional approximants
residue signed tail
route/status handoff
```

Gauge: orient an edge toward the carrier retaining more proof-critical payload,
with lower unnamed debt as the tiebreaker.  A Hamiltonian retention path is:

```text
labelled_packet_sheaf >
low_height_wall_ledger >
relation_lattice_covolume >
minkowski_successive_minima_gate >
roth_algebraic_height_fence >
residue_signed_tail >
hyperbolic_reciprocal_signature >
automatic_gap_language >
raw_volume_or_exponent_scalar
```

The challenged assumption is that the lattice vertex must be a runner or arc.
For this sidecar the plausible vertices are relation lattices, low-height wall
classes, fixed circle sections, residue tails, convex-body facets, algebraic
targets, exceptional approximants, and proof obligations.  The preserved LRC
predicate is route/status after exact packet labels; the destroyed coordinate
is the finite low-height and algebraic-exception structure unless the sidecar is
kept.

## Proof-Task Interface

Next pull:

1. Annotate HYP-2963 and HYP-2614 support-six residual packets with
   `relation_lattice`, `covolume`, `successive_minima_profile`,
   `convex_body_id`, `algebraic_target`, `height_bound`,
   `approximation_exponent`, `exceptional_approximants`,
   `low_height_wall_class`, `deleted_anti_cosets`, `residue_signed_tail`, and
   `diophantine_exit`.
2. Prove the finite low-height wall deletion before applying Minkowski.
3. Use Minkowski only after the lattice and convex-body data are named.
4. Use Roth only after the algebraic target, height, epsilon margin, and
   exceptional approximants are named.
5. Treat any unlisted low-height wall or algebraic near miss as named residual
   debt, not as proof completion.
