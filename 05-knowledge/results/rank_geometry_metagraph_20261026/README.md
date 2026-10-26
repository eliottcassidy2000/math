# Exact source observers, planar frames and tournament class generation

2026-10-26. Research packet with explicit proofs, hostile controls and
independent audits in [AUDIT.md](AUDIT.md). The root27 global rank, a new
plane-colouring obstruction, black-core connectivity and a smaller
aggregate A000568 recurrence remain **OPEN**.

The constructive result is an unrestricted class-generation recurrence,
implemented and checked through order eight. The rank and geometry lanes
identify exactly which proposed finite reductions lose needed information;
the geometry lane also recovers a realizable alternative frame action.

## Results and scope

| Lane | Established result | Remaining obligation |
| --- | --- | --- |
| [Source observer](rank/README.md) | Exact affine pullbacks characterize finite valuation-observer closure; concrete alternative sources match a long actual root27 prefix | A rank for the one fixed source7 through arbitrarily changing shapes |
| [Planar compatibility](unit_distance/RESULT.md) | Joint loop permutations demand rank3; re-pairing a planar hexagon admits a fixed-point-free residual action; full diamond boundary transfer stays positive | A unit-distance gadget that enforces incompatible boundary relations in actual colourings |
| [A000568](enumeration/README.md) | Canonical deletion selects one rooted extension per class, preserving group actions and recomputing maximal clone blocks | Eliminate the need to enumerate the actual pair-free cores |
| [Information and Frobenius defect](information/RESULT.md) | Correct class masses and unit-flow projection; scalar collision at n5; the new planar reflection has exact tournament orientation defect2 | Retain actual port incidence when reducing the recurrence |

The generated A000568 values at orders0..8 are

    1, 1, 1, 2, 4, 12, 56, 456, 6880.

At order8 the construction tests54,256 rooted extensions and accepts6,880
classes. The weighted-core route independently reconstructs all7,413
classes through order8, including4,602 pair-free cores of order8. This is
an all-order correct algorithm run to a stated finite bound, not an
infinite census or a faster scalar formula than the existing reverse engine.

One especially useful hostile is already order4: complementing the
nonpath arcs of one class can land in **any of the four classes**, depending
on the retained Hamiltonian spine. An operation on a marked presentation
does not automatically descend to an operation on its unmarked class.

## Inheritance and the live concept board

Anchor: the source7 representation of the root27 guarded orbit. Niche:
planar realization of full four-frame relations. Wildcard: old tournament
information ideas turned into exact counting measures and carrier tests.
The board has six concepts: exact source; affine carry; full boundary
relation; geometric realization; automorphism action; clone/port ownership.

The closest proofs are the source-fixed centre switches and finite-bank
limits in [Session66](../session66_root27_aliquot_graph_interfaces_20261026/session66-root/CONNECTIONS.md),
the [previous root27 packet](../root27_machine_bridges_20261026/README.md),
and the signed-curl/path-cut reconstruction in the
[metagraph reconstruction note](../metagraph_reconstruction_20261023.md).
The canonical hostiles are an actual precision refill, an abstract frame
permutation without a common planar metric, and an extension whose newly
merged clone blocks change its ownership. None is bypassed by renaming a
coarse invariant as a rank.

The March codec/unlabeling and June binary-intransitivity reflections were
recovered as provenance, then tested against current canon. This exposed
and repaired [THM-895](../../../01-canon/theorems/THM-895-helmholtz-identity-perron-census.md):
for unit antisymmetric arc flow and divergence d, the least-squares
potential is d/n, not d/(2n). Its residual energy has a transitive baseline;
even a transitive triple has nonzero unit-flow curl. A strict ranking is
encoded by signs of potential differences, which discards their magnitudes.
The normal equations, not merely an assigned potential's energy, are now
checked. The correction lineage is preserved in
[MISTAKES](../../../01-canon/MISTAKES.md).

## Incoming arithmetic: an exact six/seven connection with a different carrier

During this session commit `2daf5e00af` delivered the independently audited
[norm-seven arithmetic/F1 packet](../normseven_arithmetic_f1_20261026/README.md).
Its [marked Laurent-lift theorem](../normseven_arithmetic_f1_20261026/f1/RESULT.md)
provides a precise six/seven distinction. For a matrix C invertible modulo
m, an integral unimodular Q congruent to C exists exactly when det(C) is
±1 modulo m. Every sixth-root character automorphism lifts, while exactly
one third of seventh-root character automorphisms lift. One compensating
root character repairs the determinant obstruction; identity padding does
not. These are classical congruence mechanisms applied to the marked
fibre algebra, not a numerological identification with six colours.

Crucially, the *norm-seven* quotient here uses **quadratic** characters:
its root order is2, not7. Torus coordinate repair also does not imply a
positive affine-monoid automorphism or a whole-surface F1 structure. The
incoming boundary correction changes the complete central divisor's
rational-point possibilities to zero/eight/infinite; the eight singular
rational points can fail to lift to the smooth resolution. This changes
which arithmetic boundary data must be retained in any proposed gluing.

The separately requested norm-seven session completed during integration.
Its [resolved-fibre result](../normseven_descent_moments_20261026/RESOLVED_FIBRES.md)
gives two twists with the same residue field Q(sqrt(13)) and identical
good-prime fibre zeta functions, but minimum resolution-lifting degrees
two and four. The distinction is a tangent-conic obstruction at3. The
[parent review](NORMSEVEN_REVIEW.md) checks its proof and independently
replays the packet. This is about the specified closed fibres, not entire
surface zeta functions or the interior integer-square problem.

Here are the exact parallels and their limits:

- A root-fibre coordinate automorphism needs a determinant sidecar to
  extend to the ambient torus. A frame permutation needs a Gram realization
  to extend to one planar palette. These are different extension problems;
  no reduction between their solution sets is supplied.
- A geometric reflection can act freely on residual frames but has
  minimum defect2 on any attempted strict tournament relation on those
  same six ports. The signed pair-orbit defect is an exact common test,
  not a colouring proof.
- Canonical class generation retains the whole group action, not merely
  its order. Source-observer closure retains affine pullbacks, not merely
  the previous valuation depths. Both make the allowed next operation
  part of the state specification, without equating their dynamics.

## Strongest next proof targets

1. For source7, find an arithmetic bound on the adaptively pulled-back
   centre family which is not a finite-bank statistic or a restatement
   of guarded survival. The hostile families vary the source and therefore
   leave this precise avenue open.
2. Use the realizable re-paired hexagon action to design an actual gadget
   whose **entire** boundary-colouring relation implements the required
   transport. Check all cross-unit contacts and accidental identifications.
   A symmetry of one drawing does not itself impose a colour relation.
3. Seek an aggregate transfer rule on pair-free core types with the actual
   vertex action and clone mergers retained. The catalog now gives exact
   hostile tests and parent/deletion receipts for such a proposed rule.

## Reproduction

From the repository root:

    python 05-knowledge/results/rank_geometry_metagraph_20261026/verify.py
    python agents/check_docs.py

Each lane's normal and optimized execution must have identical parsed
receipts. Finite universes, positive controls and hostile controls are
listed in its note. The theorem arguments are written mathematical proofs
and independent reviews, not Lean certificates.
