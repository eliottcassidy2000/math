# Binary digraph matrices and the tournament boundary

2026-11-10. **PROVED scoped constructions and counting identities;
FINITE-EXACT implementation controls. Tournament spectrum completeness
remains OPEN.** The proposed implication that H=7 or21 forces a loop or a
mirrored00 pair is **REFUTED** under ordinary adjacency. This was already
settled in Sessions234–235; their exact classifications are recovered here.

This session adds a compact strong oriented constructor, a universal
endpoint-contraction bridge from tournaments to arbitrary digraphs, an
explicit loss of shared paths under tournament-resolution isomorphism, and
weighted refinements of the inherited enumeration engine. No external
priority claim is made for these elementary mechanisms.

## What a square binary matrix means

For A an n by n binary matrix, use the ordered Hamiltonian-path total

    H(A)=sum_(p a permutation of the vertices) prod_(i=1,...,n-1) A[p_i,p_(i+1)].

Every vertex occurs exactly once; reverse orders count separately when
valid. Diagonal entries never occur, so adding or removing loops leaves H
unchanged. This is not the number of walks, cycles, or isomorphism classes.

For each unordered pair, ordinary adjacency permits four states:
00 means no arc,10 and01 mean one direction,11 means both directions.
With N=binomial(n,2), the loopless families are:

| Family | Allowed pair states | Labelled matrices | Isomorphism sequence |
|---|---|---:|---|
| Tournament |10,01|2^N|[A000568](https://oeis.org/A000568)|
| Oriented graph |00,10,01|3^N|[A001174](https://oeis.org/A001174)|
| Semicomplete digraph |10,01,11|3^N|A001174|
| General digraph |00,10,01,11|4^N|[A000273](https://oeis.org/A000273)|

A tournament is exactly an oriented graph that is also semicomplete:
`A_ii=0` and `A_ij+A_ji=1` for every distinct i,j. Thus forbidding reciprocal
arcs is not enough; missing pairs must also be forbidden. Allowing all
diagonal bits gives2^(n^2) labelled matrices. Their isomorphism count is
[A000595](https://oeis.org/A000595); graph isomorphism permutes rows and
columns by the same vertex permutation.

## The inherited seven and twenty-one classification

The following ordinary adjacency matrices have zero diagonal, no00 pair,
and exactly7 and21 paths respectively:

    H=7:             H=21:
    0111             01111
    0011             00110
    1001             00011
    1000             10001
                     11000

Their orders4 and5 are minimal even among arbitrary digraphs. The first
matrix's paths are0123,1203,1230,1302,2013,2301,3012. These are inherited
[Session234 DG1–DG5](../session234_digraphs_20261109/PROOFS.md), not discoveries
of this session. A semicomplete H7 graph has exactly two reciprocal pairs;
a semicomplete H21 graph has one, two or three. Loops are irrelevant.

The complete all-order generator is an ordered sum of strongly connected
cores, with arbitrary transitive singleton gaps:

- H7: one of three strong order-four cores.
- H21: one of32 strong cores (16 of order5 and16 of order6), or one H7
  core and one of the two H3 cores, in either order.

Cross-component arcs all point forward, so Hamiltonian counts multiply.
For example the numbers of loopless semicomplete isomorphism classes have
ordinary generating functions

    U7(x)=3x^4/(1-x)^2,
    U21(x)=(16x^5+16x^6)/(1-x)^2+12x^7/(1-x)^3.

The [general fixed-count grammar](../session234_digraphs_20261109/GENERAL_GRAMMAR.md)
proves a finite strong-core alphabet for each semicomplete total. This
does not extend to oriented graphs: Session235 supplies arbitrarily large
strong oriented graphs of any fixed H>=2.

## Interpreting eleven as an empty edge

If both00 and11 mean empty, the effective adjacency is

    C(A)_ij=A_ij(1-A_ji), i!=j; C(A)_ii=0.

Its image is the oriented graphs. If raw00 is forbidden instead, raw11 is
the unique encoding of an effective missing pair and the same oriented
family results. If raw00 is assigned both arcs, the separate decoder is
`J-I-A^T`; see [Session235's three-convention dictionary](../session235_arc_pair_models_20261109/PAIR_MODELS.md).
These rules must be fixed before computing a count.

For example the raw matrix `01110/00111/10011/00001/10000` has ordinary
H11, but cancellation H7. It has no loop and no raw00 pair. The correct
necessary condition in cancellation semantics is an **effective** missing
pair, which can be written raw00 or raw11. The ordinary H7 example above
instead cancels to H1. These outcomes are inherited Session235 results.

Our [fibre refinement](enumeration/NOTE.md) retains the symmetry action on
missing pairs. An oriented graph with m missing pairs has2^m labelled raw
encodings, but its unlabelled preimages are orbits under its automorphism
group. For directed C4, four labelled encodings become three types with
ordinary counts4,4,8 and identical cancelled count4.

## A compact constructive extension

[Construction C1–C3](construction/PROOF.md) proves that every h>=3 has an
explicit strong oriented realization on

    ceil(log2(h-1))+2 vertices.

Write h-2 in binary, use a transitive k-chain with k=bitlength(h-2), and add
`chain -> a -> b`. Add `b -> i` exactly at the1-bits. The highest bit keeps
the graph strong. Paths through a→b compress to the inherited binary
semicomplete construction; exactly one additional path avoids that arc.
There are precisely `k-popcount(h-2)` missing pairs.

Equivalently, start with a tournament of count2^k+1 and delete selected
outgoing arcs of b, keeping its arc to the first chain vertex. Each deletion
subtracts its own binary weight. No mixed terms occur because one path
cannot use two outgoing arcs of b. At order k+2 this family supplies every
integer in `[2^(k-1)+2,2^k+1]`. This is not a minimum-order assertion or a
tournament realization: missing pairs are essential. Oriented-spectrum
completeness was already known in Session235; the logarithmic strong
construction and exact hole budget are the extension here.

## Why general digraphs belong inside tournament calculations

[CT1–CT3](completion/PROOFS.md) retain the ordered entry and exit of each
selected path block. Contracting such blocks computes the number of paths
using all selected arcs. The quotient can have00 or11 even when the source
is a tournament. Conversely, every m-vertex digraph is such an endpoint
quotient of a2m-vertex tournament with a designated matching.

Thus every nonnegative integer occurs as a **constrained** tournament path
response. Literal lifts of the inherited H7 and H21 digraphs have matching
responses7 and21 but unrestricted tournament totals79 and2621. Forgetting
the matching would confuse two different statistics.

Even retaining every resolved tournament isomorphism class is insufficient
to reconstruct a semicomplete total: two strong order-four digraphs both
resolve to the same pair of H5 tournament classes, yet have H6 and H8.
Their resolved path sets intersect in4 and2 paths. The missing information
is the identity of paths shared between resolutions, restored by named
matrices or joint attachment responses.

## Predicting related enumeration sequences

[The weighted Burnside polynomial](enumeration/NOTE.md) counts00,11 and
single-direction pairs separately, refining the inherited scalar palette
theory. Substitutions recover A000568, A001174, A000273 and A000088; a
diagonal variable gives A000595 and its loop-degree refinement A384105.
Tournament specialization excludes every vertex permutation with an even
cycle, because it would reverse some unordered pair. Symmetric00/11 states
permit those symmetries in the larger families.

Marking one ordered Hamiltonian path removes all automorphisms. If x,y,z
weight absent, reciprocal and single pairs, then for n>=1

    sum_[D] H(D)/|Aut(D)| x^m y^b z^s
      = (y+z)^(n-1) (x+y+2z)^(N-n+1).

For cancellation counts the first factor becomes z^(n-1). These exact
identities separate unlabelled class counts from path-weighted counts.
Neither the class count nor its mean path total determines which individual
H-values occur. The tournament conjecture that7 and21 are the only missing
positive odd values remains OPEN.

## Evidence and reproduction

The implementation controls use all small binary words and matrices,
independent literal paths and subset DP, and explicit hostile examples;
they do not repeat the inherited large censuses. See [the recovery board](RECOVERY.md),
[independent construction audit](history/COMPACT_AUDIT.md), and
[independent final review](history/FINAL_AUDIT.md).

Run `python3 -B 05-knowledge/results/digraph_matrix_boundary_20261110/reproduce.py`
from the repo root to verify saved normal/optimized evidence, source hashes,
manifest and local links. Add `--rerun` for fresh executions. The proofs,
rather than extrapolation from these finite checks, supply all-order claims.
