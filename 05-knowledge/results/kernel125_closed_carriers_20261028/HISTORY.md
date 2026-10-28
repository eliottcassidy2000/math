# What changed since the first enumeration work

**Status: CITED repository history / PROVED elementary repairs /
FINITE-EXACT replay / OPEN extensions.** Historical dates identify provenance,
not truth. The code comparison starts from commit
b8672af00bcc6d440a9a00a6915325faf39fcdff.

## March prototypes: useful direction, two invalid producers

The early scripts already used the right broad strategy: reduce labelled
tournaments to permutation cycle types, sum pair-orbit fixed counts, and
avoid repeatedly scanning all labelled tournaments.
[a000568_speedup.py](../../../04-computation/a000568_speedup.py), headed
March7 S46f, explored direct sums, sequential processing, modular work and
gcd coordinates. The March16 S90dm
[v3 prototype](../../../04-computation/a000568_v3_s90dm.py) cached powers
and gcd values. The March17
[turbo prototype](../../../04-computation/a000568_turbo.py) continued
the integral accumulation direction.

Their advertised speedups are not adopted as current benchmarks. Two
producers failed before timing comparisons were meaningful:

| Source at the pinned base | First failure | Mechanism and repair |
|---|---|---|
| speedup Strategy A | a(5)=10 instead of12 | An odd cap was combined with an even remaining size; starting at their minimum restarted the recursion on even parts. Restore odd parity at every recursion. |
| v3 | a(2)=2 instead of1 | All integer partitions were admitted, including forbidden even cycles. Restrict every chosen cycle length to odd values. |

These are errors in particular historical scripts, not a refutation of
the classical formula or later exact engines. The v3 claim of avoiding
n! also contradicted its implementation. Its formula/documentation,
partition-count display and rounded asymptotic comparison were repaired
along with the cycle filter. Strategy A and T_v3 are replayed through30;
the other speedup.py strategies and the whole historical corpus are
not certified by this audit.

The repairs are recorded in [MISTAKES.md](../../../01-canon/MISTAKES.md).
[legacy_audit.py](legacy_audit.py) extracts the original and repaired
functions without executing their top-level benchmarks. A separate
cyclic-action Burnside calculation counts pair orbits by averaging
permutation powers, rather than reusing the gcd exponent. At orders0..8
literal ordered-pair tracing also checks all cycle types, including even
ones, and rejects inconsistent orientations. Full receipts are in
[legacy_repairs.json](legacy_repairs.json).

## April pruning: small in the real metric is not zero modulo a prime

The April5 [S25 reflection](../../../07-reflections/three-computational-threads-s25.md)
was a useful prompt to exploit concentration near the identity, but its
normalization and modular inference were wrong. For one3-cycle and
n-3 fixed vertices, the exponent deficit is2n-4 and the conjugacy-class
size is2*binom(n,3). Its contribution relative to the identity is exactly

    2*binom(n,3)*2^(-(2n-4)).

At n=10 this is15/4096, not binom(n,3)/8. Its asserted persistent factor
1.714 in the successive-count approximation is also contradicted by
exact ratios n*a(n)/(2^(n-1)*a(n-1)); the receipt retains these as exact
fractions at10,20,30. A tiny relative real error may still change every
low residue bit. The old body is preserved under an explicit correction.

The rigorous successor is the exact reverse valuation potential in
[production counting](../lca_ports_reverse_20261025/counting/README.md).
It permits discarding a transfer only after proving its total prefix
and suffix binary valuation exceeds the requested precision. It does
not turn a real asymptotic truncation into an exact modular computation.

## June to October: retain what future cycle lengths can test

[THM-514, the unit-tail core kernel and prime-tail ladder](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
summed singleton tails over equal response cores; its834-versus444793
example at100 concerns outer core states, not total bit complexity.
[THM-516, the phi reframe](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md)
expressed the gcd interaction through divisor coordinates. The later
[future-boundary engine](../tournament_future_boundary_20261020/README.md)
made equality of future responses an exact state criterion.
The [direct contracted DAG](../session64_lca_reverse_valuation_20261026/session64-counting/contracted_ordinary.py)
then removed forced zero-cost skips without conflating active stages.

K6 in [THEORY.md](THEORY.md) reuses that exact structure for three
weight channels. At m=200 the graph has495657 vertices and2219922
edges. Building it once and sharing odd factors replaces6659766
separate factor evaluations by2219922. The observed total for construction
and all three evaluations was20.587 seconds. This is an execution receipt,
not a benchmark against a tuned C/GMP implementation or a larger-order
record. Earlier repository work already reaches large scalar orders.

The interpretation of SC(2m), SC(2m+1) is classical
[A002785](https://oeis.org/A002785). Our contribution is the common
integer-gauge implementation, not a newly discovered sequence identity.

## Class generation is a different cost from scalar evaluation

The immediate predecessor built the full pointed root-switch operator
and recovered ordinary classes from its identical-row blocks. That gave
all191536 order-nine representatives, with15449472 re-root queries
and a large parent-deletion cache. K5 preserves its exact pointed
encoding but asks only whether the chosen root is minimum.

| Order9 method | Re-root tests | Retained full root-switch rows | Parent cache cap |
|---|---:|---:|---:|
| Prior full operator | 15449472 | 1716608 | unbounded in that run |
| Minimum pointed key, no degree prefix | 5986841 | 0 | 65536 |
| Minimum key with degree first | 278461 | 0 | 65536 |
| Neighbor-degree refinement integrated from Session87 | 24005 | 0 | 65536 |

The degree-first method rejected1405799 candidates by degree before any deletion
isomorphism. It skipped2049256 higher-degree roots during remaining
candidate scans. It still used the actual parent action on masks, and
retained full accepted representatives. The complete order-nine bank
matches the predecessor by an exact bijective isomorphism audit.
This is a direct improvement to an existing metagraph-derived generator,
not a scalar Burnside sum disguised as an enumeration.

All-size reconstruction follows the proof by induction; execution here
stops at9. The expensive order-ten census and arbitrary native line
transport are not newly implemented. A concrete next optimization is to
enumerate only masks satisfying the degree filter, before constructing
extensions, but it must preserve one representative per parent-mask orbit.
This is OPEN implementation work, not a claimed result.

### Incoming work integrated before publication

At origin/main 7a0bab9bfa, the concurrent
[Session87 synthesis](../session87_four_vertices_enumeration_20261028/README.md)
arrived with the same odd-generator repair, a stronger root invariant,
and a faster complete-prefix scalar implementation. We retained both
correction routes when resolving the shared source edit. Its additional
hypergraph repair is preserved and is outside our independent audit.

The neighbor-degree selector was incorporated and tested through9 here.
It eliminates another108479 candidates after the initial degree filter,
uses24005 competitor-key queries and12229 exact parent lookups, and
preserves the entire class bank. Its observed generation time was9.842
seconds, excluding the exact full-bank comparison. Its rejection
certificate can be smaller without quotienting away a single accepted
child type.

Session87's degree-bucket counter improves the existing complete-prefix
baseline and reports controlled benchmarks. Our shared three-channel
experiment has a different output and graph construction. We do not
compare their wall times as if they had the same workload. The recovered
K4 partial-spread result also agrees with K4 above: injectivity on a
promised set can survive a nontrivial linear kernel, but unrestricted
addition can restore the information requirement.

## Cross-topic maps and their decisive boundaries

| Source -> target | Map and preserved predicate | Lost data / sidecar | Cheapest decisive test |
|---|---|---|---|
| Actual C5^3 charts -> F25 | a -> a0+omega*a1+omega^2*a2; equal unmarked-hub response | Common phase; known cross polarity under global converse | Five Hamilton guards restore the common phase. |
| P7 attachment charts -> semilinear F8 maps | u -> alpha^u; intertwines all21 chart actions | Pinning an owner reduces21 to3; unpinned translations cannot be discarded | One nonzero translation moves the pinned vertex. |
| Chart trees -> joint ports | Path products through the LCA reconstruct all pair gains | One overall chart gauge; retain full nonabelian gains when unpinned | All21^3 four-cell configurations, not just Frobenius phases. |
| Square-class arithmetic -> binary observation spaces | Prime-valuation parity and kernel equations preserve square-total tests for the declared bank | Actual magnitudes, repeated-prime valuations and other arithmetic operations | Adjoin a coordinate probe; K7 computes which formerly invisible directions become visible. |
| Fixed graceful or square-sum equations -> affine path gains | x_v=s*x_u+w, with s=+1 for chosen signed differences, s=-1 for prescribed sums | Root interval, all-different constraints, chosen signs/lengths/squares | The inherited root-feasibility test plus an interior collision; endpoint data alone fails. |
| Full tournament data -> minimum-root certificate | Complete pointed keys give exact acceptance and automorphism order | Other re-root answers are unnecessary after a rejection witness, but the accepted graph remains | Exact class bijection and the81/27 hostile. |
| Odd-cycle DAG -> related scalar sequences | Reweight J to beta*J+(beta-1)*v2(k!)+gamma*k | Individual classes and attachment alignments are not outputs of this scalar map | Independent partition sums through30 and native SC census through8. |

The affine-label row uses the already proved
[previous T4](../relation_tower_frobenius_ports_20261028/THEORY.md).
It does not choose a graceful edge-length permutation or prove a
square-sum Hamilton path exists. The square-class row preserves a
specified parity predicate; Liouville parity can be linear there, while
the Mobius zero/nonzero condition additionally needs squarefreeness.
No quotient is allowed to silently discard that promise.

This supplies a precise relation among five, seven and fifteen here:
three independent five-rotations produce125 charts, their relative
component is a two-dimensional F5 space, and Paley7 realizes a different
21-element semilinear action. The fifteen coordinates in a Hamming code
and the fifteen vertices of R remain distinct objects. K7 explains how
to search for their observation quotients without identifying them by
their cardinalities.

## The 81/27 hostile remains active

All three generator receipts reconstruct the12-vertex control. The two
tournaments have the same67 labelled cyclic triples and462915 Hamilton
paths each, but automorphism orders81 and27. The identity vertex order
is a legal Hamilton path only in the first. Thus cyclic triples and
total path count do not replace the chart action or named legality.
The stronger incoming witness also compares transport rather than just
stabilizer orders; see
[Session86 closed carriers](../session86_relations_attachment_carriers_20261028/CLOSED_ATTACHMENT_CARRIERS.md).
