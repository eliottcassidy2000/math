# Exact ports across growing layers, and a toric microcosm

2026-10-25. **PROVED scoped results, independently reviewed, with
FINITE-EXACT controls.** Unrestricted A000568 compression, full black-core
connectivity, Collatz, and RH remain **OPEN**. No equivalence between those
problems and P versus NP is established.

The requested endpoint-three layer and the first unit-gap layer had already
advanced on origin/main. We inherited those results before selecting the
new targets. The current outcome is:

| Lane | Result | Scope |
|---|---|---|
|[Endpoint detours](endpoint.md)|Two complementary three-black ledgers close the closed-middle endpoint-four family with incoming ports I={x,y,c}; a further ledger moves a growing prefix|Every n>=11 for the second-port family; explicit chord and n>=2h+8 guards for the prefix|
|[Anchor transitions](anchor.md)|Canonical rank-raising augmentation gives weighted finite-state transfers, with recanonicalization between steps|Every fixed feedback-vertex rank; general automata proved to exist, not materialized|
|[Parent forest](anchor.md#4-the-unrestricted-parent-forest-and-its-finite-columns)|Every class has a canonical parent decreasing order and feedback rank by one; maximum transitive size alpha stays fixed|Tree at TT_l has height at most 2^l-1-l; all tournament classes occur|
|[Zero-heavy layers](gap.md)|Exact class formulas including P3-plus-matching defects and the missing-gap/unit-gap contact|k>=3 and 0<=delta<=2k-1; k=2 through delta3 separately FINITE-EXACT|
|[F1-facing construction](f1_bridge.md)|Nested finite sets realize commuting toric powers and generalized Pythagorean orbit counts; the exceptional set grows exactly|Explicit arithmetic/counting map; blue-path predicate preserved only in the stated inherited subfamily|

## What the smaller state must retain

The regular pentagon disproves the shortcut that an ordinary canonical
child record automatically projects to an ordinary canonical parent.
The repair is a canonical **augmentation** owner, followed by exact
recanonicalization before the next step. Singleton disagreement sets tell
exactly which clone pairs appear after deletion.

The first shared-endpoint mass layer now has the formula

    Z_(k,k+1)=(2k-1)!! [binom(3k+1,k+1)-3k-2].

At k=3 its all/strong counts are 2985/650; at the next deficit they are
6060/1046. The new four-unit contact has a decisive hostile: its internal
strong T4 alone cannot say which middle vertex was the endpoint. The
retained external matching partner distinguishes the roles. The quotient
keeps that partner, so the formula does not silently merge them.

There is also a rigorous limit on the desired scalar recurrence:
A000568 grows at least as 2^(n choose2)/n!, and therefore is not
P-recursive. A fixed-order linear recurrence with polynomial coefficients
cannot be the unrestricted answer. Nonlinear and growing-state recurrences
remain available; the parent forest gives one exact structural framework.

The complete alpha=3 forest is a small example: nine classes, with a
Hamilton-path-count spine 1 -> 3 -> 13 -> 45 -> 189 ending at Paley7.
Those numbers label actual class vertices; they do not determine the
parent relation or prove a black detour.

## A positive arithmetic extension and its exact boundary

Take nested finite sets B subset A, of sizes m>n. Pair the two copies of
B in X=A disjoint-union B, leaving d=m-n fixed states; q=m+n is the total
size. On X squared the diagonal involution has d^2 fixed points. The counts

    qd, (q^2-d^2)/2, (q^2+d^2)/2

give the Euclid triple. Taking (A^r,B^r) realizes m,n -> m^r,n^r and
commutes under multiplication of exponents. This is a literal bridge to
the toric conic UV=W^2 with its prime-indexed Frobenius lifts.

The old microcosm has d=1. Powering makes
d_r=(m-n) sum_j m^(r-1-j)n^j, so the exceptional set cannot be discarded.
For example the third power of (3,4,5) gives (63,16,65). It does not give
the blue-count triple at q=63. The arithmetic extension works precisely
because it retains a coordinate the original special family fixed.

The RH-facing note separates this explicit construction from the open
global positivity problem. Its proposed interface is the correctly
normalized Weil form on all compact-support horizons, including the
archimedean term and moment constraints. It supplies a classification of
proof duties and tests for genuine reductions, not a claim that famous
open problems have already been shown equivalent.

## Incoming themes and next targets

The prior [owner theorem](../zero_heavy_ports_20261020g/zero_heavy.md),
[missing-gap theorem](../zero_heavy_next_20261020h/README.md), and
[endpoint-four frontier](../endpoint_four_frontier_20261021/README.md)
are the closest mechanisms. Recent Session58 work contributes the
[support/divisibility map](../session58_hourglass_complement_geometry_20261025/projective-support-bridge/PROJECTIVE_SUPPORT_BRIDGE.md)
and [absolute-geometry audit](../session58_hourglass_complement_geometry_20261025/absolute-geometry/ABSOLUTE_GEOMETRY_BRIDGES.md).
All preserve specified predicates while exposing a lost arithmetic or
incidence coordinate. We use the META-PATTERNS cards “Search the statement
before the method,” “Retain local incidence through counting,” and “Type
every analogy and every implication”; no new universal card is promoted.

The board is endpoint ports, canonical owner languages, integer gap roles,
finite alpha columns, toric power maps, and global correlations. Anchor:
detours and counting. Niche: exact overlap/missing-gap ownership. Wildcard:
the growing exceptional-set interpretation of F1-adjacent arithmetic.

The next precise tasks are the residual endpoint-four endcaps, a compact
realization of the rank-two augmentation machine, and deficit delta=2k
where a second missing gap becomes possible. General connectivity and
uniform owner compression are not inferred from the completed layers.

See [AUDIT.md](AUDIT.md) for independent reviews, correction lineage,
computational universes and reproducibility. Run `python
05-knowledge/results/ports_owners_f1_20261025/verify.py --replay` from the
repository root. The two `endpoint_probe.py` / `endpoint_syndrome.py`
scripts preserve exploration only; the five named checkers supply the
required receipts.
