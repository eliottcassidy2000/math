# Independent finite-classification audit

Date: 2026-10-19 (UTC).

**Result: verified.** The finite rank-at-most-two classification in
`THREE_FREE_COLOR_CLASSIFICATION.md` agrees with a separately implemented,
standard-library-only exhaustive enumeration.

The independent checker is `verify_rank2_independent.py`; its complete,
deterministic output is `independent_rank2_receipt.json`. The receipt includes
all 96 disconnected support triples, both of their states, and SHA-256 hashes
of the four original files inspected. No original file was modified.

## Exhaustive counts

- All 63 nonempty simple edge sets on four labeled vertices were ranked.
- Exactly 25 have rank at most two: 6 singletons, 15 pairs, 4 triangles.
- All 25^3 = 15,625 ordered support triples were classified.
- 1,621 triples have no rainbow spanning tree.
- 13,908 triples have a connected, nonempty width-two state graph.
- Exactly 96 triples have a disconnected state graph.
- Each of those 96 graphs has exactly two vertices, both isolated.
- The underlying spanning trees give 24 star cases and 72 path cases.

## Independence and coverage

The checker does not import or execute the original checker or its helper.
It computes graph rank using Boolean vertex-set reachability rather than the
original edge-union partition routine. It identifies all 16 spanning trees
as connected three-edge bit masks. A rainbow state is accepted precisely
when the union mask of its three chosen contacts belongs to that tree set;
this also excludes repeated contacts automatically.

For state components, the checker unions states through shared
(color, contact)-value buckets. This yields exactly width-two adjacency:
two three-coordinate states can differ in at most two coordinates if and
only if they agree in at least one. Each bucket is a clique, so joining its
members to one representative preserves the exact component partition.

Independently of the state-graph enumeration, the checker constructs the
predicted obstruction set by taking each of the 16 labeled spanning trees,
its three two-edge subsets, and all six assignments of these distinct
subsets to the colors. It verifies **set equality** between these 96
predicted triples and the 96 disconnected triples found by enumeration.
Thus both necessity and sufficiency are exhaustively checked, not only a
one-directional test or a comparison of counts.

The original `classify_rank2.py` was also replayed successfully, and its JSON
output agrees with `classification.json`. Repeating the independent checker
reproduces its receipt byte for byte.

## Mathematical scope

This audit certifies the finite three-color quotient classification under
the stated rank-at-most-two restriction. The contraction/lifting discussion
is consistent: replacing an actual edge by another realization of its same
colored component pair preserves the quotient tree, so each type-state
fiber is connected; type moves lift while preserving unchanged actual
edges. A completion cannot use a loop after contracting the fixed forest.

This audit does not independently establish the earlier rank-at-least-three
connectivity criterion cited by the document. It also does not establish
that a cyclic obstruction can or cannot be realized by integer absolute
differences and a graceful fixed forest, nor does it independently replay
the reported interval checks through m = 8. Those remain separate claims.

## Reproduce

Run `python verify_rank2_independent.py` from this directory, or supply its
full path from any working directory. The program prints the JSON receipt
and exits successfully only after all assertions pass. Do not run Python
with assertion-disabling options (`-O` or `-OO`).
