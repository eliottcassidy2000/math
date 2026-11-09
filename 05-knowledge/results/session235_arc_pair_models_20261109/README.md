# Session235: reciprocal arcs, missing pairs and the tournament boundary

**COMPLETE, independently audited.** 2026-11-09. **PROVED:** pair-state
identities, all-order parity and spectrum boundaries, and transport of the
inherited complete semicomplete core grammar. **FINITE-EXACT:** the complete
labelled oriented census through order six. This packet continues
[Session234](../session234_digraphs_20261109/README.md).

## What changes when reciprocal pairs become missing pairs?

A tournament has **exactly one direction for every pair**. Forbidding
reciprocal arcs alone gives an oriented graph, which may have missing pairs.
Those missing pairs allow seven Hamiltonian paths from order five and
twenty-one from order six. Exactly one missing pair suffices, sharply, at
every attainable order. Singleton padding supplies all larger orders.

Treating raw 11 as empty is an alternative encoding, not ordinary adjacency.
Because 00 was unspecified, the [pair dictionary](PAIR_MODELS.md) separates
three choices: forbid 00; also interpret it as empty; or interpret it as both
arcs. On the raw no-00 sector, all three give the same oriented graph.

With the last choice, decoding is the involution `tau(A)=J-I-A^T`. It swaps
empty and reciprocal pairs and fixes single directions. For every digraph,
`H(tau(A)) = H(A) mod 2`, a form of **Berge's complement-parity theorem**.
The [proofs](theory/PROOFS.md#pb4-complementreverse-parity-for-every-digraph)
give both a spanning-linear-forest identity and a square-zero-ring proof.
Parity alone does not preserve H, and on tournaments this particular
congruence is tautological rather than a proof of Redei's oddness theorem.

## Exact answers for the same raw array

Deleting both arcs of each reciprocal pair is an isomorphism-equivariant
bijection from semicomplete graphs to oriented graphs. It preserves the
entire automorphism group, although it can split strong components and
change the Hamiltonian-path count.

Transporting Session234's complete core grammar proves at **every order**:

- Every semicomplete array with ordinary H=7 decodes to H=1.
- Semicomplete arrays with ordinary H=21 decode precisely to
  **{1,3,5,7,15}** over all orders. The first four occur at order five, and
  all five occur at each order at least six.
- Ordinary H=21 and decoded H=7 occur together in exactly n-4 isomorphism
  classes and (n-4)n! labelled forms for each n>=5.

[Core transport and exact joint generating functions](census/CORE_TRANSPORT.md)
include the complete 37-core table and the argument that transports source
SCC words even when decoded blocks are no longer strong.

## Finite census and the all-order boundary

The census covers every one of the 3^(n choose 2) labelled oriented graphs
at each order n<=6. For H=7 it finds 1,440 labelled graphs at order five
(12 classes) and 255,600 at order six (355 classes). For H=21 it finds
29,520 at order six (41 classes). These are finite class censuses.

An independent direct-scatter algorithm checked **all 14,348,907 order-six
entries**, not just samples or aggregate totals. The semicomplete side has
a separate inherited full-marginal check plus endpoint checks; it was not
independently rescanned entry by entry by the second algorithm.

Over all orders, oriented graphs attain **every nonnegative integer**.
Moreover every fixed h>=2 occurs on strong oriented graphs of unbounded
order. The cycle-and-hub construction proves this also for 7 and 21, so the
finite alphabet of actual strong semicomplete cores cannot simply be
transferred to arbitrary oriented H fibres. Other compressed grammars are
not ruled out.

## Proofs, artifacts and reproduction

- [Pair-state dictionary and encoding counts](PAIR_MODELS.md)
- [All-order identities and strong-core obstruction, PB1-PB7](theory/PROOFS.md)
- [Finite census, witnesses and joint histograms](census/REPORT.md)
- [All 408 H7/H21 classes through six, explicit paths and completions](census/witnesses.json)
- [All-order joint fibres and generating functions](census/CORE_TRANSPORT.md)
- [Independent audit and exact verification scope](audit/AUDIT.md)
- [Portable full replay](audit/replay_all.py) and [deterministic receipt](audit/replay_all.json)

From any working directory, run Python 3 on `audit/replay_all.py` using the
appropriate path. It requires a C++17 compiler and keeps scratch products in
a temporary directory. It regenerates and compares 36 persistent artifacts,
including all 24 census CSVs and the complete independent oriented scatter.
The inherited Session234 generator and order-six marginal remain explicit
inputs, not claims of an independent reproof of their classification.

The recovered final packet was replayed successfully against current main,
with the same deterministic receipt as the original completed audit.
Berge parity, SCC decomposition and ordinal-sum multiplication are known
mechanisms, not novelty claims. See also current canon
[THM-5071](../../../01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md)
for the Grinberg-Stanley weighted-cycle refinement.

**The Chan-Pak tournament spectrum completeness conjecture remains OPEN.**
Changing the graph category does not settle it; neither do these results
classify all oriented graphs of a given positive H.
