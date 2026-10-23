# Same block sizes do not identify the arithmetic graphs

Status: PROVED finite graph computation with an explicit edge list. This tests one natural interpretation of the analogy, not every possible enriched map.

The three prime sets have sizes5,4,6, matching the sizes of the user's three selected induced blocks in G25. But applying the same square-sum relation to the primes themselves does not reproduce those blocks:

- On {17,19,23,31,47}, the only internal edges are17–19 and17–47: a three-vertex path plus two isolated vertices, not the connected branched five-vertex block.
- On {29,41,59,71}, the internal edges are29–71 and41–59: two disjoint edges, not P4.
- On {2,3,5,7,11,13}, the internal edges are2–7,3–13,5–11: three disjoint edges, not P6.

On all15 primes together, the square-sum graph is instead a connected tree with14 edges and five leaves3,11,19,31,71. Its complete edge list is

2–7,2–23,2–47,3–13,5–11,5–31,5–59,7–29,13–23,17–19,17–47,23–41,29–71,41–59.

It therefore has no Hamiltonian path: a tree with a spanning Hamiltonian path must itself be a path, whereas this tree has five leaves. Replacing each prime by its ordinal label1..15 creates a different relation; the ordinary interval G15 has a Hamiltonian path, so ordinal relabelling does not preserve arithmetic square sums.

Removing the four Monster-only primes leaves the Baby prime graph with two components, of sizes8 and3:

{2,3,7,13,17,19,23,47}, and {5,11,31}.

The missing primes restore the path23–41–59–5 joining those components and the pendant path7–29–71. This is a precise finite graph structure worth retaining, but no implication from the group-theoretic prime partition to this particular square-sum adjacency has been proved.

Audit history: an initial informal hand-list missed the cross-block edges2–23,2–47 and7–29 and consequently predicted the wrong component count. The exhaustive integer-square checker rejected that prediction before publication. The complete edge list above and the independent audit, not the hand-list, are the evidence.

This is a useful loss ledger. A bijection can preserve the partition cardinalities, but cannot preserve the displayed within-block square-sum edges. The genuine group-action and Frobenius interpretations in the companion Monster note use additional structure, not a cardinality-only identification. The exact port representation of G25 likewise needs incidences and endpoint connectivity, not just the5/4/6 block sizes or3/4/5 edge counts.
