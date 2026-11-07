# Arithmetic graphs, comparison defects, and path-cover transport

[Proof and boundaries](PROOF.md). The ordered path-cover polynomial and gate operator extend to every loopless digraph, but a gate preserves every missing or bidirected comparison. It never turns a nontournament into a tournament.

For the literal directed `+2`/doubling graph on `1,...,N`, `N>=3`, there is no Hamiltonian path but exactly four ordered two-path covers. Its gate is strong, has exactly four Hamiltonian paths, and no Hamiltonian cycle. Symmetrizing the same arithmetic graph gives an explicit Hamiltonian path and even word counts; completing it to the numeric tournament gives exactly one. These are different transports.

A sharp order-three hostile has `P(q)=q^3`, identical to a transitive tournament, while retaining one missing and one bidirected pair. Every nested gate observation agrees forever. Even the full cover invariant cannot recognize the comparison alphabet.

Run `python3 verify.py` and `python3 -O verify.py`; both match [RESULT.json](RESULT.json). The finite audit covers all 4,166 labelled loopless digraphs through order four, arithmetic controls through order ten, and the explicit observer hostile. [Sources](SOURCES.json) are pinned to Session210. The proofs are all-order; no global spectrum completeness or literature-priority claim is made.
