# Session68 — classification with the missing decoration retained

**Status:** audited elementary consequences of cited modular decomposition,
with separate exact computational controls. No closed classification of all
prime tournaments or of the complete metagraph is claimed.

This first checkpoint contains the decorated ternary family and literal-line
boundary. The arithmetic/Frobenius/CRT lane and its random-source recovery
are reserved for a separately audited second checkpoint.

The 12-vertex 81/27 example becomes an extensible classification test:
the quotient and its child sizes stay fixed, but changing a child's ordinary
isomorphism type destroys an ancestor's rotational stabilizer. A canonical
decorated module tree retains exactly that missing coordinate. Choosing
compatible child models gives the full split automorphism extension; an
arbitrary substitution partition need not be intrinsic and can miss symmetries.

## Results and boundaries

| Result | Status and route |
|---|---|
| Canonical decorated trees classify ordinary tournaments; coherent child gauges split the exact automorphism extension | CITED structural theorem plus elementary proof, [classification](session68-classification/CLASSIFICATION.md), Sections 1–3 |
| Iterated cyclic ternary substitution classifies every ordinary realization of one C3 hypergraph at each `n=4·3^d`, with every automorphism order `3^(3^d+e)`, `0<=e<=(3^d-1)/2` | PROVED polynomial recurrence and complete support, [classification](session68-classification/CLASSIFICATION.md), Sections 4–5; palette representatives do not exhaust labelled realizations |
| All members share labelled cyclic triples and all path-cover counts, while their ordinary marked-Hamilton-path fibre sizes have ratio `3^((n-4)/8)` | PROVED transport and group action, [classification](session68-classification/CLASSIFICATION.md), Section 6 |
| The regular uniform members are not closed under literal metagraph lines: every partner changes the cyclic-triple count by an odd number | PROVED from the inherited endpoint-flux law, [line boundary](session68-root/METAGRAPH_BOUNDARY.md) |
| Full endpoint Hamilton-path matrices for the two order-12 witnesses, and Hamilton-word complement checks across labelled tournaments of orders 1–5 | FINITE-EXACT independent controls, [line receipt](session68-root/line_boundary.json) |

The ternary stabilizer polynomial is

    P_0(z)=2,
    P_(d+1)(z)=(P_d(z)^3-P_d(z^3))/3+z P_d(z^3).

Its coefficients count all ordinary realization types of the fixed hypergraph
by quotient-stabilizer exponent. Independent quotient reflections normalize
every labelled realization into a palette representative; only oriented
quotient rotations identify two normalized representatives. At depth 2, the
palette has 512 labelled tournaments, the full labelled fibre has `2^22`,
and there are 24 ordinary classes. These are three different counts.
This is an iterated cyclic-child construction, not the ordinary binary
necklace of length `3^d`: at depth 2 it has 24 types, not 60. Every coefficient
between the endpoints is positive and even, and both extreme coefficients
equal 2. The orbit-mass identity `3^((3^d-1)/2)P_d(1/3)=2^(3^d)` is checked
independently against literal actions at small depth.

## Inheritance, concurrent work, and actual increment

The native work freezes input commit
`fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8`; native files and receipts are
preserved byte-for-byte. Their source pin describes the evidence used,
not the publication parent.

During integration, concurrent commit
`d761e4cf0ff92264b043586a6160b274dae0b29c` supplied
[tournament composition classification](../tournament_composition_classification_20261026/README.md):
the sharp 4p threshold, the unique bad order 12 fibre, Sylow-adapted towers,
decorated-tree classification, path-cover invariance, and marked join hinges.
Those results and its THM-4803 correction are inherited/concurrent, not new
claims of this packet. Subsequent commit `c9a8823436e` supplies
[complete carousel fibres](../cyclictriple_composition_20261026/README.md),
including necklace/Möbius groups, converse pairing, a multiplicity-retaining
modular-flip spectrum, and pointed growth. Its modular-flip operation is not
the literal Hamilton-path-complement operation studied here.

Session68's contribution is the explicit all-depth
stabilizer polynomial and group spectrum, and the literal-line
nonclosure/flux sidecar. Neither matching symbols nor
matching counts identify two operation carriers.

## Verification

The [independent first-checkpoint acceptance](session68-audit/CHECKPOINT1_AUDIT.md)
pins the reviewed proof/checker bytes. Its separate checker uses full
128-realization orbit action, W2 orbit BFS, direct Hamilton-word DFS, and all
385 path covers of a six-vertex hostile. These routes do not import producer
implementations. No formal-assistant verification is claimed.

From this packet directory:

    python3 -I -B verify_packet.py --replay
    python3 -I -B -O verify_packet.py --replay
    python3 -I -B verify_mutations.py

The verifier checks the exact file universe, sizes and SHA256 hashes before
copying the packet into a temporary directory and running each checker in
normal and optimized isolated Python. It compares output bytes with frozen
receipts and never rewrites them. Native scripts use explicit failure gates,
not removable assertions. No network, repository import, or third-party
package is needed for the finite controls. Integrity checks are not proof
of the mathematics: the accompanying proofs and independent audit remain
separate evidence.

**OPEN:** classification of all modular-prime tournaments, a closed exact
run/allocation/port grammar for literal metagraph composition, and the
repository's unrelated LRC(14) frontier. A canonical recursive carrier does
not eliminate those obligations.
