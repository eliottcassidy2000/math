# Tournament/graph trijections and the cut-Cayley aggregate

**Status: PROVED elementary / inherited results explicitly credited / FINITE-EXACT controls.**

Read `TRIJECTIONS_AND_TREE_GAUGES.md` for the complete proofs, operation maps, and scope ledger.

## Main results

1. An arbitrary fixed-orientation constraint graph F with m edges and c components meets exactly a fraction 2^(-(m-n+c)) of labelled switching classes, with 2^(c-1) representatives in each nonempty fibre. It is a complete one-to-one chart exactly when F is a spanning tree.
2. On four vertices, the free-host P4 and triangle-plus-isolate cubes are both eight-state cut-quotient charts. The star-host cube covers only four cut classes, twice each. The explicit good-chart transition is `(a,b,c) -> (a,b+c,a+b)` over F2.
3. Under that transition, path complement becomes one triangle-edge flip, not triangle-graph complement. Switching preserves the cut class but can change H from 5 to 1 inside this literal chart.
4. Directed Hamiltonian path markings, path-complement triples, and ordered tournament pairs with directed spanning-path agreement graphs form a genuine incidence trijection. This gives a lossless graph-coordinate route into the inherited score-fibre Gram construction.
5. Aggregating all labelled path-complement incidences by cut classes produces a weighted Cayley graph with generators `[E\P]`. For n>=3 it is connected except at n=1 mod 4, when it has exactly two parity components. It is bipartite exactly at n=3 mod 4; the exceptional two components are nonbipartite. At n=4 it is K8 minus a perfect matching, with weight4. From n=8 onward its degree is n!/2 and every positive weight is2.

The last graph is an aggregate quotient, not an equitable or graph-automorphism quotient. Its paths need not lift without switching representatives, so its connectivity does not settle original, black, or regular-return metagraph connectivity.

## Controls

`verify_trijections.py` is standard-library-only and uses explicit live checks. `verification.json` and `verification.optimized.json` are the identical normal and optimized receipts. The full labelled/constraint controls run through order5 as detailed in the proof note; all undirected Hamiltonian generator sets and their span/difference ranks are checked through order8.

`SOURCE_PROVENANCE.json` pins the inherited source checkout and source hashes. `startup.txt` records the read-only startup route.
