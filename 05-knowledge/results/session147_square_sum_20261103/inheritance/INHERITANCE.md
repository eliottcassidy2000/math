# Square-sum and metagraph inheritance (2026-11-03)

Source repository: `eliottcassidy2000/math`, commit `f139791ef3f194128283655cde6c1e77e2504337`. This note is research triage, not an amendment to canon. No source/index changes. The central mistakes ledger was not accessed. The three-list identification below is independently verified, not recovered as a statement of the user's intention.

## An exact 3-4-5 square-sum decomposition

Define simple graph Q_n on {1,...,n}, with distinct vertices adjacent iff their sum is a perfect square, and delete S={13,14,15}. Direct enumeration, reproducible by checking squares in [2,2n], gives Q_15-S as three path components:

* 10–6–3–1–8 (5 vertices, 4 edges),
* 12–4–5–11 (4 vertices, 3 edges),
* 2–7–9 (3 vertices, 2 edges).

So component orders are 5,4,3 and 5²=4²+3². Component edge counts are 4,3,2, not 5,4,3. This is a concrete Pythagorean numerical analogy; no general Pythagorean law or link to tournament metagraph is thereby established.

For n=18 the same S produces exactly the three requested vertex sets. First becomes 10–6–3–1–8–17, second stays 12–4–5–11, third has edges {2–7,7–9,7–18,9–16}. In particular the last listed order 2,7,9,18,16 is a component listing, not a square-sum path: 9+18=27 and 18+16=34. Vertices 16,17,18 attach respectively to 9,8,7. The graph Q_18 has three leaves (16,17,18), so cannot possess a Hamiltonian path. The observed list match is proved, but its provenance or whether it expresses the original author's intended meaning is unconfirmed.

In Q_15, Q_15-S has 12 vertices and 9 edges; Q_15 itself has 15 vertices and 15 edges. Source `01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md` states independently audited: Q_15 has unique Hamiltonian path 8,1,15,10,6,3,13,12,4,5,11,14,2,7,9; the only unused edge is {1,3}. It proves the degree law and cites the zigzag threshold theorem, explicitly distinguishing mechanism from numerology. Do not conflate Q_15-S components with a direct decomposition of that unique Hamiltonian path: {1,3} is unused in the latter.

## Ninja/terminal canon: meaning and scope

`01-canon/theorems/THM-4850-ninja-pairs-of-the-square-sum-graph-start-at-42-and-the-base-49-numeration.md` defines a weak B pair to be Hamiltonian paths of Q_n and Q_(n+1), both beginning 1,3, same positional parity for shared vertices, with even-order member ending 8 and odd-order member ending 10. It reports existence iff n>=42 using exhaustive small cases and cited bases and recursion; dropping prefix/end conditions gives Gerbicz nice pairs. It expressly marks identifying the owner's “ninja pair” with this mathematical term as INTERPRETATION; don't silently call that identification proved. Different normalizations have different thresholds: free nice pairs exist at {15,16,25,26,27} ∪ [29,∞), and Gerbicz-normalized ones for n>=41. The unique base-49 representation uses b∈[42,2081] and digits 24,...,72. Consecutive coded objects are not asserted related by simple insertion.

`01-canon/theorems/THM-4851-the-prehistory-of-the-square-sum-numbers-unique-cycles-at-32-and-33-terminal-spectrum-with-the-hole-at-39.md`: Q_N has Hamiltonian path exactly N∈{1,15,16,17,23}∪[25,∞); Hamiltonian cycle exactly N>=32 (cited Gerbicz). Its terminal cycles through 8–1–3 (even) or 15–1–3 (odd) exist exactly N∈{34,...,38}∪[40,∞): exhaustive through 62, independent bounded search for 63..2620, and cited recursion beyond. 39 is a real exception (50 cycles but none terminal). Reports independent audit with corrections. Do not describe bounded computation alone as an all-n theorem.

## Tournament metagraph conventions and independent enumeration

`00-navigation/METAGRAPH-ATLAS.md` explicitly takes G_n vertices to be isomorphism classes of n-tournaments: counts at n=3,4,5 are 2,4,12. It calls wiggly the distance-one flip layer and waggly all Hamming layers. An exploratory enumeration suggests simple edge counts 1,5,30 for n=3,4,5 under one-edge flip adjacency with loops and parallel edges suppressed. This calculation was not independently audited here and is not claimed as a result. The loop/multiplicity/distance/duality convention is essential. The raw match 18 versus Q_18, and Q_15-S orders 3,4,5 versus index levels 3,4,5, do not by themselves give a graph correspondence.

## Strong constructive seam and remaining proof target

`05-knowledge/results/square_sum_spectrum_overnight_20261031h/HANDOFF.md`, plus its README and PYTHAGOREAN-SEAMS.md, reports an exact two-seam cycle-count product for selected sparse subcubic square-sum carriers with same-square-class marked ports, using square scalings and a reverse two-edge-cut bijection; tested classes 1,2,5,10. It explicitly warns these are not full induced graphs or initial intervals. Port shape compatibility can require duplication/square-scaling; this is not a one-copy heterogeneous join theorem. Next meaningful target is a verified map from the 3–4–5 path components and their attachment ports to a precise graph operation, and/or a phase/interval construction preserving induced edges and ownership; equality of integers or sizes is not such a map. Source flags the latter and one-copy joins as open.

## Limitations

Attempted literal full-tree and narrow `git grep` searches blocked on promisor-object retrieval (proxy CONNECT aborted); they returned no authoritative provenance for the three lists. Do not infer from the unsuccessful search that no source exists. No external source was consulted. Treat all other speculative analogies as conjectural, and preserve the distinction between proved, cited, finite checked, and interpreted.
