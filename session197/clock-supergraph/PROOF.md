# The minimum ternary clock cannot reach its +4 child by any induced extension

## Exact result

Let

`T = [2,12,145,149,99,143,175,19]`.

Rows encode outgoing edges by zero-based bit positions. This is the strong eight-vertex seed of the returning `60*3^m-1` transitive-clone clock, with H(T)=179.

**There is no tournament of Hamilton count 541 containing T as an induced subtournament, at any order.**

This is stronger than the separately proved mixed three-port failure. It has no port restriction at all: even arbitrary new vertices with arbitrary orientations to every old vertex cannot produce `3*179+4=541` from this exact induced seed. The +2 child 539 exists by adding a vertex beating exactly old vertices 0,1,4.

The assertion is about this graph, not merely its count. Section 4 gives another strong graph with count 179 and a common physical four-port interface that produces both 539 and 541 in one birth. That is a finite, explicit rebasing escape; no returning grammar is claimed for that replacement graph.

## 1. Exact all-order reduction

Suppose a tournament U contains T induced and H(U)=541. Since T is strong, its vertices lie within one strong component C of U. Induced-count monotonicity gives H(C)>=179. Tournament Hamilton counts multiply over the ordered strong-component chain. Because 541 is prime, H(C)=541 and every other strong component has count one. Only C is needed for the contradiction; it still contains the entire fixed induced seed.

Busch's strong-order bound gives `5^(|C|-1)<=H(C)^3`. The exact comparison

`5^11 <= 541^3 < 5^12`

forces |C|<=12. Therefore any alleged unrestricted extension has an induced witness obtained by adding at most four vertices to T.

The invoked theorem is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, Theorem 1, DOI 10.37236/1141:
https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf
The stated piecewise bound is stronger than the safe inequality used here.

## 2. Complete unrestricted extension census

Start with the fixed labelled seed on vertices 0 through 7. At each step add a new vertex in EVERY orientation to the entire current graph, including a universal source and a universal sink. Weak intermediate graphs are not omitted. Reject only graphs of count greater than 541; induced-count monotonicity makes this cap safe for every intermediate induced subgraph of a possible final witness.

Deduplicate by fixing all eight seed labels pointwise and minimizing the adjacency tuple over every permutation of the added vertices. This is exact isomorphism relative to a fixed embedded seed, not a degree signature or approximate invariant. It need not identify isomorphisms that move the seed: retaining extra representatives is harmless. Every possible induced copy of T can first be relabelled to the fixed seed labels.

Every extension of a parent identified under this seed-fixed isomorphism appears among the unrestricted orientations of its representative. Hence this generation is complete at each order by induction. There are at most four added vertices, so the canonicalization checks at most 24 permutations.

The exact inventories are:

| Order | All tested extensions | Accepted before isomorphism | Seed-fixed classes with H<=541 |
|---|---:|---:|---:|
| 8 | initial seed | initial seed | 1 |
| 9 | 256 | 33 | 33 |
| 10 | 16,896 | 136 | 68 |
| 11 | 69,632 | 311 | 105 |
| 12 | 215,040 | 564 | 144 |

None of the retained graphs has count 541. The finite census plus Section 1 proves the all-order theorem. It is not an extrapolation from a bounded-depth unsuccessful search.

## 3. Consequences and limits

The same exclusion holds for the converse of T, since reversing all arcs preserves Hamilton counts and induced containment. It also holds for any graph containing either of those cores, including neutral singleton padding of the clock source.

Thus a compiler that keeps this seed induced cannot supply both first-step ternary children. Adding more ports, re-exposing forgotten vertices, or allowing arbitrary orientations of new vertices cannot fix this specific failure. Rebuilding the graph, changing old edges, or choosing a different count-179 representative is outside the exclusion.

No claim is made about the +4 children of later clone multiplicities, or about arbitrary count-179 tournaments. The only excluded target here is 541 for the specified induced seed and its stated consequences.

## 4. A physical four-port rebasing escape

Take instead

`G = [2,64,11,51,39,135,61,95]`,

with common exposed ordered ports `[0,1,2,5]`. This graph is strong and H(G)=179. Both following one-birth moves retain all four old ports:

- Left birth losing only to old vertex 2: `L(mask4,drop4)`. The resulting rows are `[2,64,267,51,39,135,61,95,251]`, with H=539.
- Right birth beating exactly old vertices 0,1,5: `R(mask11,drop4)`. The resulting rows are `[2,64,267,307,295,135,317,351,35]`, with H=541.

Every other incidence is the appropriate default of that side. Both children are strong because each newborn has both an in-neighbor and an out-neighbor in the strong parent. The verifier checks the three counts both by ordinary subset dynamic programming and by literal permutation enumeration.

This is one common actual exposure, not two unrelated chosen pairs, and gives both scalar children from the same actual parent graph. It is only a one-generation branching certificate. It does not assert that either child returns to a class with both moves available, or that a physical birth word changes T into G. The point is precisely that equal Hamilton count does not imply equal extension reachability.

## 5. Reproduction

Run `python3 verify.py` and `python3 -O verify.py` in this directory. A C++17 compiler is required. The C++ generator compares forward and backward endpoint recurrences on all 301,824 extension candidates. It writes a deterministic complete representative bank to temporary scratch. Python independently checks every retained graph's count, embedded seed and exact canonicalization, and independently regenerates the entire one-vertex bank. The receipt records the bank SHA-256, all inventory counts, and both rebasing moves.

No large certificate file is a prerequisite: the complete bank is regenerated and hashed. All checks use explicit exceptions or runtime errors rather than removable assertions. The original port-pruning snapshot and earlier clock packets are unchanged.
