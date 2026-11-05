# Independent structural audit: complete C3 bridge theorem

**PASS for the all-order mathematical reduction**, with no correction to the frozen proof. The independent finite-census regeneration is owned by the separate bridges reviewer and must be read with this report before the final combined theorem is marked fully audited.

Frozen source bindings:

- C3_COMPLETE.md:a919deb88bc78a84b8d0cf40199fc7966051541f1e7027be3666bb9f0713b444.
- verify_c3_complete.cpp:e49f0b9ff6e29f7b7260d6516ccb0fafe03ea4bbcf06dccc0664b7fddead806f.
- verify_c3_complete.py:88096aede3e208dfa66a45f02ae576d7459c8322841663b8627e630e24547799.
- c3_complete_certificate.json:5461864f55486797ef825342841ef30e5258f90cc74b194087573b5fb2ce341c.

This report does not claim to independently rerun or regenerate the numerical bank. The structural prerequisites C1_COMPLETE, START_TWO, BRIDGE_TWO, C2_COMPLETE and the separating-vertex product were independently reviewed in the preceding ribbon audit reports.

## All-order growth and completeness

For a proper strong induced S in a strong T, a mixed outside vertex gives a strong single-vertex extension and an injective cycle insertion. If there is no such vertex, the outside vertices split into those dominating S and those dominated by S. Both types exist, and an edge from the latter to the former exists by strongness. The resulting ear S→q→p→S is strong.

Its cycle count is exactly H(S): the two ear vertices must be consecutive q,p, and the remaining vertices form one directed Hamilton path of S. Different cuts of any old Hamilton cycle give |S| distinct Hamilton paths, and different old cycles cannot yield the same path because closing the path recovers the cycle. Thus C(S+q+p)=H(S)>=|S|C(S). Iterating these extensions proves strong-induced cycle monotonicity.

Starting from a strong S of order at least4 inside a C<=3 target, an ear would create at least4 cycles, contradicting that monotonicity. Consequently a mixed singleton extension is always available and every intermediate count remains at most3. The strong-four existence proof also checks: extend a triangle by a mixed vertex, or retain two of its vertices and the ear pair. Its unique isomorphism class is indeed[2,4,9,3]. Therefore generation from that one seed by all mixed singleton extensions of count at most3 is complete at every order. No existence of a strong deletion of an arbitrarily chosen vertex is assumed in this generation argument.

## Finite-search contract

Every nonempty proper out-neighbor mask represents exactly a mixed extension. The insertion-gap sum is a safe lower bound: deleting the newborn recovers both the old cycle and its insertion edge, so the counted images are disjoint. Rejecting only when this bound exceeds3 cannot lose a desired graph.

Exact cycle recursion capped at4 suffices to distinguish the retained counts1,2,3 from rejection. Canonicalizing over all actual Hamilton cycles and all rotations is exact: an isomorphism maps that complete set bijectively; equal complete adjacency encodings give the converse. These are valid completeness/canonicalization rules, independently of the code's eventual numerical result. The bridges reviewer is independently testing that result and implementation with a different count method and canonicalizer.

## Bridge orientation and the twelve-vertex cutoff

For a C3 bridge split, the exact cycle product is s_a(A)e_b(G)=3. One factor is1. The endpoint-one structural lemma yields an outdegree-one vertex on the G side or, after reversal, an indegree-one vertex on the A side. Reducible sides are handled correctly by their ordered SCCs; singleton factors do not hide an extra incoming or outgoing arc at the selected surviving vertex. Thus the correct conclusion is a disjunction, not both orientations.

If a bridge-free C3 graph existed, choose one of minimum order. The separating-vertex product makes it2-strong. If its order were at least13, every deletion would be strong with1,2,or3 cycles. The completed C1/C2 theorems, excluding the regular-five exception by order, or minimality at count3 would supply a degree-one vertex of at least one orientation in every deletion.

Let L be full outdegree2 vertices and R full indegree2 vertices. Minimum semidegree2 shows each deleted vertex is beaten by L or beats R. For l=|L|, at most2l−binom(l,2) arcs leave L; hence l<=5 when nonempty and |L union N+(L)|<=l(7−l)/2<=6. The empty case contributes0. Reversal gives the same bound for R union N−(R). These two sets cover every vertex, so the order is at most12, a contradiction. Possible overlap only improves the inequality; no independence or disjoint-endpoint assumption is made.

This is a bound on the smallest hypothetical atom. Combining it with a complete verified absence through12 proves that there are no atoms at any order. The proof does not incorrectly claim that the cover argument by itself bounds every possible atom.

## Scope

The six-vertex asymmetric-degree control is correct: r→C3→v→r has a start-three root whose only indegree-one vertex is that root; the added b→r destroys that indegree-one condition while leaving a sole outdegree-one vertex. Thus the C2 two-orientation equivalence and rigidity conclusions are not silently extended to C3.

Once the separate finite certificate is verified, the approved conclusion is that every strong tournament with exactly3 Hamilton cycles has a strong bridge and therefore at least one degree-one orientation. The class-count row is still a bounded inventory through12, not an all-order enumeration formula. Nothing here proves the full Hamilton-path spectrum conjecture or scalar coverage of the low-cycle classes.
