# Independent source-bound audit: last one-apex internal wiring

22 October 2026. Verdict: PASS. No mathematical correction required. This review makes no novelty claim, no full-metagraph bridge claim, and no repository changes.

## Sources and reproduction

The exact reviewed RESULT.md, catalogue.py, catalogue.json, run-output.txt, and the companion session26-blue/RESULT.md are bound by SOURCE_SHA256SUMS. The subject generator was imported and run without modifying its files: every field and entry of the generated n=3,...,8 JSON equals the supplied catalogue.json (catalogue-rerun.txt).

A separate implementation, enumerate.cpp, literally enumerates all 2^(binom(n-1,2)) fixed-path masks for every n=3,...,8, including all 2,097,152 order-eight masks. It does not generate apex words or use Hamilton-path generation. On each triple i<j<k, cyclicity is exactly the forward-arc pattern 110 or 001, with the two consecutive fixed-path arcs forced where applicable. It intersects the vertex sets of all cyclic triangles. This intersection is exactly the set of deletion-transitive apices, because tournaments are transitive precisely when they contain no cyclic triangle. For a triangle-free tournament the intersection remains all vertices, correctly including the transitive class. Empty intersections are discarded early.

The resulting selected.txt records all selected masks and apex sets. audit.py independently reconstructs their adjacency matrices, derives deletion orders by repeatedly removing a dominating vertex, canonicalizes all available apex words and their reverse-complements, and forms literal kappa/reflection orbits. It checks every catalogue statistic, every internal edge's actual masks/endpoints/color, and every node's ordinary and merged tiling masses. All comparisons pass. Full detailed results, including all 24 order-eight masks, their apices and their labels, are in audit-results.json.

Reproduce with:

    g++ -O3 enumerate.cpp -o /tmp/last_apex_enumerate
    /tmp/last_apex_enumerate > selected.txt
    python audit.py

The Python script has fixed absolute input/output paths under the reviewed directory.

## Structural proof

The three-window hitting argument is sound. Any deletion apex from each endpoint must jointly hit all six windows. For an ordered pair a<b, hitting the first and last windows forces a<=3 and b>=6; avoiding a window strictly between them forces b-a<=3. Thus a=3,b=6. A singleton cannot suffice. This also excludes a transitive endpoint: choosing its arbitrary deletion apex would contradict the same constraint. Holding either endpoint's apex fixed and varying the other proves both apex sets are singletons, at opposite positions 3 and 6.

For apex 3 in T and apex 6 in T', their common six-vertex base has only the matching agreement arcs 1->2,4->5,7->8. If an agreeing pair has a vertex between its endpoints in a transitive order, reversing the two cross arcs while retaining that pair produces a cyclic triangle in the other tournament. Hence each pair is an adjacent block and all cross-block order reverses. The retained constraints 5->6->7 and 2->3->4 respectively give 45<78 and 45<12 in T's base order. Exactly two block orders remain: 45,12,78 and 45,78,12.

In the first order, vertex 6 can occupy any of the three cuts between 5 and 7 (before, inside, or after block 12). In the reversed companion order 78,12,45, vertex 3 has just the cut between 2 and 4. In the second order, these counts are 1 and 3. The insertion data determine disjoint incident arc sets except the unconstrained 3--6 arc. Both choices of this arc work; all path arcs are satisfied. Thus there are exactly (3+3)*2=12 masks for this apex assignment, and twelve with the assignments interchanged. Literal enumeration confirms the exact same 24.

Any anti-isomorphism must fix the unique apex and reverse the remaining transitive seven-vertex order, fixing its middle vertex too. Two distinct fixed vertices are impossible for an anti-isomorphism. All involved endpoint classes are therefore non-self-converse, so these internal edges are black.

## Actual multigraph verification

The six specified orbits, with minima 2247,2511,2503,4095,6343,6599, all have four elements and match the catalogue exactly. Kappa is bitwise complement of the non-path mask; reflection is its stated coordinate permutation. These commute. Reflection cannot equal kappa on a mask, because the coordinate (1,n) is reflection-fixed, so no omitted anti-fixed orbit type affects the color accounting.

The endpoint graph is exactly:

    0001010 -- 1000010 -- 0010110
    1000110 -- 0010010 -- 1001010
    0011010 == 1010010

The last pair carries two DISTINCT black orbits (minima 6343 and 6599), not one edge counted twice. There are 47 other isolated nodes and no loops. Thus this is two three-vertex paths and a double-edge component, on eight touched vertices among 55 nodes. The four path edges are bridges of this induced multigraph only. Neither parallel edge is an induced bridge. No conclusion about full-metagraph bridges follows.

## H labels and boundary

Independent subset dynamic programming counting directed Hamiltonian paths yields the following H values for the eight canonical merged-word representatives:

    0001010: 9       0010010: 15
    0010110: 15      0011010: 19
    1000010: 51      1000110: 45
    1001010: 53      1010010: 61

Each representative has a unique apex, hence is rigid: an automorphism fixes that apex and every vertex in the unique transitive deletion order. Thus these H counts equal ordinary tiling masses; converse doubles each merged mass here. All labels match the JSON.

Boundary edges were counted directly as distinct actual kappa/reflection orbits with exactly one selected endpoint, independently of subtraction. The result is 18 blue and 1928 black. This agrees with (3898-18)/2 - 2*6 = 1928. Each internal black edge consumes two selected black incidences; loops, where present at smaller orders, likewise consume two incidences but remain one edge in the edge/loop columns.

## Small catalogue completeness

The source matrix routine correctly constructs every tournament with a specified deletion apex on a transitive base. Every selected tournament has such an apex, so appears among the binary words. Each fixed-path presentation is a relabeling along a directed Hamiltonian path of one of those rooted representatives. The visit recursion enumerates every permutation satisfying exactly those path constraints, without pruning valid paths. Mask deduplication removes repeated root/path representations without removing different literal tilings. Canonical minimization over all transitive deletion apices is a complete isomorphism invariant within this family; converse acts by reverse-complement. The internal-edge lookup tests the literal companion, not merely a class-level approximation.

All complete small-order statistics are independently verified, including loops at orders 4,5,6 and their single-edge counts. The independently obtained totals are exactly the six rows of RESULT.md; no assumptions about rigidity or quotient size were used to enumerate masks or edge orbits.

## Higher-order independence

The companion's disjoint consecutive-triple argument proves tau(T)+tau(T') >= floor(n/3). For two endpoints with tau<=1, n>=9 makes the right-hand side at least three, a contradiction. Consequently this induced stratum has no edges or loops for n>=9. This proof is independent of the finite catalogue. This audit does not re-audit all other closed formulas in the companion note.
