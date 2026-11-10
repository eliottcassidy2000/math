# Independent audit: equivariant graph-corona bridge

Verdict: **PASS**, with one minor example boundary correction applied by the producer (the graph-cycle example requires odd q>=3). No theorem-level defect found. Chan–Pak completeness remains open.

## Mathematical audit

The source construction is genuinely equivariant, not merely a graph-invariant sketch. For an odd-order group acting freely on graph vertices, an edge stabilizer cannot fix either endpoint nontrivially, and cannot interchange endpoints because an odd-order permutation has no 2-cycle. The edge-vertex action is therefore free. The same argument applies to every unordered pair of constructed vertices. All five/six private vertex types inherit the free graph action. Choosing an incidence order only at each vertex-orbit representative and transporting it is well-defined, including when the graph has multiple vertex orbits or edges connecting those orbits.

Two long cycles share at most one edge-vertex, so their prescribed unordered pairs are disjoint. Each private triangle adds exactly one new pair. Consequently the forced pair count is 2e+6m+o, the constructed order is e+5m+o, and the completion exponent is exactly (choose(N,2)-2e-6m-o)/|A|. These count labelled completions of this prescribed partial tournament and this fixed action, not isomorphism classes or tournaments with exactly A as full automorphism group.

The positive completion tax is inherited correctly from S237. A shortest directed cycle through the private vertex a_v contains a triangle through it: the chord from its first vertex to its third either skips the second vertex or closes a triangle. Such triangles differ from the chosen private triangle and all other selected cycles. Private labels and supports distinguish the m extra cycles. The argument proves a coefficientwise lower bound, never equality or a count-preserving reduction. Choosing the extra triangles equivariantly is permitted by trivial vertex stabilizers.

The even-group obstruction concerns faithful ordinary tournament actions. By Cauchy an even-order faithful group gives a nonidentity involutory permutation, which reverses the arc of a transposed pair. It does not prohibit nonfaithful actions or signed converse actions. The proof correctly distinguishes selected-carrier graph planarity, O(2)/SO(3) faithful geometric representations, and Cayley graph planarity. In particular the regular C3 x C3 edgeless example is a valid planar selected-carrier example with noncyclic odd symmetry; it is not a planar geometric realization of the group.

The path action is free on ordered spanning Hamilton paths, yielding divisibility. The stronger exclusion of subgroup orbit counts 7 and 21 is explicitly inherited from THM-4987, whose statement includes arbitrary actual subgroups A<=Aut(T). This audit checked that dependency's statement and mechanism but does not freshly re-prove its all-order certificate dependency THM-4812. S190's actual groups proof Sections 1–2 supports the geometric distinctions. S237 REPORT Propositions A–B support the corona and tax. The bridge is an elementary synthesis of credited mechanisms; no external literature priority is certified.

## Independent executable controls

`check.py` imports no producer code, uses explicit exception gates rather than assert, and rebuilds the graph action, partial arcs, pair orbits, and completions independently.

- 73 graph/action cases: all graphs on 1, 2, and 3 vertices under the trivial group; every invariant graph for regular C3, C5, and C7 actions; every invariant graph on two free C3 vertex orbits (32 cases); every invariant graph on the regular C3 x C3 set (16 cases).
- Checks literal vertex count, exact forced-pair count, completion exponent, freeness, all group automorphisms, induced corona intersections, and extra triangles.
- Exhausts all 16 completions of the five-vertex single-graph-vertex carrier. Independently enumerates directed odd cycles with cyclic rotation normalization and counts disjoint packings; independently counts Hamilton paths by subset DP. OCF and the coefficientwise tax agree for each completion. H histogram: 9 occurs 5 times, 11 once, 13 six times, 15 four times.
- A C2 regular pair verifies the even obstruction.
- Normal and optimized runs both PASS and produce byte-identical JSON. Injecting reversal of a prescribed arc causes a nonzero exit in both modes at the explicit prescribed-arc gate.
- Producer verifier replayed independently in normal and optimized modes: both PASS, byte-identical outputs. The producer's 17-case universe is narrower and remains correctly stated.

These finite checks complement the proof; they do not establish the theorem at arbitrary group order, and do not enumerate nonabelian group actions. The all-order argument does not assume commutativity.

## Reproduction

Run `python3 check.py` and `python3 -O check.py` from this folder. Add `--inject-failure` to either command for the expected failing control. All files in this audit directory are owned by the independent reviewer. No producer files were edited or committed by this reviewer.
