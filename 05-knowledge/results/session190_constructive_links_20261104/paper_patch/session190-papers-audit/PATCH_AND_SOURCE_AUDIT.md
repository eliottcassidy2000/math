# Independent constructive patch and source audit

## Verdict

PASS for the derived boundary-signature transfer lemma and the H -> H1 demo. Both narrow source issues are real and their proposed repairs are valid. This audit does not reproduce the papers' global theorems, computer searches, or configuration libraries.

The useful extra simplification is a two-anchor deletion certificate: delete old vertex 3 to cover every new deletion in {1,3,5,7,9}, and delete old vertex 8 to cover every new deletion in {8,10,11}. Two old anchors are minimal for the displayed labeling.

## 1. The signature is complete for the claimed purpose

The degree vector and pairing of degree-1 terminals contain everything needed for gluing a spanning linear forest into a Hamiltonian cycle. A degree-0 terminal is an isolated patch component and all its cycle edges lie in the exterior. A degree-1 terminal is a connection point between patch and exterior. A degree-2 terminal has no exterior cycle edge, so its location along a particular patch path need not be recorded. After suppressing patch path interiors, the exterior sees exactly the same matching before and after replacement. Restoring those paths preserves both degree 2 everywhere and the number of cycle components, in particular connectivity.

Every nontrivial fragment component is a path with two terminal endpoints: internal patch vertices are required to have degree 2, and cycle components are explicitly excluded. This also explains why the terminal degree vector cannot be omitted: it distinguishes an isolated terminal from one traversed internally by a path.

The exterior-size hypothesis is sufficient and is used at the correct place. With at least two exterior vertices outside the terminal set, a cycle spanning the glued graph, or any one-vertex deletion of it, still visits an exterior-private vertex. Its patch restriction therefore cannot contain a closed cycle. The restriction is a permitted fragment in every transfer. This is a sufficient hypothesis, not a claim of minimality; weakening it requires separately handling cycles entirely inside the patch.

All deletions are covered without reversing an inclusion:

- Exterior-private deletion uses S(P) subset S(Q).
- Terminal deletion uses S(P-t) subset S(Q-t), with the same surviving labeled terminal set.
- New internal deletion v uses an old internal w with S(P-w) subset S(Q-v).
- Non-Hamiltonicity uses the reverse full inclusion S(Q) subset S(P).

No planarity assumption is required for this transfer argument. Disjoint patch/exterior edge sets and no shared vertices outside the terminals are essential as stated.

## 2. Actual H and H1 patches

I visually checked Figure 2, Durias et al., PDF/printed p. 5. Take vertices 0,1,...,7 clockwise around the outer 8-cycle starting at the upper-left corner. In H, 8 and 9 are the upper and lower interior vertices. In H1, 8 is upper interior, 9 central, and 10,11 lower-left/lower-right. This exactly matches the supplied edge lists. The four terminals are the corners 0,2,4,6. Remark 2.1 on that page supplies the required attachment restriction; ordinary subgraph containment is insufficient.

The independent checker `independent_boundary_audit.py` uses exhaustive edge subsets and a disjoint-set cycle test. It does not import the supplied enumeration code. It verifies:

- 14 distinct full states on each patch, with equal state sets.
- 34 actual spanning linear forests in H and 68 in H1, grouped into those states.
- Every one of the 28 supplied full-fragment witnesses is valid and has its recorded signature.
- Every terminal deletion has exactly 3 states in each patch, with the required inclusions.
- All internal-deletion counts and every listed coverage relation match the supplied JSON.
- All internal deletion sets are nonempty, so coverage is not vacuous.

The 14 states are precisely the 12 choices of one isolated terminal and a different degree-2 terminal, the remaining two terminals being paired, plus the two noncrossing perfect matchings (0,2)(4,6) and (0,6)(2,4), all terminal degrees 1. The crossing matching (0,4)(2,6) does not occur. This matches Proposition 2.2 on p. 6 and Figures 22-23 on p. 24.

The counting argument is structural, not merely empirical. In both patches all terminals lie in the larger bipartition class, whose surplus is 2. Each forest path, including a singleton terminal, has surplus 1. Hence every full fragment has exactly two components. Two singleton components would leave internal vertices forced into excluded cycles, so the only possibilities are a singleton plus a spanning path through the other terminals, or two nontrivial paths. Planarity excludes the crossing pairing.

The reduced internal-deletion certificate is:

- New 1,3,5,7,9: old anchor 3, with 4 old states.
- New 8,10,11: old anchor 8, with 5 old states.

No single old anchor is present in all coverage lists. The two classes correspond to the two bipartition colors. The audit JSON saves full and deletion states with explicit edge witnesses, extending the original artifact's full-state witnesses to all deletion states.

The supplied script was also rerun separately; its output equals the supplied JSON structurally.

## 3. Shared rings in the 4CT implementation

The issue is confirmed against the exact definitions:

- Configuration/free-completion definitions: PDF p. 5 / printed p. 4. Ring vertices may be identified in the actual graph, so repeated ring positions are explicitly possible.
- Degree-4 singleton D-reducibility: PDF p. 6 / printed p. 5, Lemma 2.1.
- Membership in the 8202-configuration set and non-touching definition: PDF p. 9 / printed p. 8, Lemma 3.2 and following paragraph. Non-touching means disjoint with no joining graph edge. The stated consequence there is noncrossing rings, not disjoint rings.
- The incorrect disjoint-ring inference and single-host-ring array: PDF p. 50 / printed p. 49, proof of Lemma 12.1, visually confirmed.

The octahedron example correctly refutes that inference. Its opposite degree-4 vertices are non-touching singleton configurations with the same four-vertex ring. Four distinct colors around that ring form a proper coloring of the remaining 4-cycle and are not 0-extendible to either singleton.

One precision upgrade is available because Lemma 12.1 literally says “some (not all)” non-colored configurations. Use the bipyramid over an 8-cycle v0,...,v7, with nonadjacent poles x,y joined to every equator vertex. The four even-index equator vertices form pairwise non-touching degree-4 configurations. Leave them uncolored, color x,y by 1,2, and color v1,v3,v5,v7 alternately 3,4. All four rings have all four colors and hence are non-0-extendible; every pair shares x,y. Choosing any two configurations is a proper subset, meeting even the literal parenthetical. This partial coloring also comes from a valid recursive reduction: remove the four even vertices and add the cycle v1-v3-v5-v7-v1, obtaining a smaller octahedron with precisely this proper coloring. Thus overlap is not an artifact of an impossible recursive input.

The incidence-list repair is valid. The total number L of ring positions is at most 18r even when many rings share a graph vertex. For a fixed complementary color partition, find the disjoint Kempe chains, then iterate chain by chain through its vertices and their ring incidences. A timestamp per ring deduplicates the (ring,chain) pair in O(n+m+L) total time; chain-by-chain order makes the timestamp implementation explicit and avoids a quadratic matrix or sorting. Retain all positions to evaluate the effect on each ring coloring, while counting each distinct chain once per ring for probabilities.

Each ring has at most 18 distinct incident chains. Store one improving assignment and maintain its live indicator and number of undecided chains. While processing a chain, scan its incidence list, compare the two conditional-expectation contributions, choose the better option, and update those rings. Dead rings can be skipped but their incidences still cost only O(L) total. The weights 2^(-u) can be scaled by 2^18 to exact integers with O(log n)-bit total scores. This gives the claimed O(n+m+L), hence O(n), per move without ring disjointness. The repair changes the local data structure, not the combinatorial probability argument.

## 4. Modular correction

Durias et al., PDF/printed pp. 20-21, derive 2 f4 + 3 f5 = 2 v - 4 correctly. Reduction modulo 3 gives v - f4 - 2 = 0 (mod 3), or v = f4 + 2 (mod 3). Both the subsequent general congruence and the corresponding first congruence in the proof of Proposition 4.1 contain the same incorrect '-1'.

For f4=5 the corrected result is v=1 (mod 3), agreeing with 28,31,34,37,40 in the table and the proposition's stated congruence. The final implication from the printed incorrect intermediate line is itself arithmetically invalid, but replacing '-1' with '-2' makes the paragraph consistent. This is a narrow textual repair; it supplies no new evidence about the exhaustive search or broader order-34 conjecture.

## Reproduction and scope

Run `python independent_boundary_audit.py` in this directory. Output is `independent_boundary_audit.json`, including complete fragment witnesses for all tested deletions. PDF renders retained here document the three visually checked source pages. No startup helper, installations, broad repository search, forbidden notes, or generic commit diffs were used.
