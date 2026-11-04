# Narrow source issues for independent review

These notes identify two local textual issues in the uploaded versions. They are not a claim to have verified or refuted either paper's main theorem. The constructive research synthesis is in a separate file.

## 1. Shared rings in the 4-coloring implementation

**Source:** Inoue, Kawarabayashi, Miyashita, Mohar, Thomassen, Thorup, *The Four Color Theorem with Linearly Many Reducible Configurations and Near-Linear Time Coloring*, uploaded arXiv:2603.24880v3. Lemma 3.2 and definition of non-touching: printed p. 8 / PDF p. 9. Lemma 12.1 and its implementation: printed p. 49 / PDF p. 50.

**Exact relevant hypotheses:** Configurations are vertex-disjoint and have no graph edge joining them. The configurations' free-completion rings have length at most 18. In Lemma 12.1 the selected configurations are not yet colored, and their ring colorings are not 0-extendible. The proof asserts: “Because the active reducible configurations are non-touching, their rings are strictly vertex-disjoint.” It then proposes one host-ring ID per actual graph vertex.

**Counterexample to that inference:** Let G be the octahedral graph, with opposite vertices x,y and equatorial cycle a,b,c,d. Its edges are ab,bc,cd,da and all eight edges from {x,y} to {a,b,c,d}. G is a planar triangulation. Each singleton {x},{y} is a degree-4 D-reducible configuration from the paper's set, and they are non-touching because xy is absent. Both rings are the same equatorial 4-cycle. Color a,b,c,d with four distinct colors; then neither singleton's ring coloring is 0-extendible. Thus even the active-ring hypotheses do not force disjoint rings.

**Proposed local repair, proved at the data-structure level:** Replace the single host-ring array with incidence lists. For each actual graph vertex v, store all ring positions (i,j) whose ring vertex maps to v. Repeated positions and shared vertices are allowed. If there are r rings, total incidence L is at most Σ_i |R_i| ≤ 18r. For any chosen complementary color-pair partition:

1. Find its Kempe chains in O(n+m) time.
2. Traverse each vertex and its ring-position incidences. Deduplicate each (ring, chain) pair using timestamps. This creates at most L ring-chain incidences in O(n+m+L) time.
3. For each ring, enumerate its at most 2^18 choices of chain swaps to select an improving pattern, using the paper's precomputed levels. This is O(r) with a fixed large constant.
4. Process each chain once. Scan the live rings incident with it to compare the two conditional-expectation choices and update the affected live indicators and remaining-chain counts. Each ring-chain incidence is visited O(1) times. Weights 2^{-u_i}, 0≤u_i≤18, can be scaled to integers.

The resulting cost is O(n+m+L)=O(n) for a planar graph and r=O(n). Ring disjointness is unnecessary. This repairs this specific complexity argument; it does not independently establish other algorithmic or computer-check claims in the source.

## 2. Modular arithmetic in the quadrilateral/pentagon count

**Source:** Durias et al., *Local modifications in planar hypohamiltonian graphs*, uploaded arXiv:2610.21210v1, printed/PDF pp. 20–21, including Proposition 4.1.

The paper correctly derives 2f4 + 3f5 = 2v − 4 and f5=(2v−2f4−4)/3. It then writes v−f4−1≡0 and v−1≡f4 (mod 3). Those congruences do not follow. The correct relation is

v − f4 − 2 ≡ 0 (mod 3), equivalently v ≡ f4 + 2 (mod 3).

For f4=5 the correct consequence is v≡1 (mod3). Thus the tabulated values 28,31,34,37,40 and the statement of Proposition 4.1 agree with the corrected calculation. This is a narrow arithmetic typo, not evidence against that proposition's stated conclusion or its search result.
