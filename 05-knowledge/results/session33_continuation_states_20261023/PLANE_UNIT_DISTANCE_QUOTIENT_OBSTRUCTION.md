# A plane unit-distance graph can locally cover a complete graph of arbitrarily high chromatic number

Status: PROVED all-dimensional construction and obstruction; FINITE-EXACT rational128-point witness. This is a deliberately negative transfer result for the Hadwiger–Nelson direction, not an improved bound on the chromatic number of the plane. Coding-theoretic covers and the graph-homomorphism inequality are classical mechanisms; this note supplies an explicit joint witness connected to the Fano tournament chart.

## 1. An exact rational geometric realization of the Fano flip cube

For i=1,...,7 take rational unit vectors

v_i=((1−i²)/(1+i²), 2i/(1+i²)).

For m in {0,1}^7 put p_m=Σ_i m_i v_i in R². All128 positions are distinct, and

|p_m−p_n|=1 iff m,n differ in exactly one coordinate.

The forward implication is the nontrivial finite part: the checker examines all8128 unordered pairs by exact rational arithmetic. Thus the induced unit-distance graph on these points is exactly Q7, with448 edges and a two-coloring by bit parity. “Plane unit-distance” is geometric; this does not say Q7 is a topologically planar graph.

Using the seven nonzero F2³ normals from the Fano note, s(m)=Σ_i m_i a_i gives a graph covering onto K8. At each point its seven neighbors map bijectively to the seven other quotient vertices. Yet the source has chromatic number2 and the quotient has chromatic number8.

Moreover K8 itself cannot be a unit-distance graph in R². Eight pairwise unit-separated points would give seven difference vectors with Gram matrix having diagonal1 and off-diagonal1/2, which has rank7. Their Euclidean embedding dimension must be at least7.

The quotient therefore preserves every local adjacency star while destroying both the two-coloring and realizable planar unit-distance geometry. Its closed triangles lift to open cube paths with nonzero Hamming-code displacement, exactly the holonomy detected in the Fano chart.

## 2. This is an unbounded family, not an isolated eight

Inheritance: `01-canon/theorems/THM-432-unit-distance-product-crossover-hamming-tie.md`, section A, already proves the more general generic-angle Cartesian-product construction. The cube case is its specialization to repeated K2 factors. The short proof below is a self-contained rederivation, not a new existence theorem. No edge-density, kissing-number or extremal claim elsewhere in that note is imported.

**Proposition.** For every d>=1, Q_d has an induced unit-distance embedding in R² with distinct vertices.

Proof by choosing unit directions inductively. Assume all subset sums of v_1,...,v_d are distinct and every difference with at least two nonzero coefficients in {−1,0,1} has norm different from1. Choose a new unit vector v. A new collision has the form w±v=0, where w is one of finitely many old signed sums. Since w≠0 whenever its coefficient vector is nonzero, each such equality excludes at most two points of the unit circle. A new undesired unit distance has |w±v|=1 with w≠0, equivalent to w·v=∓|w|²/2, which excludes finitely many directions. The case w=0 is precisely a desired new cube edge. Avoid the finite union of excluded directions. Old relations are unchanged. This gives the next embedding. Rational points are dense on the circle, so rational directions can be chosen at every finite stage. QED.

Now choose d=2^r−1 and index the directions by all nonzero vectors a in F2^r. The syndrome map m→Σ_a m_a a is a regular graph cover Q_d→K_(2^r). Hence a two-colorable plane unit-distance graph can have a locally bijective quotient of arbitrarily high chromatic number.

No fact about an infinite simultaneous embedding is needed: the statement is finite for each d.

## 3. What a valid color-forcing transfer would need

A graph homomorphism G→Q pulls a coloring of Q back to G, giving χ(G)<=χ(Q), not χ(G)>=χ(Q). Thus our high-chromatic quotient cannot certify a high lower bound for its source. Local bijectivity does not reverse this inequality.

A homomorphic section Q→G would reverse the inequality and force equality, but none exists here: K8 contains triangles while Q7 is bipartite. More generally the lifted odd-cycle obstruction blocks a section. A geometric coloring argument needs its actual unit-distance constraints or a legitimate lower-bound-preserving map, not merely a quotient's chromatic number.

This connects three prior themes through actual maps:

- Fano/Hamming: the syndrome forgets a deck-coordinate.
- Tournament continuation: specified line operations descend, but closed words can change the representative.
- Hadwiger–Nelson: after forgetting that coordinate, local adjacency survives while the global coloring problem changes drastically.

It does not relate the numerical gaps7,21 to a plane coloring bound, nor equate graph-cover holonomy with the different affine/clock obstruction used in LRC.

Reproduction: `python3 check_plane_cube_cover.py > plane_cube_cover.json`, also under `-O`. All8386 explicit tests agree exactly; no floating-point tolerance or randomized search is used.
