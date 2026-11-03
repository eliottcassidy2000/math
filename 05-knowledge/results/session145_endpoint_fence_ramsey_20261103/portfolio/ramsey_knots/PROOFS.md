# Sparse Ramsey barriers and endpoint-gluing loop counts

Session 145, 2026-11-03. **PROVED elementary deductions / FINITE-EXACT controls / CITED external bounds.** No claim of a new unrestricted Ramsey bound or knot invariant, and no historical priority claim.

## Findings

1. For a red/blue/green complete-graph coloring avoiding red and blue triangles and a green K_k, the largest possible order is exactly 2(k−1) if the union of the red and blue graphs is a forest; exactly 3(k−1) if that union is a pseudoforest; exactly 4(k−1) if one of the two graphs is a forest and the other is a pseudoforest; and exactly 5(k−1) if both are pseudoforests.
2. For r disjoint paths with 2r named endpoints and any fixed endpoint matching P, the loop-generating polynomial over every matching Q used to close them is

   sum_Q z^{c(P,Q)} = z(z+2)(z+4)...(z+2r−2).

   In particular 2^{r−1}(r−1)! closures produce one component. This exact component observer forgets knotting, over/under data, and framing.
3. The corrected Paley-11 endpoint data are recovered and independently recounted: (a,b)=(628,1101), a+b=1729, H=95095, and C=5505 unrooted directed Hamiltonian cycles. The number of Hamiltonian words containing any specified directed arc is 17290, not 1729.
4. Current primary material inspected on 2026-11-03 supports 45≤R(3,3,5)≤57 and 61≤R(3,3,6)≤86, with 86 explicitly a computer-assisted preprint bound by Dauth. Both exact values remain open. The April 2026 survey's upper bound 91 is not the current frontier; see SOURCES.md.

## Inheritance and object discipline

The closest attachment mechanism is the endpoint-block polynomial carrier in `primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md`, alongside Session86's distinction between a full named weighted relation and contextual quotients. The closest topological warning is `kernels_games_detours_20261020f/walks_knots.md`: a residue/permutation quotient forgets integer braid winding. The Ramsey inheritance is the saturated-port packet and its neighbourhood-counting correction. The 1729 inheritance is Session54's corrected endpoint operator, not an interpretation of an old scalar ratio.

The six live objects are: oriented Hamiltonian paths; signed endpoint-block polynomials; endpoint matchings; actual embedded tangles; sparse unions of two forbidden colors; Paley endpoint counts. Their observers are different. In particular a polynomial's monomial support being a forest does **not** say the original tournament or a Ramsey color graph is a forest. The endpoint-block dimensions 7,28,126 are rational linear dimensions of averaged observations, not counts of perfect matchings or knot types. No all-order minimality of that carrier is assumed in this file.

## A. Exact restricted Ramsey thresholds

A pseudoforest is a finite simple graph in which every connected component has at most one cycle. Let C be a coloring of every pair of n vertices red, blue or green, and put H=R∪B. A green clique is exactly an independent set of H. Thus avoidance of a green K_k says alpha(H)≤k−1.

**Lemma A1 (four chromatic bounds).**

- If H is a forest, chi(H)≤2.
- If H is a pseudoforest, chi(H)≤3.
- If R is a forest and B is a pseudoforest, chi(H)≤4; the same holds with the colors exchanged.
- If R and B are pseudoforests, chi(H)≤5.

**Proof.** A forest is bipartite. A unicyclic component can have its cycle colored with three colors and its attached trees colored outward; tree components need two. For the third case, every nonempty induced subgraph on s vertices has at most (s−1)+s=2s−1 edges. Its average degree is below four, so it has a vertex of degree at most three. Repeated deletion and reverse greedy coloring use at most four colors. In the fourth case the corresponding edge bound is 2s, so there is a vertex of degree at most four and reverse greedy coloring uses at most five colors. Forest and pseudoforest properties survive restriction, which is essential to these degeneracy proofs. QED.

Each proper c-coloring partitions V(H) into c independent sets. Hence alpha(H)≥ceil(n/c), and n≤c(k−1) in each respective class.

**Sharpness.** Take k−1 disjoint blocks, and color all cross-block pairs green. Within a block use the following complete red/blue colorings, with vertices numbered from zero:

- K2: red edge 01, no blue edges.
- K3: red path 0–1–2, blue edge 02.
- K4: red path 0–1–2–3; blue path 2–0–3–1.
- K5: red cycle 0–1–2–3–4–0; blue complementary cycle 0–2–4–1–3–0.

Every red/blue graph has no triangle; the relevant forest/pseudoforest promises hold globally because blocks are disjoint. An independent set of H contains at most one vertex per complete block and can select one from each, so alpha(H)=k−1. These attain all four upper bounds for every integer k≥2. QED.

If the restricted Ramsey threshold is defined as the first n where every coloring **in its stated class** has one of the forbidden cliques, its exact values are 2k−1, 3k−2, 4k−3, and 5k−4. The maximum avoiding orders at k=5 are 8,12,16,20; at k=6 they are 10,15,20,25. These do not change the unrestricted Ramsey numbers. Rather, even the most permissive of these sparse architectures cannot reach the known unrestricted witness orders 44 and 60.

**Hostile boundaries.** A bichromatic triangle is a pseudoforest union and violates the stronger forest factor two. The displayed K4 is the union of two forests, but is not itself a pseudoforest and violates factor three. The displayed K5 has two pseudoforest color classes and violates factor four; neither color class is a forest. “At most one cycle per component” must not be mistaken for “at most one cycle in the entire graph.” Conversely the red/blue triangle-free hypotheses are unnecessary for the upper estimates, but are necessary to verify that the sharp witnesses belong to the Ramsey-avoiding universe.

**Connection map.** Source: an actual complete three-coloring, plus the declared sparsity promise on its actual edges. Target: a proper coloring of H and its independent-set certificate. Map: forget the distinction between red and blue, then use Lemma A1. Preserved predicate: an independent set is a green clique. Lost data: the two separate triangle predicates and all attachment multiplicities. Required sidecar: retain the original two colors to check membership of the avoiding universe. Cheapest decisive tests are K3, K4 and K5 above.

### A2. A genuine bridge from attachment-polynomial supports

The inherited forest-support theorem says a monomial in one fixed attachment Hamilton-count polynomial is supported on a linear forest of the old vertices. In a product of r such observations on the same old source, every support F is contained in a union of r linear forests. Reduction on the Boolean cube can delete repeated edges; it cannot add edges outside that union.

**Support-coloring lemma.** Such a graph F has chi(F)≤2r for r≥1. Consequently, if a proposed Ramsey construction takes its non-green edge set to be F, its green graph contains a clique of order at least ceil(n/(2r)). Avoiding a green K_k requires

   n≤2r(k−1).

**Proof.** For every nonempty induced vertex set S, each forest contributes at most |S|−1 edges. Thus e(F[S])≤r(|S|−1)<r|S|, so some vertex has degree at most 2r−1. Reverse greedy coloring proves the claim. The independent-set and green-clique conclusions then follow as in A1. QED.

For a support construction reaching the known witness orders 44 at k=5 or 60 at k=6, this necessary bound forces r≥6 in either case. It gives neither a sufficient construction nor a lower bound for the averaged carrier dimension. Assigning the support's edges to red and blue while avoiding monochromatic triangles is a further obligation. Without the explicit assumption that actual non-green edges are that support, the lemma makes **no claim** about the source tournament or any Ramsey graph. In particular a signed coefficient or its cancellation is not a positive graph constructor. At r=2 the K4 pair of paths above shows that the factor four is sharp even with the triangle-free color split.

## B. The exact gluing quotient for paths and tangles

Start with r disjoint nontrivial paths, each with two distinct named endpoints. Collapse each path to its endpoint pair; the result is a perfect matching P on the 2r ports. A closure specifies another perfect matching Q, realized by fresh internally disjoint connector paths. Its abstract components are the connected components of the multigraph P∪Q. A common edge of P and Q is retained twice: it gives a legitimate two-edge cycle in the collapsed multigraph, not a graph-theoretic loop and not zero components.

All vertices have degree two, counting multiplicity, so each component alternates P and Q edges. Write c(P,Q) for the number of components. If p,q are the fixed-point-free involutions encoding the matchings, then

   c(P,Q) = number_of_cycles(p∘q)/2.

Indeed stepping twice around an alternating component separates its vertices into two cycles of p∘q, including the doubled-edge case where the two vertices are fixed separately. This provides an independent count.

**Theorem B1 (loop polynomial).** Put F_0(z)=1. For any perfect matching P on 2r ports,

   F_r(z):=sum_Q z^{c(P,Q)} = product_{j=0}^{r−1}(z+2j).

**Proof.** Relabeling ports takes any fixed P to any other, so the sum is independent of P. Select a P-edge ab. If Q contains ab, removing a,b removes one closed component, contributing zF_{r−1}. Otherwise Q pairs a with x and b with y for distinct remaining x,y. Remove a,b and replace ax,by with xy. This gives a matching Q' on the remaining ports and does not change the number of components: the alternating segment x–a–b–y is shortened to x–y. Conversely choose any of the r−1 edges xy of Q' and either orientation of its ends to recover Q. This contributes 2(r−1)F_{r−1}. Thus F_r=(z+2r−2)F_{r−1}; the product follows. QED.

At z=1 the product is (2r−1)!!, the total number of matchings. The coefficient of z is 2^{r−1}(r−1)!, the number producing one component. The coefficient of z^r is one, because maximal loop count occurs precisely when Q=P.

**Corollary B2 (exact contextual minimum).** For the grammar allowing every named matching closure Q and observing its number of components, two states have identical observations exactly when they have the same P. Sufficiency is the definition of c. If P≠P', choose Q=P: the first loop count is r and the second is smaller. Hence this deterministic quotient has exactly (2r−1)!! states. This is not a rational linear-dimension statement, not a bit count for arbitrary encodings, and not a claim about a grammar restricted to planar/noncrossing closures.

**Endpoint blocks.** Each signed path-block coordinate Q_pi from the inherited averaged carrier has a matching between the two endpoints in each nonsingleton block. Every concrete internal Hamilton path order in that coordinate has the same abstract matching. Thus B1 can evaluate its connector-component statistic independently of its internal order. The full coordinate retains the named block vertex sets and signs; this matching quotient discards them and is not sufficient for future Hamiltonian attachment responses. Singleton blocks are not two-ended paths and must be treated separately, not given fictitious distinct endpoints.

## C. What the gluing quotient cannot say about knots

For actual tangles consisting of exactly r arc components and **no pre-existing closed components**, retain their embedded strands, boundary identification and chosen connector embedding. Under this explicit scope restriction, the matching determines the number of components after gluing because following each strand gives exactly P∪Q. More generally, if the input tangle also has m closed components, the total after closure is **m+c(P,Q)**; the carrier must retain m as well as P. Adding a disjoint closed circle is an immediate hostile to a claim using P alone for unrestricted tangles. It does not determine ambient isotopy, linking, crossing signs, integer twisting or framing. A rotation system for an abstract graph alone does not restore missing over/under information in a spatial knot diagram.

A concrete hostile is the pair of two-strand braids sigma_1 and sigma_1^3 with the usual braid closure. Both have the same endpoint transposition and both closures have one component. Their closures are the unknot and trefoil, respectively. The first loses its only crossing by a Reidemeister-I move. They can be distinguished without relying on names: Fox 3-coloring counts are three and nine.

For completeness, at a crossing the rule is 2b=a+c in F_3, with b the overstrand color. The operation a*b=2b−a has

   a*a=a,   (a*b)*b=a,   (a*b)*c=(a*c)*(b*c).

These identities give unique local coloring transport under Reidemeister I, II and III, so coloring count is invariant. One positive two-strand braid crossing acts on colors as B(a,b)=(2a−b,a), with a stated choice of which strand crosses over. Closure colorings are fixed points of B^m. Over F_3, B has exactly three fixed points and B^3 is the identity, so B^3 has nine. Thus the same endpoint matching and same component count do not determine knot type. The topological move interpretation is classical; the algebraic identities and fixed-point counts are reproduced here and checked independently.

There is therefore a valid loop-count bridge between path matching and tangles, but no justified assignment of a knot to an unembedded tournament or to the scalar 1729. A knot-sensitive carrier must retain an appropriate diagram/tangle class; a framed target additionally needs framing. No claim that Fox coloring itself is a complete carrier is made.

## D. The recovered 1729 coordinate has a precise observer

Let P11 be the Paley tournament on Z/11 with arc x→y when y−x is a nonzero square. Its square-affine group acts transitively on directed arcs. The corrected prior endpoint theorem yields, and the present two independent counts confirm,

   E=628 A+1101 A^t,  H=55(628+1101)=95095.

Thus 1729 is the sum of the two endpoint channel multiplicities and the number of free square-affine Hamilton-word orbits. It is not the count of words using one specified arc. Every Hamilton word uses ten arcs. Arc transitivity gives equal incidence at all 55 arcs, and double counting gives

   words containing a fixed arc = 10H/55 = 17290.

The checker also counts the specified arc 0→1 directly, so this value is not verified only through the same double count.

Closing a path requires its end to beat its start. The reverse endpoint channel therefore controls directed cycles:

   C=5·1101=5505,

with rotation identified and reversal not identified. This is a genuine attachment predicate that needs endpoint compatibility. The sum 1729 alone loses the split (628,1101). Session54 retains the stronger quadratic operator and its discriminant; this session does not reproduce its all-prime theorem or infer modular-form coefficients from numerical coincidences.

The underlying undirected graph of P11 is K11, with 55 edges, not a pseudoforest. Each individual Hamilton path support is a linear forest; confusing that support with the tournament is exactly the invalid transfer blocked in A. The equality 1729=7·13·19 creates no map to a Ramsey coloring or an embedded tangle.

## Reproduction and scope

`verify_structures.py` exhausts all three-colorings through five vertices and checks the sparse inequalities on those with no red/blue triangle. It verifies explicit extremal families for k=2,...,7, compares two loop algorithms on all matching pairs through eight ports, checks loop polynomials through twelve ports, checks Fox identities, and reconstructs the Paley-11 endpoint and specified-arc counts. `independent_checks.py` uses different graph and path algorithms and imports no producer code. Both run with explicit exception gates under normal Python and Python -O. The all-order proofs above do not rely on extrapolating these finite checks.
