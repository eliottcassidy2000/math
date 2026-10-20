# Session 16: path-count spectra, regular return geometry, and torus winding

20 October 2026. Proof-first continuation of the tournament metagraph program.

**Global status:** It remains open whether the ordinary Hamiltonian-path counts of finite tournaments are exactly the positive odd integers other than 7 and 21. The two omissions are established; no theorem below proves completeness or creates another global omission. Ordinary counts H, path-orbit counts H/Aut, regular-only spectra, and values at a fixed order are different objects.

This session supplies structural classifications and obstructions, not a larger brute-force search for the first global gap. Classical inputs and independently checked finite examples are distinguished in the linked packets. No external priority claim is made.

## 1. A complete intersecting-cycle subclass, and a limit on that strategy

Call a tournament **intercyclic** when it has no two vertex-disjoint directed cycles, equivalently no two disjoint directed triangles. Let S_1 be the all-order path-count spectrum of tournaments which become transitive after deleting at most one vertex. The complete subclass classification is

    S_intercyclic = S_1 ∪ {3·2^t+1 : t≥2} ∪ {9·2^t−3 : t≥1}.

The proof imports the two-vertex transversal theorem of Bang-Jensen–Bessy, then classifies the orientation words over a transitive remainder by cross-intersecting inversion families. It counts the remaining odd cycles directly. Both additional rays have explicit strongly connected witnesses.

Consequences:

- Every tournament with H=29 has triangle-packing number exactly two. An explicit strong six-vertex witness establishes attainability. The proof rules out the whole intercyclic class, not just a finite bank.
- Divisibility by 29 is insufficient: the strong one-apex word 1100100 has H=87=3·29, while all its cycles contain its apex.
- The simplest one-apex grammar cannot realize 13, although unrestricted tournaments can. Its sharp strong-order minima and maxima are proved by a strict adjacent-exchange formula.
- A general marked feedback set of size k has an exact finite-state path counter of dimension at most k!·3^(k+1), with no omitted insertion-orphan term.
- Fixing both k and the number R of orientation runs produces only O((log X)^R) different values up to X. Positivity, commuting idempotents, and a surjection expansion prove this obstruction without a large tournament census.

See [the feedback-word packet](feedback_word_hamilton_spectrum_20261020/README.md) and its exact intercyclic theorem. The full global problem still needs constructions beyond this intersecting-cycle subclass.

### Quantitative binary-product refinement

For a reduced one-apex binary word b, let N be its binary integer and R its reversed-complement binary integer. The exact path formula gives

    H(b) − NR/2^m = 1 + number_of_10_transitions − backward_weight,

with absolute error at most 1+m/2. Uniform reduced words therefore have normalized counts H/2^m converging to the product of two independent Uniform[1/2,1] variables. The argument uses independent finite prefixes at opposite ends, not a false assertion that the two finite integers are independent.

Summing the word sublevels across lengths cancels the dyadic phase exactly:

    #{reduced words b : H(b)≤X} = X log 2 + O(√X).

This counts presentations, not different values. Complement-reversal and the exact two-ray classification imply

    |S_intercyclic ∩ [1,X]| ≤ (log 2/2)X + O(√X).

Thus the lower relative asymptotic density of odd integers outside the intercyclic spectrum is at least 1−log 2, approximately 30.6853%. This is not a density statement about globally forbidden tournament counts. An attainable value outside this subclass requires two disjoint triangles. No existence of a numerical-spectrum density is asserted.

## 2. Two path-complements reveal exact polyhedral incidence

Preserving a directed Hamilton path P while reversing all other arcs changes the score vector to

    d(K_P(T)) = (n−1)1 − d(T) + e_start(P) − e_end(P).

For n>1, one such move exits regularity. Two successive moves return to the original labelled score vector exactly when their ordered endpoint pairs agree. Every nontrivial two-step return reverses at least six arcs; the six-arc case is classified by an interchange of two interior path blocks.

On five labelled vertices the resulting regular return graph is exactly **two icosahedra**, separated by orientation parity. Each has 12 vertices, 30 edges, and 20 triangular faces. All 24 labelled states belong to one ordinary tournament isomorphism class, so the geometry is not preserved by simply forgetting labels.

The faces are not a visual analogy: each is the neighbourhood of an actual intermediate tournament. Let B record source/intermediate incidence in one labelled score fibre. The path witness for an entry is unique, so the return-presentation matrix satisfies

    W = BBᵀ,    W(T,T)=H(T).

At order five, each parity component satisfies

    W = 10I + FFᵀ = 15I + 2A,

where F is vertex/triangular-face incidence and A is the icosahedral adjacency matrix. Each node has 15 identity-return presentations and two presentations to each of its five neighbours.

A precisely specified uniform-path/then-uniform-legal-continuation walk is reversible with measure H(T). On a fixed-order regular fibre it is lumpable under relabelling and converse; its ordinary and merged stationary masses are proportional to H/Aut and κH/Aut, respectively. These are exactly the established tiling-fibre masses. At order seven the regular ordinary-class weights are 57:25:9.

This does not make deterministic complement descend to isomorphism nodes, prove general irreducibility, or identify the derived graph with the original black/blue metagraph. The full endpoint-resolved incidence is essential.

See [regular return geometry and Gram kernels](regular_return_gram_icosahedra_20261020/README.md), including the labelled illustration and stationary-mass proof.

## 3. Four complementary-graph games

The game is a one-move complementary-graph problem. Player 1 chooses red R, and Player 2 searches its blue complement. A spanning walk permits repeated vertices and edges and is equivalent to connectivity. In variant A/B, Player 1 wins when red lacks property A and blue lacks property B.

With Player 2 free to choose a favourable starting vertex, the H/H, H/walk, and walk/H variants are possible exactly for n≥4. The walk/walk variant is impossible: the complement of a disconnected graph is connected. Their minimum red-edge costs are n−1, n−1, and 2n−5, respectively, except the last cost is 6 at n=6 and 10 at n=8.

Requiring Player 2 to succeed from every prescribed starting vertex changes the thresholds for defeating that guarantee in the two Hamilton-output variants to n≥3 and the corresponding minimum costs to n−2. The rooted path-cover criterion for graph joins characterizes these endpoint effects exactly.

Both colours can also be connected and nontraceable, first at n=7: the construction works at every larger order, and the exclusion through six has independent exact checks. The game is not a tournament problem unless an orientation map preserving the desired predicate is supplied; a naive red/backward orientation fails.

See [the complete game packet](complement_path_games_20261020/README.md).

## 4. A real torus bridge, with its missing coordinates retained

The triangular-lattice quotient (a,b)↦a+3b mod7 is a cellular toroidal K7. Its lift is three-colourable, although its quotient needs seven colours: periodic identifications can manufacture a chromatic obstruction. This is not a seven-vertex Euclidean unit-distance graph.

Two orientations give the regular carousel and Paley tournaments. Both have Euler-tour length 21, but H=175 and H=189. A winding triangle reversal yields the third regular class with H=171. Degree alone, Euler length, and total Euler flux do not determine H.

The exact transport is through the directed line graph:

    H(L(D)) = |E(D)| · E_cyclic(D)

for a finite strongly connected balanced digraph with at least one distinguishable arc. The line graph is generally not a tournament. Every Euler tour has the same total edge-cocycle homology; Hamilton cycles instead have a nontrivial winding polynomial. The relation to Hamilton paths is H=nC+R, with R counting nonclosing paths. Paley's R=21 is not an example with total H=21.

A checkerboard torus embedding forces even Euler homology, obstructing a primitive essential knot slope. Graph automorphisms need not preserve the marked ambient knot type: the Paley multiplier permutes slopes representing an unknot, T(2,7), and T(3,7).

All 2,640 labelled regular orientations of this particular toroidal K7 have 25 Euler-flux sectors, of sizes 660 once and 273,42,14,1 six times each. The structural face-flip connection is cited from Zhang–Qian with hypotheses checked; the exact sector counts have complete independent certificates.

The appearance of 273 is real but does not identify the natural move graphs: a 273-state face-flip sector has 770 edges; the 13·21-state two-path matching graph has 914. Their equal cardinality does not imply isomorphism.

See [Euler/torus transport](euler_torus_hamilton_transport_20261020/EULER_TORUS_BRIDGES.md) and [voltage-sensitive geometry](torus_voltage_and_halfplane_forcing_20261020/README.md).

## 5. The original Fibonacci identity and the geometric boundary

Incoming repository work already gives the matching interpretation of

    1,2,6,15,40,104,273,... = F_j F_(j+1) = Σ_(i=1)^j F_i²,

including the explicit 63+210 matching partition at 273 and the marked torus framing of T(13,21). These inherited results are not claimed anew. Marking and winding are necessary; the same abstract directed triangle can be embedded along arbitrarily many torus-knot types without changing its H=3 or Euler length three.

For a Toroflux-style mathematical model, treat the material wire as a closed embedded curve and its motion as an isotopy unless a cut or self-intersection is permitted. Such a motion preserves knot type; subdividing the same wire into more graph edges changes combinatorial Euler length without changing that type. A toroidal appearance alone does not supply the marked winding pair, and no elastic-motion model of the physical toy is claimed here.

Unconditionally, every open or closed half-plane has the full plane's unrestricted unit-distance chromatic number, by finite-configuration translation and compactness. This is not a measurable-colouring transfer. The recent external six-colour manuscript was checked for source definitions and provenance, not independently built or fully proof-verified; no new plane chromatic bound is asserted here.

Finally, a finite-boundary counting proof rules out duplicating the square lattice by finitely many translations, or the triangular lattice by its lattice-preserving isometries, even with arbitrary infinite pieces. This classic nonparadoxicality result has a precise action and scope. It neither extends Banach–Tarski to a lattice nor substitutes quotient identification for a piecewise isometric bijection.

## Remaining mathematical obligations

Global odd-spectrum completeness; a sufficiently rich constructive grammar beyond intersecting cycles; general regular-return connectivity; black-metagraph connectivity; and structural maps between the different 273-state families remain open. The next useful step must retain the endpoint, multiplicity, winding, or cycle-packing coordinate that the proved obstructions show cannot be discarded.
