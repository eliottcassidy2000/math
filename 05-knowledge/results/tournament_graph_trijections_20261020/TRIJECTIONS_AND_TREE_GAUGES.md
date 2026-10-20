# Three graph shapes, three exact correspondences, and the tree-gauge obstruction

2026-10-20. **PROVED, elementary:** the observer-coordinate formulas; general constrained-slice multiplicity/codimension theorem; tree-chart transition and complement-transport formulas; literal order-four chart; directed-agreement incidence correspondence; weighted cut-Cayley aggregate, its exact component/bipartiteness boundary, and generator injectivity from order eight. **INHERITED PROVED:** the spanning-tree transversal, Hamiltonian-path orbit fibre, switching-class path sum, and score-fibre Gram construction. **FINITE-EXACT:** the explicitly bounded controls in `verification.json`. No new natural unlabelled tournament/all-graph bijection, graph-colouring transfer, or original/black/regular-return metagraph-connectivity result is asserted.

## 0. Inheritance, scope, and the live board

The source question concerns tournaments as intransitivity/invariance and a “trijection” connecting tournament bijections with graph bijections. It singles out the three four-vertex, three-edge shapes P4, K1,3, and K3 plus an isolate, with a fixed Hamiltonian path associated to P4.

The closest proved mechanisms are:

- [THM-1430, tiling/class/metagraph dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md), Sections 1–3: the fixed-path cube is a spanning-tree transversal for labelled switching classes; its tournament-isomorphism fibres have size H/|Aut|. That file expressly credits the earlier [THM-1415, canonical star quotient](../../../01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md). Neither identification is claimed new here.
- [THM-1467, oriented-spanning-tree switching sum](../../../01-canon/theorems/THM-1467-oriented-spanning-tree-switching-sum.md): every oriented spanning tree is realized exactly once inside each labelled switching orbit; summing H over that orbit gives n!.
- [Defect-boundary web](../defect_boundary_20261020/web.md), Sections 4–6: the 12/4/4 graph census, complement versus tetrahedral opposition, cut/circuit distinction, the binary bicycle-space warning, and the affine observer correction. Those facts are inherited, not discoveries of this packet.
- [Tournament future boundary](../tournament_future_boundary_20261020/README.md), Section 6: the same three shapes carry different kernels for the characteristic-zero **signless** edge-sum map. That is not the binary cut-map rank used below.
- [Weighted SCC/Lyndon codec](../tournament_even_bijection_20261011c.md), Sections 1–4: a specified rank-dependent bijection of tournament **types** with **automorphism-even** graph types. Its graph class is not all simple graphs and “even” is not Eulerian. The [converse-context packet](../tournament_converse_context_20261011d.md), Sections 2–4, has a separate leaf-pair-marked star/P3-plus-isolates operation. Neither map is the observer chart or tree-gauge transition below.
- [Regular-return Gram addendum](../regular_return_gram_icosahedra_20261020/gram_return_addendum.md): literal common-intermediate incidence, raw H on the Gram diagonal, and the two labelled icosahedra. Section 7 below gives an actual graph-coordinate encoding of that incidence rather than identifying objects by their counts.

Canonical hostile: a directed triangle and a transitive three-tournament are switching-equivalent but have H=3 and H=1. Corrected near miss: all three eight-element graph powersets admit abstract bijections, but not all are transversals for the same cut quotient. Least-used sidecar: the prescribed orientation on a constraint graph, including its cycle-parity syndrome and relative component switch.

The live board is observer order / cut quotient / graphic-matroid rank / Hamiltonian marking / operation conjugacy / Gram incidence. The inherited META-PATTERNS cards used are “Type every analogy and every implication,” “Retain local incidence through counting,” and “Compute the repair quotient before testing the residual defect”; no new general card is claimed from this single lane. Anchor: an exact trijection. Niche: deficient constraint graphs. Wildcard: whether the graph encoding really reaches the existing icosahedral incidence. The answer is yes for incidence encoding, no for identifying the 20 three-edge graphs with the 20 faces of an icosahedron.

## 1. First correspondence: every observer gives a lossless tournament/graph chart

Fix a labelled vertex set V of size n>=2 and a total order pi on V. Let E be the unordered pairs, N=binom(n,2), and let R_pi be the transitive tournament pointing from earlier to later in pi. For a tournament T define the simple graph

    G_pi(T) = {uv in E : T disagrees with R_pi on uv}.

Every subset of E is allowed, and its inverse reverses precisely those arcs of R_pi. Thus

    labelled tournaments on V <-> simple graphs on V

is an exact bijection **after pi is specified**. It preserves every labelled tournament datum when the observer is retained, because it is invertible. It does not turn ordinary graph invariants into tournament invariants.

If rho is another observer, let I(pi,rho) be the graph of pairs on which the orders disagree. Then

    G_rho(T) = G_pi(T) xor I(pi,rho),                    (1)
    I(pi,sigma) = I(pi,rho) xor I(rho,sigma).            (2)

These are pairwise disagreement identities. Changing the observer changes no tournament arc; it translates its graph coordinates. Conversely, translating a graph while keeping the observer fixed **does** reverse tournament arcs. These operations must not be conflated.

Several exact commuting operation laws follow:

- Tournament converse T^op corresponds to ordinary graph complement: G_pi(T^op)=E\G_pi(T).
- Reversing a specified edge set D corresponds to symmetric difference with D.
- Vertex switching on a set S corresponds to symmetric difference with the cut delta(S).
- Simultaneous relabelling obeys G_(g pi)(g T)=g G_pi(T).

With the observer held at the reference order 1<...<n, the last action is **affine**, not just graph relabelling:

    G(g T) = g G(T) xor I(reference,g reference).        (3)

Its correction is the inversion graph of the relabelled observer. In particular

    |G(g T)| = |G(T)| + parity(g) mod 2.                 (4)

This is the orientation-parity transformation underlying the two sign sheets of the regular-return construction.

Any graph bijection f transports uniquely to the tournament bijection G_pi^{-1} f G_pi. That is an exact commuting square, but by itself is only conjugation. A substantive claim must still name the property or operation preserved by f.

### Cyclicity is not a graph shape

For i<j<k write a=x_ij, b=x_jk, c=x_ik for the backward-edge bits. The triple is a directed cycle exactly when a=b!=c. Its integer indicator is

    c + ab - ac - bc.                                  (5)

Summing (5) over triples recovers the tournament's cyclic-triangle count. Ordinary graph triangle count, edge count, or degree sequence does not do so. For the 20 backward graphs with exactly three edges at n=4, the finite exact H distributions are:

| Backward-graph shape | Labelled graphs with H=1 | H=3 | H=5 |
|---|---:|---:|---:|
| P4 | 2 | 4 | 6 |
| K1,3 | 2 | 0 | 2 |
| K3 plus isolate | 2 | 0 | 2 |

Even inside one of the three shapes, transitivity can change. The observer is essential, not cosmetic.

### No natural all-graph quotient is being asserted

There is no S_n-equivariant set bijection between all labelled tournaments and all labelled simple graphs for n>=2: the empty graph is fixed by every permutation, whereas no tournament is fixed by a transposition exchanging the endpoints of one of its arcs. At n=4 there are four tournament isomorphism classes but eleven graph isomorphism classes. This obstruction does not contradict the separately scoped automorphism-even graph-type codec cited above.

## 2. Exact normalization theorem for an arbitrary constraint graph

Work over F2 in X=F2^E, using the fixed observer's backward bits. The cut space is

    B = { delta(s) : s in F2^V },
    delta(s)_uv = s_u+s_v.

The kernel of delta is the constant vectors, so dim B=n-1. A labelled switching class is a coset x+B. Do not further quotient these classes by vertex relabelling; that would be a different set.

Let F be an arbitrary spanning constraint graph on V, including its isolated vertices, with m edges and c connected components. Fix the desired bit assignment a on E(F). Define the slice

    A(F,a) = {x in X : x restricted to E(F) = a}.

Its unconstrained coordinates are the subgraphs of the **free host** E\E(F), so it contains 2^(N-m) tournaments.

**Theorem 1 (constraint rank, cycle obstruction, and exact fibres).** Put beta=m-n+c, the cycle rank of F.

1. A switching class x+B meets A(F,a) if and only if x+a has even edge sum on every cycle of F. It suffices to check a cycle basis.
2. Every nonempty intersection has exactly 2^(c-1) distinct tournaments.
3. Exactly 2^(N-m-c+1) labelled switching classes meet the slice, a fraction 2^(-beta) of all 2^(N-n+1) classes.
4. The slice projects bijectively onto all labelled switching classes **if and only if F is a spanning tree**.

**Proof.** Required switches solve

    s_u+s_v = x_uv+a_uv,           uv in E(F).           (6)

Summing around a cycle proves necessity. If all cycle sums vanish, choose one root in each component, assign its potential arbitrarily, and propagate along a spanning forest. The cycle conditions make all omitted-edge equations hold. There are exactly two root-potential choices per component, hence 2^c potential vectors. Two vectors produce the same cut of K_n precisely when they differ globally by the constant one vector, giving 2^(c-1) distinct cuts and distinct tournaments. This proves (1)–(2).

Partition the slice, of size 2^(N-m), by its equally sized nonempty orbit intersections. Division gives (3). For a bijection one needs c=1 and beta=0, equivalent to connectedness and m=n-1; this is exactly a spanning tree. Conversely a spanning tree satisfies both. QED.

The binary graphic-matroid rank is n-c, not the number of prescribed edges. Adding a cycle creates a parity obstruction; disconnecting the constraint graph creates a switch ambiguity. The equality m=n-1 alone detects neither. This is the strongest survivor of the proposed common trijection: tree charts work, and the precise failure is measured by (beta,c-1).

This theorem is an elementary extension/repackaging of the inherited spanning-tree transversal, not a novelty claim for cut/cycle duality.

### Sharp four-vertex hostile

Each of the three free hosts has three edges and hence eight subgraphs. But their **constraint complements** behave differently:

| Free host | Prescribed constraint graph F | c | beta | Switching-class coverage |
|---|---|---:|---:|---|
| P4 | P4 | 1 | 0 | All 8 classes, once each |
| K3 plus isolate | K1,3 | 1 | 0 | All 8 classes, once each |
| K1,3 | K3 plus isolate | 2 | 1 | Exactly 4 classes, twice each |

For example, prescribe the triangle 23,24,34 to be forward and leave the star at vertex 1 free. A switching orbit is admissible exactly when

    x_23+x_24+x_34 = 0 mod 2.                           (7)

The opposite parity gives no representative. In an admissible orbit, switching vertex 1 toggles all three free star edges but preserves every triangle constraint. The two representatives therefore cannot be distinguished by the switching quotient.

No choice of a fixed orientation a on that triangle repairs this: it merely selects the other affine parity sector. A reversible chart can instead replace one triangle constraint by an edge to the isolate, obtaining a spanning tree. Alternatively it must retain both the cycle syndrome and the relative-component gauge data instead of pretending neither exists.

This binary rank distinction does not contradict the earlier characteristic-zero **signless** incidence calculation, where all three shapes had rank three. Over F2 plus and minus coincide; the triangle then has rank two. The field and the operator are load-bearing.

## 3. Second correspondence: switching classes and two graph hosts

Let F be a spanning tree. There is a unique cut which makes all its backward bits zero. Call the resulting representative N_F(x). Choose a root r, set s_r=0, and propagate

    s_v = sum of x_e over the unique F-path from r to v.

Then N_F(x)=x+delta(s). It depends only on x+B. For an edge e=uv outside F, let C_e^F be the fundamental cycle formed by e and the unique u–v path in F. Its free coordinate is

    N_F(x)_e = sum_(f in C_e^F) x_f.                   (8)

Thus the graph on E\F consists exactly of fundamental-cycle parities. This gives a complete labelled switching invariant, not a complete tournament-isomorphism invariant and not a cyclic-triangle indicator.

For two spanning trees F,H, define the graph-chart transition

    L_(F->H)(z) = N_H(z),            z restricted to F=0. (9)

It is an invertible linear map between their free graph-coordinate spaces, with inverse L_(H->F). It obeys the exact cocycle law

    L_(H->J) L_(F->H) = L_(F->J).                       (10)

Indeed successive normalization remains in the same switching coset, and the J-normalized representative is unique. Consequently there is a genuine three-way bijection

    labelled switching classes
       <-> subgraphs of E\F
       <-> subgraphs of E\H.                          (11)

Its invariant is the switching coset. Its loss relative to T is the cut representative: n-1 bits must be retained to reconstruct an arbitrary original labelled T. Its dependency is the two chosen trees and desired orientations. No arbitrary canonical order is smuggled into an unlabelled claim.

### Literal P4-to-triangle map at order four

Use path constraints 12,23,34 and star constraints 12,13,14, all forward. The path's free host is the P4 with edges 13,14,24. Write their backward bits as (a,b,c). Choose switching potentials

    (s_1,s_2,s_3,s_4)=(0,0,a,b).

The star's free host is the triangle on 2,3,4 plus isolate 1. Its bits on 23,24,34 become

    (u,v,w) = (a,b+c,a+b).                             (12)

The inverse is

    (a,b,c) = (u,u+w,u+v+w).                           (13)

All sums in (12)–(13) are in F2. This is the requested concrete graph-to-graph bijection compatible with the tournament switching correspondence. It is not merely a matching of eight objects by rank.

The third eight-state host, the star, can of course be matched abstractly to either graph powerset. But Theorem 1 proves that such a match cannot simultaneously be the same switching-quotient chart with its triangle complement prescribed. Its natural map has four double fibres and four missed classes.

### What is lost: even H is not invariant

Already at n=3, switch the middle vertex of the order 1->2->3, with 1->3. The new tournament is 2->1->3->2, so H changes from 1 to 3 in the same switching class.

There is a literal hostile inside (12). For (a,b,c)=(0,1,0), the path-normalized tournament has the four-cycle 1->2->3->4->1, with 1->3 and 2->4, and H=5. Its star-normalized image is (u,v,w)=(0,1,1), namely the transitive order 1->4->2->3, with H=1. The bijection preserves the cut class and changes H, strong connectivity, scores, and cyclicity. These are not promised invariants.

## 4. Which operation really commutes with the graph bijection?

On a fixed directed Hamiltonian path P, let K_P reverse every arc outside P and preserve the path. In the P-normalized chart it is exactly graph complement **inside the free host** E\P:

    z -> z + 1_(E\P).                                (14)

This is not full tournament converse, which also reverses the path arcs.

Under a change to tree H, equation (8) gives the transported operation

    y -> y + d_(P,H),
    d_(P,H)_e = |C_e^H \ E(P)| mod 2,    e outside H.  (15)

**Proof.** N_H is linear, so N_H(z+1_(E\P))=N_H(z)+N_H(1_(E\P)). Apply (8) to the final term. QED.

For the literal order-four map (12), complementing the three path-free bits changes (a,b,c) to (a+1,b+1,c+1), and therefore

    (u,v,w) -> (u+1,v,w).                              (16)

So **path complement becomes a single triangle-edge flip**, not complement of the triangle host. Equivalently L(1,1,1)=(1,0,0). This is an explicit, nontrivial commuting operation law. Demanding “complement on both graph hosts” would be a false diagram. Full tournament converse, path complement, observer reversal, and chart change are four different maps.

More generally every fixed flip set D transports by the fundamental-cycle syndrome N_H(D). Only the chosen coordinate chart determines whether that syndrome is one edge, many edges, or zero.

## 5. The fixed Hamiltonian fibre is not the 20-graph layer

Fix the directed path 1->2->...->n. Its n-1 arcs are prescribed and the other N-(n-1)=binom(n-1,2) arcs are free. Therefore its fibre contains exactly

    2^binom(n-1,2) tournaments.                        (17)

At n=4 this is eight, whereas the middle layer of all graph coordinates has binom(6,3)=20 elements. The latter is the family of all three-edge subgraphs of K4, not the family of completions of one path.

The off-path host happens to be another P4 only at n=4, apart from the trivial n=1 case. Equality of its edge count with that of the path requires

    binom(n,2)=2(n-1), or (n-1)(n-4)=0.

At n=4 the complement of a P4 is indeed a P4; the two paths partition the six edges. This exception does not make all three three-edge shapes trees. There are 16 spanning trees of K4: 12 paths and four stars. A triangle plus an isolate is dependent and has graphic rank two.

The inherited Hamiltonian incidence double count is

    sum_(labelled T) H(T)=n! 2^binom(n-1,2),
    sum_(tournament types [T]) H(T)/|Aut(T)|=2^binom(n-1,2). (18)

The second equation counts path-marked tournaments modulo simultaneous relabelling. Normalizing the marked Hamiltonian word gives the fixed-path cube; the action on Hamiltonian words is free. The same cube is also a cut transversal by Theorem 1, but those are **different projections**. Ordinary tournament isomorphism, cut switching, and forgetting a Hamiltonian mark are not interchangeable quotients.

## 6. Third correspondence: path witnesses and directed agreement graphs

For labelled tournaments T,S let A(T,S) be their undirected agreement graph: a pair is an edge when both tournaments orient it the same way. In observer coordinates,

    E(A(T,S)) = E \ (G_pi(T) xor G_pi(S)).              (19)

**Theorem 2 (literal path-complement incidence).** The following objects are in bijection:

1. A tournament T with a marked directed Hamiltonian path P.
2. A triple (T,P,S) with S=K_P(T).
3. An ordered pair (T,S) whose agreement graph is a spanning path directed consistently from one end to the other in T.

The map sends (T,P) to (T,K_P(T)); its inverse recovers P as the agreement graph with its directions in T. Such directions determine a unique ordered path. Since the chosen path is also directed in S, the incidence is symmetric under exchanging T,S. All labels and directions are retained.

This is an actual three-way incidence, not an unexplained trijection among unlabelled species. Its graph predicate contains both the underlying path and its orientation. A star agreement graph and a triangle-plus-isolate agreement graph cannot qualify; nor does a P4 with an inconsistently oriented middle vertex qualify.

At n=4, agreement is P4 if and only if disagreement is P4. For each of the 12 possible disagreement paths, its complementary path has eight possible orientations, exactly two of which are consistently directed. There are eight choices for the remaining arcs. Thus each disagreement path supports 2*8=16 ordered incidences, for a total

    12*16 = 192 = 4!*8.                               (20)

Neither the four star supports nor the four triangle-plus-isolate supports yields a path-complement incidence. This supplies an operation-based reason to select P4 from the three-shape layer rather than a numerical resemblance.

## 7. Actual connection to the existing Gram/icosahedral construction

Fix a labelled score vector d, let F_d be its tournament fibre, and use the inherited binary incidence matrix B[T,S] for S=K_P(T). Equations (19) and Theorem 2 re-express each nonzero entry as:

- the symmetric difference of two observer graphs is the complement of a spanning path;
- the recovered path is directed in T;
- T has the required labelled score vector.

This is lossless if the observer, labels, and direction decoder are kept. Forgetting them is not a map to ordinary unlabelled graph shapes that preserves B.

For a path from s to t, direct degree accounting gives

    d(K_P(T))=(n-1)1-d(T)+e_s-e_t.                     (21)

Two path complements return to the same score vector precisely when the two paths have the same ordered endpoints. Their net edge reversal set is

    E(P) xor E(Q),                                    (22)

so it has even size. Under (1), observer changes cancel in the XOR of two tournament encodings, which is why the agreement/disagreement object is observer-independent as a graph on the fixed vertices.

The inherited unique-witness Gram identity therefore remains literally

    W=BB^T,     W[T,T]=H(T).                           (23)

No extra multiplicity or unidentified path has entered the graph encoding. At regular order five, the previously proved face columns and private columns give

    W=10I+FF^T=15I+2A

on each of the two labelled icosahedral components. The map in (19) preserves the actual incidence and hence its Gram matrix. It does **not** identify the 20 three-edge graphs on four vertices with the 20 triangular faces of either icosahedron. Those have different defining sets, labels, and incidence; an equal cardinality alone supplies no bridge.

## 8. A structural payoff: the cut-partition aggregate is a Cayley graph

Let A_n be the graph on **all labelled tournaments**, joining T to K_P(T) for every directed Hamiltonian path P. For n>=3 its edges have unique path witnesses by Theorem 2, and there are no loops. Partition its vertices into labelled switching classes q in Q=X/B. Define

    C[q,r] = #{ ordered incidences (T,S) : T in q, S in r, S=K_P(T) }.

This is an **aggregate incidence quotient**, not an equitable quotient: individual vertices in one switching class can have different H and different outgoing neighbour counts. The cut action need not be an automorphism of A_n. The aggregate nevertheless has a clean exact structure.

**Theorem 3 (Cayley aggregate and all-order connectivity boundary).** For n>=3 and every undirected Hamiltonian path P, put g_P=[E\P] in Q.

1. The aggregate is the weighted Cayley matrix

       C[q,r] = 2 #{undirected Hamiltonian paths P : r=q+g_P}. (24)

   Every row sum is n!, no generator is zero, and the support graph is regular. All positive weights are equal at a fixed n.
2. Its support graph is connected unless n=1 mod 4. In that exceptional case it has exactly two connected components, distinguished by edge-weight parity of any tournament coordinate representative.
3. It is connected and bipartite when n=3 mod 4. For even n it is connected and nonbipartite. When n=1 mod 4, each of its two components is nonbipartite.
4. For n>=8, distinct undirected Hamiltonian paths give distinct generators. Consequently the support degree is exactly n!/2 and every aggregate edge weight is 2.

**Proof of (1).** Fix a switching class q and an ordered Hamiltonian vertex word p. Its edges form a tree with a prescribed orientation. By Theorem 1 exactly one member T of q realizes p. The image K_p(T) has class q+[E\P], independent of that representative. Summing over all n! words gives (24), with the two orientations of every undirected P counted separately. Distinct words cannot produce the same ordered tournament pair by Theorem 2. Thus the sum is an incidence count, not an overcount by hidden witnesses.

If E\P were a nontrivial cut, then P would be the disjoint union of two complete graphs, contrary to connectedness. The trivial cut would require P=K_n, impossible for n>=3. So there are no zero generators. Relabellings act linearly on the displacement quotient X/B and transitively on undirected Hamiltonian paths. Therefore every generator fibre has the same size, which proves equal positive weights. The state-coordinate action on tournaments is still affine as in (3); no contradiction is involved.

**Proof of (2)–(3).** Let X_even be the even-edge-weight subspace of X. The span of all differences of Hamiltonian path edge vectors is precisely X_even. Inclusion is immediate because each path has n-1 edges. For the reverse inclusion, compare

    (a,b,...,c) and (b,...,c,a).

Their edge-set difference is {ab,ac}. Every pair of adjacent edges arises this way. Such pairs span the entire even-weight subspace because the line graph of K_n is connected: pair any two edges by telescoping along an adjacent-edge chain. Since g_P+g_Q=[P+Q], the generator-difference span is exactly the image of X_even in Q.

For even n, a one-vertex cut has the odd size n-1. Hence every coset has an even representative, the image of X_even is all Q, and the generators span Q. For odd n, every cut has even size |S|(n-|S|). Thus parity descends to Q and the image of X_even has codimension one. Every generator has parity

    |E\P| = binom(n-1,2) = (n-1)/2 mod 2.              (25)

If n=3 mod 4 it is odd, so adding any generator to the difference span gives all Q. If n=1 mod 4 it is even, so the span is exactly the even-parity hyperplane. Components of a Cayley graph on an F2-vector space are the cosets of its generator span. This proves (2).

For bipartiteness, a Cayley graph on its generator span is bipartite exactly when some linear functional takes value one on every generator. One direction colors by that functional. Conversely, if no odd closed walk exists, the parity of the number of generators in a representation of a vertex is well-defined and gives such a functional. A functional constant on all generators vanishes on their differences. For even n, or within either n=1 mod 4 component, those differences already span the component, making a constant-one functional impossible. For n=3 mod 4 the parity functional supplies it. This proves (3).

**Proof of (4).** Equality g_P=g_Q means P xor Q is a cut. Its maximum vertex degree is at most four. If n>=9, a nontrivial cut has a side of size at least five, and every vertex in the other side has degree at least five, impossible. For n=8 the only remaining cut shape is K4,4, with every degree four. But an endpoint of P has degree at most 1+2=3 in P xor Q. Thus the cut is zero and P=Q. Each undirected path has two directed orientations, giving degree n!/2 and weight two. QED.

The threshold in (4) is sharp: the exact order-seven census has 2,520 undirected paths but only 1,260 generators. This finite boundary is supplementary; the theorem for all n>=8 is proved by the degree argument.

### The three small shapes determine the complete order-four quotient

At n=4, all four stars are cuts, so they lie in class 0. All four triangles-plus-isolate are complements of stars, so they lie in class [E]. Every P4 lies in one of the other six classes, two paths in each. Here is a short multiplicity proof, not an inference from equal counts. For P=12,23,34, the only distinct three-edge graph in its switching class comes from the balanced cut between {1,4} and {2,3}, giving P'=13,23,24. The other two balanced cuts intersect P in one or three edges and produce five or one edges, respectively; odd cuts cannot preserve odd weight. Relabelling proves the same two-element fibre for every P4.

Complementing P adds [E], permuting those six nonzero, non-[E] classes. Hence the cut aggregate has every possible adjacency except equality and the antipodal pairing q<->q+[E]:

    support = K8 minus a perfect matching = K_(2,2,2,2),
    every aggregate edge weight = 4.                  (26)

This is the graph of the four-dimensional cross-polytope, but only the explicitly proved adjacency is used. Its eight vertices are labelled switching classes, not the 20 three-edge graphs and not the 24 regular five-tournaments.

Finite exact generator data are:

| n | dim Q | distinct generators / degree | positive aggregate weight | component count |
|---|---:|---:|---:|---:|
| 3 | 1 | 1 | 6 | 1 |
| 4 | 3 | 6 | 4 | 1 |
| 5 | 6 | 15 | 8 | 2 |
| 6 | 10 | 180 | 4 | 1 |
| 7 | 15 | 1260 | 4 | 1 |
| 8 | 21 | 20160 | 2 | 1 |

At order five its two components each have 32 vertices. They are **not** the two 12-vertex regular-return icosahedra. The node sets, allowed steps, and retained scores differ.

The connectivity theorem is a real payoff of retaining the right quotient and the path-incidence multiplicity. It does not lift to the original path-complement graph from a specified source: consecutive quotient edges may require different representatives of the intermediate switching class. Switching between them is an additional operation that the original path process need not allow. It therefore gives no black-metagraph or regular-score-fibre connectivity theorem.

## 9. What a valid “trijection” can mean

There are three clean answers, with different preserved structures:

1. **Observer conjugacy:** T_n <-> G_n, transporting a named graph bijection into a tournament bijection. Keeps every bit plus the observer; no ordinary graph invariant is automatically intrinsic.
2. **Tree-gauge trijection:** labelled cut classes <-> P4-host graph states <-> triangle-host graph states at order four, by (12). Keeps switching data; loses the representative and H. The star-host third chart is obstructed exactly by Theorem 1.
3. **Hamiltonian incidence trijection:** marked tournament/path <-> path-complement triple <-> directed-agreement pair, by Theorem 2. Keeps the actual witness; it is the route into the inherited Gram construction.

The SCC/Lyndon automorphism-even graph-type codec is a fourth, different construction. It can also transport operations, but its atom ranking, factorization context, and canonicalization data must not be replaced by these tree charts. In particular its literal star/P3-plus-isolates toggle concerns P3, not the three-edge P4 here.

## 10. Reproducibility and stopping boundary

Run from this directory:

    python3 verify_trijections.py > verification.json
    python3 -O verify_trijections.py > verification.optimized.json
    cmp verification.json verification.optimized.json

The verifier uses standard-library enumeration and explicit RuntimeError checks, never optimization-disabled assertions. Both modes agree byte-for-byte. It checks every labelled tournament at orders 2–5; every constraint graph and every prescribed pattern through order four; every order-five constraint graph with patterns 0, all-one, and one fixed alternating mask (duplicates removed); every spanning-tree normalization through order five; every ordered tree-chart pair through order four, and transitions from every order-five tree to the fixed path and star. It also checks observer affine actions, cyclic-triangle formulas, path-witness uniqueness, switching-class H sums, all 192 order-four path-complement incidences, and the literal aggregate Cayley identity through order five. A separate direct fundamental-triangle encoder enumerates all undirected path generators through order eight, checking both their linear span and their difference span, the component/bipartiteness ranks, and the displayed multiplicities. The output records the finite universes and H histograms.

The all-order proofs do not rely on extrapolating those finite tests. The packet makes no claim to settle H-spectrum completeness, all-order regular-return connectivity, black-metagraph connectivity, or a natural unlabelled graph/tournament bijection. The new cut-aggregate connectivity theorem has its stated proof and scope; it cannot be substituted for those open targets. The useful endpoint is an exact dictionary, a sharp quantified obstruction, and a genuine all-order quotient structure, not a forced identification of all three shapes.
