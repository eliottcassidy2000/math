# Finite algebraic plane witnesses: critical extraction and the six/seven distinction

**Status.** **CONDITIONAL** on the user-accepted statement that the Euclidean plane has no proper five-colouring. The compactness and real-closed-field transfer are **CITED classical** and already recorded in [finite_algebraic_witness_20261007](finite_algebraic_witness_20261007.md). The critical-extraction proofs and certificate interface below are **PROVED elementary**, without a novelty claim for those ingredients. The implementation is **FINITE-EXACT** on the stated small universe. It supplies a four-chromatic Moser spindle control, **not a six-chromatic graph**, and does not audit the accepted OpenAI theorem.

## 1. Inheritance, types and the requested conclusion

Write $U(X)$ for the graph on $X\subseteq\mathbb R^2$ whose edges are exactly the pairs at Euclidean distance one. Such a graph need not be planar as an abstract graph: its straight unit edges may cross. The four-colour theorem for planar graphs therefore does not provide an upper bound here.

The [official OpenAI scope statement](https://github.com/openai/math/blob/main/lean/docs/158.md), accompanying [the supplied paper](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf), explicitly concerns arbitrary colour classes and states a seven-colouring with all boundary points included. In this note the lower bound is an accepted hypothesis $H:\chi(U(\mathbb R^2))\ge6$, not a newly verified theorem. Together with the classical seven-colouring it leaves the alternatives six and seven.

The closest proved repository mechanism is the earlier exact real-root/coordinate/refutation checker. The decisive hostile is a restricted coordinate field known to be five-colourable: searching more points in that field cannot yield six colours. [THM-4558 — residue colourings of field planes](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md) records the repaired field bounds. MISTAKE-576/579 replaced the false class-number roadmap; no such roadmap is used here. The least-used sidecars are the chosen real embedding and the original labels retained through deletion.

The live concepts are compactness, algebraic realization, critical deletion, finite proof checking, and exact geometric boundary ports. The incoming [planar_triangle_carrier](planar_triangle_carrier_20261011.md) concerns planar Jacobian/Hamiltonian maps and generic algebraic curves. Its use of “planar” supplies no graph-colouring transfer.

## 2. The classical finite-algebraic chain

Graph compactness says that a graph is $q$-colourable iff all its finite subgraphs are. The original reference is [de Bruijn–Erdős, 1951](https://research.tue.nl/en/publications/a-colour-problem-for-infinite-graphs-and-a-problem-in-the-theory-/). One proof assigns a variable with values in a finite palette to every vertex: finite edge constraints have a common solution precisely when the corresponding closed subsets of the compact product of finite palettes have the finite intersection property.

Under $H$, choose a finite point set whose induced unit graph is not five-colourable. Repeated vertex deletion gives an induced graph $G$ still not five-colourable, with $G-v$ five-colourable for every vertex $v$. Adding a new colour for $v$ proves $\chi(G)=6$.

For this fixed labelled graph, let $d_{ij}=(x_i-x_j)^2+(y_i-y_j)^2$. The finite existential formula

$$
\bigwedge_{i<j}d_{ij}>0,\qquad
\bigwedge_{ij\in E(G)}d_{ij}=1,\qquad
\bigwedge_{ij\notin E(G)}d_{ij}\ne1
$$

has a real solution. Real-closed-field transfer gives a solution in the real algebraic numbers $\mathbb R_{\rm alg}$. Distinctness and nonedges are explicit; hence the whole induced graph, not just its required edges, survives. Finitely many algebraic coordinates lie in one finite real number field. This recovers a finite algebraic graph of chromatic number exactly six.

[Madore, §5.4](https://arxiv.org/pdf/1509.07023) gives the classical algebraic-field and effective-enumeration consequences. ZFC supplies the compactness framework; algebraic coordinates additionally use real-closed-field transfer. Neither premise provides a practical coordinate list or numerical search deadline.

### A complete search in principle

Enumerate finite simple graphs. For each graph, decide five-colourability by exhaustive finite search. For each non-five-colourable graph, decide the displayed existential real formula by an effective real-closed-field procedure, obtaining algebraic sample coordinates on success. Every individual test terminates. Under $H$, some candidate succeeds.

An alternative is to enumerate finite descriptions consisting of a squarefree rational polynomial, a rational interval isolating one real root, and rational coordinate polynomials in that root. Enumerating these finite strings fairly is complete: a finite tuple of real algebraic numbers lies in a simple finite extension in characteristic zero. Exact arithmetic then reconstructs every distance-one incidence. Invalid descriptions are rejected by finite tests.

The present program implements the exact candidate checker, critical extraction, and a consumer of a supplied candidate stream. It does **not** implement a general real-closed-field solver or claim its finite demonstration stream is complete. Conditional halting of the full search is a mathematical conclusion; its implementation and a practical size bound are separate tasks.

## 3. A checked critical certificate

Here **vertex-critical at $k$** means $\chi(G)=k$ and every proper induced subgraph is $(k-1)$-colourable. It does not mean that deleting every edge reduces the chromatic number.

A packet contains:

1. one polynomial, its isolating interval, and the coordinate polynomials;
2. the complete induced unit graph on those distinct points;
3. a branching refutation of $(k-1)$-colourability;
4. a proper $k$-colouring;
5. a proper $(k-1)$-colouring for each marked one-vertex deletion.

The checker reconstructs item 2 from item 1. The inherited refutation checker verifies all palette branches and checks each conflict leaf against an actual monochromatic edge. It separately checks every supplied colouring. A proper induced subgraph omits some vertex $v$ and inherits the checked colouring of $G-v$, so the one-vertex witnesses certify the full induced-minimality condition.

The extraction algorithm receives a finite non-$(k-1)$-colourable graph. In increasing original-label order, it deletes a vertex exactly when the remaining graph still rejects $(k-1)$ colours. One pass suffices: if deletion of a retained vertex was colourable at its test, any later smaller deletion graph stays colourable by restriction. At termination, non-colourability plus one-deletion colourings proves the chromatic number is exactly $k$.

The output keeps the original vertex-index tuple. This is the map back to the source geometry, and prevents abstract isomorphism alone from silently identifying different geometric ports. Exact embedding checks are repeated after extraction. An edge-minimal non-$(k-1)$-colourable **subgraph** can additionally be obtained by deleting edges; its coordinates remain algebraic, but it need not be the induced unit graph on those vertices. The supplied packet certifies the induced vertex-critical notion only.

### SAT-compatible certificate object

For $N$ vertices and $q$ colours, variable $X_{v,c}$ means vertex $v$ has colour $c$. Include one at-least-one clause per vertex, the $\binom q2$ at-most-one clauses per vertex, and $q$ clauses per edge forbidding equal colours. Satisfying assignments correspond exactly to proper colourings.

The program returns this CNF as an explicit alternative input to a SAT search. Its verification control uses the inherited branching refutation, not an unverified solver response. A larger implementation may use a compressed independently checkable SAT refutation; the source coordinates, complete edge list, and proof all remain necessary.

## 4. What changes between six and seven

Let $C=\chi(U(\mathbb R^2))$ and assume $H$ and the seven-colour upper bound.

| Case | Finite algebraic consequences | Global consequence |
|---|---|---|
| $C=6$ | A vertex-critical six-chromatic witness exists; no finite unit-distance graph requires seven colours. | A proper six-colouring exists. No particular algorithm for it is furnished here. |
| $C=7$ | A vertex-critical seven-chromatic witness exists, and contains an induced six-chromatic vertex-critical witness after further deletion. | No proper six-colouring exists. |

For the second cascade, delete any vertex from a vertex-critical seven-chromatic graph. Its chromatic number is six: it is at most six by criticality and cannot be at most five, since restoring the vertex with a new colour would give six colours for the original graph. Extract a vertex-critical six-chromatic induced subgraph from it. Every coordinate remains in the same field.

Thus a finite six-chromatic witness confirms the accepted lower bound but cannot distinguish $C=6$ from $C=7$. Searching for a finite non-six-colourable unit graph halts exactly in the seven case. A finite unsuccessful search is not an upper-six proof. Compactness is an existence principle for global colourings, not an algorithm supplied by this packet; no claim is made that every global six-colouring must be noncomputable.

### Some finite real number field already attains the full value

Enumerate all real algebraic numbers as $\alpha_1,\alpha_2,\ldots$ and put
$K_j=\mathbb Q(\alpha_1,\ldots,\alpha_j)$. These are nested finite real extensions with union $\mathbb R_{\rm alg}$. Compactness and the algebraic transfer imply

$$
\chi(U(\mathbb R^2))
=\chi(U(\mathbb R_{\rm alg}^2))
=\sup_j\chi(U(K_j^2)).
$$

Indeed every finite witness has only finitely many coordinates, so lies in one $K_j^2$. The bounded nondecreasing integer sequence $\chi(U(K_j^2))$ therefore eventually equals $C$. This supplies an unspecified finite real number field attaining the plane's full chromatic number. Here finite means finite degree over $\mathbb Q$, not finitely many elements. It gives no field, degree, effective stabilization test, or index bound. In particular, a field that already has a six-chromatic witness need not be the final field if $C=7$.

## 5. Necessary filters and boundary interfaces

A vertex-critical $k$-chromatic graph is connected and has minimum degree at least $k-1$: otherwise a colouring of $G-v$ extends to $v$. Consequently

$$
|E(G)|\ge \left\lceil (k-1)|V(G)|/2\right\rceil.
$$

Two distinct plane points have at most two common unit-distance neighbours. Double counting length-two paths gives

$$
\sum_v\binom{\deg v}{2}\le 2\binom N2=N(N-1).
$$

For $k=6$ this gives $N\ge11$; equality forces a five-regular graph on eleven vertices, impossible by the handshake identity. Thus $N\ge12$. For $k=7$ it gives $N\ge16$. These inherited elementary bounds are weak filters, not assertions of realizability at those orders.

A graph of treewidth at most $w$ is $(w+1)$-colourable. For completeness, remove redundant contained leaf bags of a tree decomposition; a remaining leaf bag has a vertex absent from its parent, whose neighbours all lie in that bag and hence number at most $w$. Deleting that vertex and repeating yields a greedy colouring in reverse order. Therefore a $k$-chromatic witness has treewidth at least $k-1$. A bounded-width proof search below that threshold cannot find it.

The complementary [boundary-palette construction](boundary_palette_compression_20261011.md) retains all realizable equality partitions on the marked neighbour set, together with lifts and exact geometric ports. Our deletion witnesses supply genuine colourings of $G-v$ to such an interface. They cannot be glued to unrelated coordinates merely because their abstract boundary partitions agree. A unit-distance neighbourhood is locally bipartite, but external constraints can force every global colouring of $G-v$ to use the full palette on that neighbourhood; the signature retains precisely that external restriction.

## 6. Reproduction and exact finite scope

[Script](../../04-computation/experiments/plane_algebraic_witness_20261011.py) and [output](plane_algebraic_witness_20261011.out):

~~~
python -B -X utf8 04-computation/experiments/plane_algebraic_witness_20261011.py
python -B -O -X utf8 04-computation/experiments/plane_algebraic_witness_20261011.py
~~~

The universe is all 1,099 labelled simple graphs on one through five vertices at palette threshold three; 672 nonbipartite inputs; independent colour truth tables for every proper induced subgraph of every extracted core; the exact seven-vertex, eleven-edge Moser spindle over $\mathbb Q(\sqrt3,\sqrt{11})$; the spindle with one isolated extra algebraic point; all $3^7=2,187$ one-hot assignments of its 21-variable, 61-clause three-colour CNF; and 17 malformed packet/search controls.

The extracted small cores have orders three and five. The spindle has a 139-node three-colouring refutation, a four-colouring, and seven checked three-colour deletion witnesses. Scaling its coordinates by two while retaining its old graph is rejected, as are coincident vertices, missing proof branches, wrong deletion colourings, an invalid root interval, and boolean/float label aliases.

Normal and optimized execution agree on **8,548 exact checks**. These controls validate the certificate mechanism. They neither construct a six-chromatic graph nor determine whether the plane needs six or seven colours.
