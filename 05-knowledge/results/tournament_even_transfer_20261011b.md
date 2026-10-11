# Tournament counts, even flows, and marked colouring reconstruction

**Status:** PROVED (the elementary finite-group identities and reconstruction below); CITED / INHERITED (the tournament equinumerosity and current equivariant-map results); FINITE-EXACT (the declared controls). No six- or seven-chromatic unit-distance witness is supplied. No equal-count theorem is presented as a colouring algorithm.

## 1. Inheritance and the two meanings of even

The closest current result is [THM-4531, twisted Mallows–Sloane principle and no natural bijection](../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md), with the detailed [current natural-bijection note](procgen_tbij_20261001_natural_bijections.md). The older [HYP-3799, two even-graph notions](../hypotheses/HYP-3799-two-even-graph-notions-and-the-tournament-equinumerosity.md) correctly separates the notions, but its historical open-status discussion does not supersede those later scoped results.

For a labelled graph $X$ and $g\in\operatorname{Aut}(X)$, define

\[
\epsilon_X(g)=(-1)^{\#\{ij\in E(X):i<j,\ g(i)>g(j)\}}.
\]

This is a character, independent of the chosen reference orientation. Call $X$ **automorphism-even** when it is trivial. Call $X$ **Eulerian** here when every vertex has even degree; disconnected graphs and isolated vertices are allowed. The two predicates differ already on three vertices: a triangle is Eulerian but not automorphism-even; a three-vertex path is automorphism-even but not Eulerian. Empty graphs satisfy both, while one edge plus an isolate satisfies neither.

Royle, Praeger, Glasby, Freedman and Devillers prove equality between the numbers of **unlabelled tournaments** and **unlabelled automorphism-even graphs** of a given order. This is a statement about isomorphism types, not labelled cardinalities or vertex-degree parity. Their proof uses fixed points under relabelling. [Primary paper, Theorem 1.1](https://arxiv.org/html/2204.01947v2).

Current THM-4531 gives the exact signed fixed-point identity, more generally for an $S_n$-invariant edge-reversal subspace $U$:

\[
\#\operatorname{Fix}_{\mathcal T/U}(g)
=\sum_{\substack{X\in U^\perp\\gX=X}}\epsilon_X(g).
\tag{1}
\]

At $U=0$, averaging (1) gives the cited theorem. At $U$ equal to the cut space, $U^\perp$ is the Eulerian cycle space: tournament switching classes are equinumerous with graphs that satisfy **both** predicates. The bridge is a signed character calculation, not a map sending each tournament to a particular graph. On an automorphism orbit, the stabilizer-character average is one or zero according as the character is trivial or nontrivial.

Finite sets of equal cardinality always admit an arbitrary bijection. The stronger word “natural” in the current repository results means a relabelling-equivariant map on labelled objects inducing a bijection on isomorphism types. For $n\ge2$, there is no equivariant map from all labelled automorphism-even graphs to tournaments: the empty graph is fixed by all of $S_n$, while no tournament is. The reverse-direction question has its separate, explicitly bounded proved ranges in THM-4531; this note makes no unrestricted reverse-direction assertion.

The exact small table illustrates the types:

| $n$ | Labelled tournaments | Labelled automorphism-even | Labelled Eulerian | Tournament types | Automorphism-even types | Eulerian types |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 1 | 1 | 1 | 1 | 1 |
| 3 | 8 | 4 | 2 | 2 | 2 | 2 |
| 4 | 64 | 20 | 8 | 4 | 4 | 3 |

The relevant corrected near miss is [THM-4069, even-graph basis dependence](../../01-canon/theorems/THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope.md), corrected through MISTAKE-495. A spanning-tree fundamental-cycle basis gives a labelled cube, but its simple graph image after quotienting by graph isomorphism depends on the tree. For $K_4$, a star basis gives $P_3$ and a path basis gives $K_3$ on the same three Eulerian graph types. The script independently reproduces this. The least-used sidecar worth retaining is the **marked spanning forest and actual edge coordinates**, rather than just the quotient type or cycle-space dimension.

## 2. The exact colouring transform for every modulus

Let $G=(V,E)$ be a finite simple labelled graph with $n$ vertices, $m$ edges and $c$ connected components. Orient each edge from its smaller label to its larger label. Let $q\ge2$ be an integer, not necessarily prime. Write

\[
Z_q(G)=\ker\bigl(\partial:(\mathbb Z/q)^E\longrightarrow(\mathbb Z/q)^V\bigr),
\]

where an edge value contributes negatively at its tail and positively at its head. These are Kirchhoff flows. For $f\in Z_q(G)$ let $s(f)$ be its number of nonzero edge values. The exact chromatic count is

\[
\boxed{\displaystyle
P_G(q)=q^{n-m}\sum_{f\in Z_q(G)}(-1)^{s(f)}(q-1)^{m-s(f)}.}
\tag{2}
\]

**Proof.** For a primitive $q$th root of unity $\zeta$ and $d\in\mathbb Z/q$,

\[
\mathbf1_{d\ne0}=\frac{q-1}{q}-\frac1q\sum_{a=1}^{q-1}\zeta^{ad}.
\]

Multiply this identity over all edges with $d=c_v-c_u$, then sum over all vertex colourings $c$. Orthogonality at each vertex kills exactly the assignments of edge frequencies with nonzero divergence. Every surviving flow contributes $q^n$, yielding (2). This is a finite identity; no convergence or positivity assumption is hidden in the expansion. In particular, it works for $q=6$.

For $q=2$, flow values are edge indicators and conservation means even degree. Thus (2) uses **Eulerian** subgraphs, not automorphism-even graphs. On a triangle the two binary flows contribute $1-1=0$. The presence of nonzero flows does not prove a proper colouring: the signed sum can cancel.

Choose a spanning forest. Its chord count is

\[
\beta=m-n+c.
\]

Assign arbitrary values to those $\beta$ chords and remove forest leaves successively, choosing the unique parent-edge value that cancels each leaf divergence. Only multiplication by $1$ or $-1$ is used, so composite moduli cause no division problem. Component conservation finishes at the roots. Conversely every flow determines these chord values. Therefore $|Z_q(G)|=q^\beta$, with a concrete bijection implemented by `complete_flow`.

This identifies precise triangular dimensions for complete graphs:

\[
|E(K_n)|=\binom n2,\qquad \beta(K_n)=\binom{n-1}2.
\]

They count different coordinates: tournament orientations have $2^{\binom n2}$ labelled states, while Eulerian edge subsets have $2^{\binom{n-1}2}$. Neither formula counts labelled automorphism-even graphs.

## 3. A positive certificate is a marked cut, not a flow

For an edge vector $d\in(\mathbb Z/q)^E$, the following are equivalent:

1. $d_{uv}=c_v-c_u$ for a proper vertex colouring $c$;
2. every $d_e$ is nonzero and every oriented cycle sum of $d$ vanishes.

Given one root colour per connected component, including isolated vertices, integrate $d$ along the marked forest. The remaining chord equations are exactly the cycle conditions. This recovers a **unique** proper colouring with those root colours. Conversely differences of any proper colouring satisfy these conditions. The exported `colour_from_cut` checks the complete supplied edge vector after integration, rather than trusting a cycle summary.

The two constraints are different: a flow satisfies a vertex conservation law; a cut satisfies cycle integrability. For the oriented triangle $(01,02,12)$, the $\mathbb Z/6$ flow $(1,5,1)$ is nowhere zero and has divergence zero, but its oriented cycle sum is $1+1-5=3\pmod6$. It is not a cut, and the colouring decoder rejects it. A flow and a cut can coincide for other graphs or moduli; no universal disjointness claim is intended.

Relabelling must reverse the sign of a nonbinary edge value whenever the reference orientation reverses. `relabel_flow` retains this marking, and the checker independently verifies conservation after every permutation of the triangle.

This is the constructive output of the transfer: an exact edge certificate plus component phases gives the actual colours. To certify a plane unit-distance graph, the original point coordinates and all relevant distance constraints must also remain attached. A signed count or an abstract coloured graph does not supply that geometric realization.

## 4. Six colours: the CRT preserves an overlap coordinate

Coordinatewise Chinese remaindering gives

\[
Z_6(G)\cong Z_2(G)\times Z_3(G),\qquad f_e=3a_e+4b_e\pmod6.
\]

Its zero test is $f_e=0$ iff $a_e=b_e=0$. Therefore

\[
\boxed{\displaystyle
P_G(6)=6^n(5/6)^m
\sum_{a\in Z_2(G)}\sum_{b\in Z_3(G)}
(-1/5)^{|\operatorname{supp}(a)\cup\operatorname{supp}(b)|}.}
\tag{3}
\]

The edge-by-edge **support union** is essential. Separate support sizes or independent support enumerators do not identify each paired term. Take two triangles $012$ and $034$ sharing vertex $0$, in edge order $(01,02,03,04,12,34)$. Let

\[
a=(1,1,0,0,1,0),\quad
b=(1,2,0,0,1,0),\quad
b'=(0,0,1,2,0,1).
\]

All three are legitimate flows in their respective moduli. The pairs $(a,b)$ and $(a,b')$ both have marginal support sizes $(3,3)$, but support-union sizes $3$ and $6$. Their weights in (3) are $-1/125$ and $1/15625$, of opposite signs. The cheapest repair is to retain marked edges, or the full joint overlap enumerator required by the sum.

Likewise a proper six-colouring corresponds to binary and ternary vertex labels whose **joint** edge difference is nonzero. Neither component must separately be a proper colouring. A triangle has no proper binary colouring but has $6\cdot5\cdot4=120$ proper six-colourings.

## 5. What this can and cannot buy for six versus seven

The transform is exact but does not promise a smaller exhaustive search. A connected vertex-critical $k$-chromatic graph has minimum degree at least $k-1$: otherwise a colouring after deleting a low-degree vertex extends. Hence

\[
\beta=m-n+1\ge\frac{k-3}{2}n+1.
\tag{4}
\]

For $k=6$ and $7$, this is at least $3n/2+1$ and $2n+1$, respectively. Thus direct enumeration of all $q^\beta$ flows can be larger than the $q^n$ vertex-colour enumeration on the critical witnesses of interest. The useful regime is a graph or fragment of low cycle rank, or a factored computation that retains exact boundary data; replacing it by the abstract dimension alone does not help.

Even imposing both even-graph predicates gives no abstract chromatic bound. For any $k\ge2$, take the complete $k$-partite graph with two vertices in every part. Its degrees are $2(k-1)$, so it is Eulerian. Its automorphism group is generated by within-part swaps and whole-part swaps. A within-part swap reverses no edge between those two vertices; a swap of two contiguous parts in the reference order reverses four interpart edges. Both signs are positive, proving automorphism-evenness. Colouring the parts and taking one vertex from every part prove chromatic number exactly $k$. This is an **abstract graph** family, not a claim that these graphs embed at unit distance in the plane.

The live concept board is consequently: **automorphism character / Eulerian cycle space / marked cuts / joint CRT support / geometric realization**. Their valid maps are:

| Source → target | Map and preserved predicate | Lost information | Required sidecar / decisive hostile |
|---|---|---|---|
| Tournament relabelling action → signed graph sum | Walsh-character trace and Burnside; unlabelled cardinality | Individual object correspondence | Stabilizers; empty graph has no tournament-valued equivariant image |
| Vertex colours → flow partition function | Finite Fourier expansion; exact number of proper colourings | A particular colouring, positivity of individual terms | Signed weights; binary triangle cancels |
| Nonzero cut → proper colours | Forest integration; all edge inequalities | Component additive phases | One root colour per component and all chord equations |
| Mod-6 flows → mod-2 × mod-3 flows | CRT; conservation | Overlap if only marginals retained | Actual edge supports; the two-triangle hostile |
| Marked cycle basis → graph-type quotient | Projection to isomorphism classes | Generator adjacency | Original basis; $K_4$ star/path yields $P_3/K_3$ |

The transfer supplies exact alternative coordinates and independently checkable finite colour certificates. It supplies neither a universal plane six-colouring nor a seven-chromatic plane witness.

## 6. APIs, exact universe, and reproduction

[Companion experiment](../../04-computation/experiments/tournament_even_transfer_20261011b.py). Public colouring/flow inputs use exact integer labels and moduli, canonical sorted edge tuples and canonical residue tuples. Booleans and floats are rejected before equality can alias valid data. `colour_from_cut` also requires the complete root-phase vector and rejects zero edge differences or failed chord equations. Resource use is exponential when the caller requests the full flow or colouring enumeration; no large-graph efficiency is claimed.

The checks are independent finite paths, not a proof by census:

- all labelled graphs of orders $1,2,3,4$ and all $q=2,\ldots,7$: 75 graphs, 450 graph/palette cases, literal colour counts versus (2), with a separate CRT evaluation of (3);
- all relabellings in the order-$\le4$ count tables, checking signed fixed points before comparing unlabelled totals;
- 15,496 proper colourings of a triangle and a disconnected five-vertex graph, reconstructed exactly from edge differences and component phases;
- the marked-basis, overlap, nowhere-zero-flow, and both-even/chromatic hostiles above, plus malformed certificate types.

Reproduce from the repository root:

```powershell
python -B -X utf8 04-computation/experiments/tournament_even_transfer_20261011b.py
python -B -O -X utf8 04-computation/experiments/tournament_even_transfer_20261011b.py
```

The saved [output](tournament_even_transfer_20261011b.out) records **19,067 exact checks**. Normal and optimized runs agree byte for byte; the LF-normalized SHA256 is `5ea9e6462d4afb81beb6c835c8dfcbfec7850807f57540a7c98362211f63c8de`. The all-length arguments are the proofs above and the cited current theorem, not an extrapolation of this finite universe.
