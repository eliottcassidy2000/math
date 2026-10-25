# Exact counts of labelled path-locked score fibres

Status: PROVED all-size identities, classical graphical-zonotope/Tutte theory, and FINITE-EXACT independent checks through n=14. No novelty claim; no claim of equivalence with A000568, tournament isomorphism classes, score multisets, or Hamilton-count classes.

## 1. Scores are exactly the lattice points

Let G be a finite simple graph on labelled vertices, and let

    Z_G = sum_{uv in E(G)} conv(e_u,e_v).

An orientation contributes e_u for an edge directed out of u, so its labelled outdegree vector is a lattice point of Z_G. Conversely, take an integral point d in Z_G and represent it by one fractional endpoint choice per edge. Choose reference orientations and write d=b+B x, 0≤x≤1, where B has columns e_u−e_v and b is integral. Consider the edges with strictly fractional x_e. If their undirected graph has a leaf, that vertex's equation Bx=d−b consists of integral contributions plus exactly one nonintegral contribution, impossible. Thus if any fractional edge remains, their subgraph has a cycle. A signed cycle vector lies in ker B; perturb x in that direction until some fractional coordinate first reaches 0 or 1. No integral coordinate changes, all coordinates stay within [0,1], and at least one fractional coordinate disappears. Repeat finitely many times. The result is an orientation with score d.

This elementary rounding proof also handles parallel edges (a pair can form a 2-cycle). It is not a claim that every lattice point of an arbitrary integral polytope is a projected cube vertex.

## 2. Exact Ehrhart / forest / Tutte theorem

For every nonnegative integer t,

    E_G(t) := #(t Z_G intersect Z^V)
            = sum_{F subset E(G), F acyclic} t^{|F|}.

All vertices are present in a spanning forest, including isolated vertices. Forests here are arbitrary acyclic edge subsets, not necessarily induced subgraphs or spanning trees. Thus the number of nonempty labelled score fibres is E_G(1)=T_G(2,1). More generally, with r=|V|−c(G),

    E_G(t) = t^r T_G(1+1/t,1)

as a polynomial identity (the displayed rational substitution is interpreted after simplification at t=0).

### Self-contained deletion–contraction proof

Translate Z_G integrally to the signed graphic zonotope sum_e [0,e_u−e_v]; integer translation preserves every integer dilation count. Choose a nonloop edge e=uv and put v=e_u−e_v. This is a primitive lattice vector. Project along v onto the quotient lattice Z^V/Zv. The projection of the deletion zonotope is the contraction zonotope: the quotient identifies u and v, and any parallel edges remain distinct generators; loops become zero generators.

For any quotient lattice point y in tZ_{G/e}, its fibre in tZ_{G−e} is a nonempty compact real interval [a,b] in a coordinate whose lattice is Z. Adding [0,tv] changes this interval to [a,b+t]. Its lattice-point count increases by exactly t:

    (floor(b+t)−ceil(a)+1) − (floor(b)−ceil(a)+1) = t.

The formula also holds if the original interval has no lattice points: a nonempty real interval has floor(b)−ceil(a)+1≥0. Summing fibres proves

    E_G(t)=E_{G−e}(t)+t E_{G/e}(t).

A loop has zero generator and can be deleted without changing E. With no edges E=1. The forest polynomial has the same initial conditions and recurrence, by separating forests not containing e from those containing e and contracting e. This proves the identity, including disconnected graphs and intermediate multigraphs. The Tutte specializations follow directly from its rank-subset definition: at y=1 only independent subsets survive.

### Classical provenance

Richard P. Stanley, “A zonotope associated with graphical degree sequences,” Applied Geometry and Discrete Combinatorics, DIMACS Series in Discrete Mathematics 4 (1991), 555–570, Theorem 2.2, gives the general lattice-zonotope Ehrhart formula in terms of independent generating sets and their maximal-minor gcds. For graph incidence vectors, independent sets are forests and each such gcd is 1. This is a classical specialization, not a new theorem.

Primary paper: https://math.mit.edu/~rstan/pubs/pubfiles/83.pdf
Author's bibliographic entry: https://math.mit.edu/~rstan/pubs/ (entry 83).
Author's later Ehrhart lecture slides, containing the general formula: https://math.mit.edu/~rstan/transparencies/ehrhart1.pdf

The proofs above and below are supplied here independently; they do not require treating the citation as an unproved black box.

## 3. Hamilton-path locks

Fix 0→1→...→n−1 in a labelled tournament. The free graph is G_n=K_n−P_n. Locked outdegrees contribute one fixed vector, so full tournament-score fibres correspond bijectively to residual-score fibres on G_n. Consequently the number of nonempty fibres is exactly the number of spanning forests of G_n, irrespective of the chosen directions of the locked path edges. The same statement holds for any fixed locked edge set, though the formulas below exploit a path specifically.

For n=5,6,7 the full polynomials are:

    E_5(t)=1+6t+15t²+19t³+11t⁴, total 52.
    E_6(t)=1+10t+45t²+116t³+176t⁴+130t⁵, total 478.
    E_7(t)=1+15t+105t²+445t³+1224t⁴+2109t⁵+1807t⁶, total 5706.

The score-fibre count is not the number of orientations (2^((n−1)(n−2)/2)), nor does it count orientations within a prescribed fibre.

## 4. Exact path inclusion–exclusion and compressed recurrence

Write F_n(t)=sum_{F forest of G_n} t^{|F|}. For a chosen subset S of forbidden path edges, its components are consecutive blocks with ordered sizes s=(s_1,...,s_k), an ordered composition of n; |S|=n−k. In a forest of K_n containing S, contracting S leaves a forest on k distinct blocks. Between blocks i,j there are s_i s_j possible original edges. No additional edge inside a block is permitted, since it would create a cycle; neither may two parallel contracted edges both be chosen. Let Q_s(t) be the weighted forest polynomial on these blocks with edge weight s_i s_j. Inclusion–exclusion gives

    F_n(t)=sum_{s composition of n} (-t)^{n−length(s)} Q_s(t).

For a nonempty set B of blocks define

    W_s(B)=1                                      if |B|=1,
    W_s(B)=(sum_{i in B} s_i)^{|B|−2} prod_{i in B}s_i  otherwise.

Weighted Prüfer coding proves W_s(B) is the sum of tree weights on B: each tree contributes prod_i s_i^{degree(i)}, and its Prüfer word supplies the factors beyond the first s_i at each vertex. Therefore

    Q_s(t)=sum_{pi set partition of [k]} prod_{B in pi} W_s(B) t^{|B|−1}.

In particular the exponent in each full inclusion–exclusion summand is n−|pi|, as required. This is a proved exact formula, not a rooted-forest determinant substitution.

### Compress compositions to integer partitions

Q_s depends only on the multiset of sizes. If a_i is the number of size-i blocks and k=sum_i a_i, then the number of ordered compositions with this multiset is k!/prod_i a_i!. Thus

    F_n(t)=sum_{a: sum_i i a_i=n} (-t)^{n−k} (k!/prod_i a_i!) Q_a(t).

Compute Q_a without enumerating set partitions. Fix one distinguished block of size r with a_r>0. Choose b_i other blocks of each size to join its component, 0≤b_i≤a_i−delta_{ir}. Write h=sum b_i, M=r+sum i b_i. Put W=1 if h=0 and W=M^{h−1} r prod_i i^{b_i} otherwise. Then

    Q_a(t)=sum_b [prod_i binom(a_i−delta_{ir},b_i)] W t^h Q_{a−e_r−b}(t),
    Q_0(t)=1.

The binomial coefficient counts choices among distinct equal-sized blocks. Memoization needs only multisets of total mass at most n, bounded by sum_{m=0}^n p(m) states. This is a genuine state-space compression, not a claim of polynomial time: each state still sums over its submultisets, and coefficients use arbitrary-precision integers. This recurrence and the subset method below are mathematically distinct exact implementations.

## 5. Independent subset/tree recurrence and path-run acceleration

For any graph G and U subset V, let A_U(t) be the forest polynomial of G[U]. With r=min U,

    A_empty(t)=1,
    A_U(t)=sum_{B subset U, r in B} tau(G[B]) t^{|B|−1} A_{U minus B}(t),

where tau(singleton)=1 and disconnected induced graphs have tau=0. This selects the tree component containing r uniquely. Precomputing all induced tree counts using the matrix-tree theorem gives a subset recurrence with O(3^n) component terms, up to polynomial overhead for polynomial arithmetic. This improves the exponential parameter from the Θ(n²) free-edge count to n, without promising practical scalability at arbitrary n.

For G_n, let s=|B| and H=P_n[B]. Then L_{G_n[B]}=sI−J−L_H. The common constant direction is annihilated by both Laplacians; on its orthogonal complement the first matrix is sI−L_H. Hence, including disconnected cases,

    tau(G_n[B])=det(sI−L_H)/s².

For s=1 the right side is 1 as well. H is a disjoint union of path runs, so its determinant factors over these runs. For a run of length l=1 its factor is s. For l≥2, use the continuant D_0=1, D_1=s−1, D_j=d_j D_{j−1}−D_{j−2}, with d_j=s−2 at internal vertices and d_l=s−1 at the final endpoint. This computes each subset tree count in O(n) integer arithmetic operations; division by s² is verified exactly in the checker. Full Bareiss-cofactor tree counts supply a third check through n=12.

## 6. Certified counts and controls

The executable standard-library checker checks:

- Direct enumeration of all orientations and separately all edge subsets with union-find acyclicity, n=1,...,7.
- Weighted composition inclusion–exclusion (compressed to integer partitions), n=1,...,14.
- Independent subset/tree recurrence with path-run continuants, n=1,...,14.
- Independent Bareiss matrix-tree cofactors inside the subset recurrence, n=1,...,12.
- Exact full-coefficient agreement, not only total agreement.

Totals for n=1,...,14 are:

    1, 1, 2, 8, 52, 478, 5706, 83824, 1461944, 29510268,
    676549450, 17361810016, 492999348348, 15345359136232.

The full coefficient rows are in certificate.json and run.out. The observed combined checker run took under one second on this execution host; certificate.json records its own perf_counter duration. This is a reproducibility observation, not a general performance guarantee.

### Hostile determinant control

Do not replace the unrooted forest polynomial by det(I+tL_G). The latter is the rooted/component-size-weighted forest enumerator:

    det(I+tL_G)=sum_F t^{|F|} prod_{C component of F}|C|.

For the triangle, the desired polynomial is 1+3t+3t², total 7; the determinant is 1+6t+9t², total 16. The checker explicitly verifies 7 versus 16. Isolated vertices remain components of size one in both conventions. A single matrix-tree cofactor does correctly count spanning trees; the error is using the all-roots determinant for unrooted forests.

## Reproduction

    python check.py

No third-party Python dependencies, floating point mathematics, random sampling, remote writes, or modifications to the cycle-move packet are involved. This packet supplies the number of nonempty labelled score fibres, complementing (but not replacing) the cycle-reversal connectivity theorem.
