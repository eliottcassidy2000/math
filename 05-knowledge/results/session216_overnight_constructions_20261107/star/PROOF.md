# Matching degree and the two-run star normal form

**PROVED by elementary arguments, with finite exact controls; independently audited.** The matching characterization was proposed by the parent Session216 worker and developed here into a complete support theorem and an explicit arithmetic factorization. No literature-priority claim is made.

## 1. Exact support of the module-run polynomial

Let G be a nonempty tournament with a marked vertex v. For the usual occupied-gap expansion

    H(G[U at v]) = sum_j c_j F_U(j),

F_U(j) counts ordered spanning covers by j nonempty directed paths; c_j counts skeletons with j occurrences of v separated by nonempty paths of outside vertices. The c_j are computed by the core-gap polynomial in the interval packet. The empty-input coefficient is c_0=H(G-v), with H(empty)=1.

Set A=N^+(v), B=N^-(v). Form the bipartite graph with parts A,B and edges ab exactly when a→b in G. Let ν be its maximum matching size.

**Theorem 1.** c_j>0 exactly for 0<=j<=ν+1. In particular the maximum run count is ν+1, and there are no gaps in its positive support.

*Upper bound.* A skeleton with j occurrences of v contains j−1 pairwise vertex-disjoint intervening outside paths. Each starts in A and ends in B, hence contains an arc from A to B. Choosing one such arc from each path gives a matching of size j−1.

*Lower bound.* For any 1<=j<=ν+1, choose a matching a_i→b_i, i=1,...,j−1. Replace v by j transitive clones v_1,...,v_j. The sequence

    v_1, a_1, b_1, v_2, ..., a_(j−1), b_(j−1), v_j

is a directed path. Insert all remaining outside vertices one at a time using Rédei's insertion argument: a vertex can be prepended, appended, or inserted between adjacent vertices whenever the orientation changes from pointing into it to pointing out of it. This never removes any existing separator between the clones. The final Hamilton path therefore has exactly j clone runs, yielding a positive j-run skeleton. The case c_0>0 is Rédei on G-v. ∎

Only positivity, not the numerical value of c_j, follows from matching number.

## 2. Complete structural characterization of two-run contexts

**Corollary 2.** c_j=0 for every j>=3 iff the bipartite A→B graph has matching number at most one. Its nonempty case is a star: every such arc shares a common endpoint. If G is strong and has at least two vertices, this graph is nonempty, so exactly two runs are possible.

Indeed if A→B is empty then B⇒{v}⇒A, and v lies in a singleton strong component. If a bipartite graph with an edge ab has matching number one, every edge meets a or b; having both an edge ab' and a'b would give a two-edge matching, so all edges share one endpoint.

If the star center x belongs to A, put L=B and R=A minus {x}. Then

    G-x = L ⇒ {v} ⇒ R,      v→x.

The only cross arcs between R and L are L→R. The incidences of x toward L and R are otherwise unrestricted. If the star center belongs to B, pass to the converse tournament first; Hamilton-path counts and all module-run coefficients are unchanged. Thus, up to converse, this one-apex normal form describes every context having a positive second coefficient and no higher coefficients. Empty L or R is allowed in the general formulas below; a positive second coefficient requires x to beat at least one vertex of L.

## 3. Exact arithmetic factorization

Take arbitrary tournaments L,R, possibly empty, one vertex v, and one apex x. Keep

    L ⇒ {v} ⇒ R,    L ⇒ R,    v→x,

and choose arbitrary directions between x and L,R. Let S_L be the vertices of L beaten by x, and I_R the vertices of R beating x. Define:

- α_L: the number of ordered spanning covers (P,Q) of L where Q is nonempty and starts in S_L, and P may be empty.
- β_R: the number of ordered spanning covers (P,Q) of R where P is empty or ends in I_R, and Q may be empty.
- τ_L: the number of ordered spanning covers (P,Q) of L where both are nonempty, P ends in S_L, and Q starts in S_L.

All paths are directed; an empty path has one choice. Write h_L=H(L), h_R=H(R), and h_Rx=H(G[R∪{x}]).

**Theorem 3.** The root-v coefficients are

    c_2 = α_L β_R,
    c_1 = 2α_Lβ_R + h_L h_Rx − h_R τ_L,
    c_j = 0  (j>=3).

Consequently for every nonempty tournament U,

    H(G[U at v])
      = α_Lβ_R Q_U(2) + (h_L h_Rx − h_R τ_L) H(U).

This is a constructive identity, not a claim that α,β,h_L,h_R,h_Rx,τ may be selected independently.

*Proof.* In any Hamilton path, remove x. The remaining vertices form either one path or two ordered paths, the first ending at a vertex beating x and the second starting at a vertex beaten by x. Each path traverses its nonempty portions of L,U,R in this order. Since U→x, a second U-run is possible only if both paths contain U; in that case the second path must start in S_L. The first path's R-portion is empty or ends in I_R. Choosing the two L-portions gives α_L choices; choosing the two R-portions gives β_R choices. This proves c_2=α_Lβ_R. There cannot be three U-runs because removing x leaves at most two paths.

For c_1, the single U-run lies either before or after x.

If it lies before x, and the path portion after x uses L, the count is α_Lβ_R. If the portion after x uses no L, all L is before x and the R vertices together with x form an arbitrary Hamilton path, giving h_L h_Rx choices. These are disjoint and exhaustive, including the possibility that x is last.

If the U-run lies after x, the portion after x must start in S_L. When the portion before x uses R, its count is α_L(β_R−h_R): the term h_R removes precisely the case of an empty first R-portion. When that portion uses no R, every R vertex lies after x. The L-portion before x is empty or ends in L minus S_L. Of the α_L choices, exactly τ_L fail this condition, so this case contributes (α_L−τ_L)h_R. These cases also include x first. The total contribution with U after x is α_Lβ_R−τ_L h_R. Adding both orientations proves the formula. Finally F_U(2)=Q_U(2)−2H(U). ∎

## 4. A freely prescribed coefficient, and its coupled offset

If L is transitive on l vertices indexed 0,...,l−1 with increasing edges, then

    α_L = sum_(i in S_L) 2^(l−i−1).

Choose the smallest vertex of the second path Q: everything before it lies in P and every later vertex independently chooses P or Q. Thus α_L ranges over every integer from 0 through 2^l−1. In particular the leading two-run coefficient α_Lβ_R admits a freely chosen binary factor.

However the offset remains coupled. For this transitive L,

    τ_L = sum_(i<j, i,j in S_L) [2^(j−i−1) + indicator(j=i+1)].

To see this, prescribe the first vertex i of Q and last vertex j of P. For i<j, the vertices strictly between them have two assignments; for i>j, a partition is possible only when i=j+1, with one assignment. The two possibilities combine into the displayed sum. In particular τ is even. Taking R empty gives

    H(G[TT_k at v]) = α_L 2^k + 1 − τ_L.

This is a recovered one-apex special family, not a newly proved completeness theorem. Specifying α does not independently specify τ.

## 5. Exact controls and interval connection

`degree_check.py` compares the support theorem with explicit core-gap enumeration for every labeled tournament of orders1..5 and every marked vertex. `check.py` independently uses literal vertex-subset path dynamic programming for all L,R tournaments of orders0..3, every apex cut, and TT_1/TT_2 root modules: 11,250 counts.

`INTERVAL_FACTORS.json` decomposes all eleven cores of the interval packet into this normal form, passing to converse only when needed. The common coefficient36 comes from several actual factorizations αβ:

- 36×1 supplies offsets −7,−3,3,5,13;
- 18×2 supplies −5;
- 6×6 supplies −1,9;
- 12×3 supplies 1,7,11.

This explains why sharing the same two-run coefficient need not mean similar surrounding graphs or a one-edge reconfiguration. It gives an arithmetic interface for further constructions, but no returning operation that expands the eleven-value width has been proved.

## Scope and inheritance

The occupied-gap expansion and path-cover substitution kernel are inherited. The maximum-matching support criterion, star normal form, and factorization are the scoped results developed here; no assertion of literature priority is made. A matching certifies available run counts and discards multiplicity information. The factorization restores endpoint/cover sidecars, rather than treating the matching number as an arithmetic count.

The separate interval packet has its own literal witnesses and proof. This structural packet is not needed to validate those counts and should be independently audited before canon promotion.

## 6. A returning append-word compiler and four persistent choices

The star interface has a useful exact update. Let L have an outcut S_L and state

    (h,α,τ,e),   e=sum_(v in S_L) end_v(L).

Append a new sink z to L (every old L vertex beats z), and put z in S_L precisely when ε=1. Then

    h' = h,
    α' = 2α+εh,
    τ' = τ+ε(α+e),
    e' = εh.                                      (A)

*Proof.* A Hamilton path of the new L ends at z, proving the h and e laws. In the covers counted by α, an old nonempty Q can receive z in either path, yielding2α; the case Q={z} gives εh. For τ, if z lies in P then its required last-vertex condition contributes εα. If z lies in Q and Q already had an old vertex, the contribution is τ. The remaining possibility Q={z} contributes εe. ∎

This is a returning physical update, not merely a recurrence fitted to counts. Distinct inputs sharing (h,α,e) retain equal (h,α,e) after every common binary word, while their τ differences remain unchanged.

Four explicit cut-labeled seeds have common (h,α,e)=(9,245,9) and successive τ values150,152,154,156:

| τ | outgoing rows of L | cut bitmask S_L |
|---|---|---|
|150|0,5,9,19,39,15,63|105|
|152|0,1,3,23,71,31,175,63|77|
|154|0,5,17,7,43,15,63|101|
|156|0,1,11,35,15,87,31|109|

For any finite binary word w, apply(A) to each seed using the same bits, and write their common α as α_w and the first seed's τ as t_w. Closing the star with R empty and any transitive module of size k>=1 gives four actual counts

    α_w 2^k + 9 − t_w − 2j,     j=0,1,2,3.         (B)

More generally, an arbitrary nonempty root module U gives counts α_w Q_U(2)+(9−t_w−2j)H(U), so the spacing is2H(U). For transitive U these are four consecutive odd integers. At the empty word and k=1 they are349,347,345,343. For words of length m, the outputs have at most k+m+9 vertices.

All outputs are strong. `MENU.json` supplies a directed Hamilton cycle for each initial singleton-module core. Appending z retains the old strong core, has an incoming edge from every old L vertex, and has the outgoing edge z→v, so the extended core remains strong regardless of ε. Arbitrary nonempty root substitution preserves strongness as proved in the interval packet.

`menu_check.py` recomputes all four seed states, validates the four core cycles, and checks every binary word of length at most2 and root-module sizes1,2 against literal vertex-subset path counts (56 composite checks). The all-word conclusion follows from the exact update(A), not from those bounded tests. The normal and optimized runs should be compared just as for the other verifiers.

The returning update acts on the cut-labeled source L, with the source/apex/root partition retained. It is not permission to reinsert the unmarked closed output as an arbitrary new source. The four-choice width is preserved, not enlarged. These arbitrary-word families do not prove that their centers fill an interval or that any sufficiently large target appears. They are a composable positive source interface beyond a homogeneous-block ray.
