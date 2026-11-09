# Independent audit of the OCF attachment kernel and connectivity threshold

Date: 2026-11-09. Verdict: **ACCEPT, with a singleton-domain qualification**.

Sources inspected: `session230-exclusion-bridges/PROOF.md`, its scripts and result files; `session230-generation-audit/REVIEW.md`; and the actual canon files for THM-4097, THM-4114 and THM-4118. This audit neither changes that repository nor imports producer code. The OCF identity is an established premise, not newly proved here. No all-odd realization or all-order first-creator genealogy is established.

## 1. All-order kernel identity

Let H count directed Hamiltonian paths, with H(empty)=1. For distinct old vertices a,b, let K(a,b) sum H(T-V(P)) over all directed a-to-b paths having an even number of vertices. An odd-cycle family in the extended tournament either omits the new vertex v or contains a unique cycle through it. In the second case, cut that cycle at v. The resulting path P has even order, starts at a vertex beaten by v, and ends at a vertex beating v. Its residual odd-cycle family has weight sum H(T-V(P)), and the removed cycle supplies exactly a factor two. This is a bijective partition, without path-deletion assumptions or multiplicity ambiguities: a directed cycle is considered modulo cyclic rotation, and v fixes its rotation.

Consequently, for the set A beaten by v,

    H(T+v_A) = H(T) + 2 sum_{a in A,b outside A} K(a,b).

The proof is valid at all finite orders. In particular K is nonnegative and integral. The actual directed arc between a and b is a two-vertex path, so

    K(a,b)+K(b,a) >= H(T-{a,b}) >= 1.

The last inequality is Redei's Hamiltonian-path theorem (also follows from OCF's empty family). Expanding y_a(1-y_b) gives constant mixed attachment derivative -2(K(a,b)+K(b,a)); the quantitative strict-submodularity bound follows. Complementing all attachment coordinates to x_a=[a->v] preserves pair curvature and changes the singleton coefficient to twice the incoming K-column sum. Thus the stated split-path coefficient bridge is correct, provided L and M retain the stated convention.

The uniqueness invoked here is uniqueness of a multilinear Boolean polynomial, **not** uniqueness of a directed-capacity matrix realizing that polynomial. No reconstruction of the tournament from the attachment response has been established.

## 2. Neutrality and first birth

A positive directed arc crossing A to its complement gives a positive two-vertex contribution. If no such arc exists, no directed path can cross the partition in that direction. Therefore the neutral cuts are exactly the terminal SCC unions. This proof is independent of the historically invalid claim that deleting a vertex from a Hamiltonian path always leaves one.

For a strong tournament with at least two vertices, every nontrivial cut is positive, and the minimum strictly increasing extension equals H(T)+2 mincut(K). In tournament language this means a strong base of order at least three. **Qualification:** the singleton is often called strong, but has no nontrivial cuts, so formulas (5) and (6) should explicitly exclude it or declare an empty-minimum convention. The substantive statements are unaffected.

The reduction to ordered source/sink mincuts is sound: every nontrivial cut separates an ordered pair, and every feasible source/sink cut is nontrivial. Computing the capacities themselves may be expensive; no polynomial construction of K is proved.

Padding any tournament with source or sink vertices preserves H. Hence every value first attained before order n belongs to S(n). If H(T)+2 mincut(K)>max S(n), every nonconstant attachment is a first creator at n+1. This is a sufficient certificate only. Values in interior holes require membership testing against the entire S(n). The certificate does not prove existence of a fertile attachment, genealogy for every child, or sufficiency of one representative per first-born value.

## 3. The 31-word exclusion is a complete all-order subclass proof

For a transitive base labeled in increasing orientation, every residual tournament has H=1. An even-order i-to-j path chooses an even-size subset of the j-i-1 interior vertices. Thus K(i,i+1)=1, K(i,j)=2^(j-i-2) for j>=i+2, and reverse entries vanish.

In a binary cut word, initial zeros and terminal ones contribute no forward cut arc; removing them does not change distances between retained positions. A nonzero cut reduces to a word starting in 1 and ending in 0. If its length m is at least seven, the endpoint contribution alone is 2^(m-3)>=16. A cut of value 3 or 10 must therefore have 2<=m<=6. There are exactly 1+2+4+8+16=31 residual words. Independently enumerating them gives:

- m=2: 1
- m=3: 2
- m=4: 4,5
- m=5: 7,8,9,11
- m=6: 12,13,15,16,18,19,22,23,25

Neither target occurs. Since the extension H is 1+2 cut, neither 7 nor 21 occurs at any order in the apex-transitive subclass. A triangle transversal of size at most one is equivalent to this subclass: deleting its vertex leaves no cyclic triangle, hence a transitive tournament. Empty or transitive degeneracies have H=1. This proof never invokes the global 7/21 exclusion, and is therefore genuinely noncircular **within its stated subclass**. It is not a proof of the general all-order omission.

## 4. Support geometry and inherited scope

For a disjoint odd-cycle family F, use its union U(F) as the squarefree monomial support and coefficient 2^|F|. Summing these monomials evaluates to H(T[S]) on the indicator of S by OCF. Inclusion-exclusion on a marked set R keeps exactly families satisfying R subset U(F) subset B union R. This proves every mixed support difference nonnegative, including the explicit family-cover formula. It differs from attachment orientation variation, whose pair derivatives are negative. A family covering R can involve multiple disjoint cycles; replacing it by one cycle covering R would be incorrect. The example of two ordinally summed triangles with one marked vertex in each gives mixed derivative four.

Minimal nonconstant supports are exactly cyclic triangles: every nontransitive tournament has one, and a triangle has no smaller cycle support. Thus hitting those supports is exactly feedback vertex deletion. The full polynomial still contains multiplicities and longer-cycle information absent from the unweighted clutter.

**Attribution:** THM-4097 already develops the exact split-path cut/field calculus. THM-4114 already proves nonnegative induced-support Mobius atoms, all mixed support differences, the opposite Boolean curvatures, and a nonnegative cut representation with auxiliary terminals. THM-4118 already warns that gcd/lattice information does not imply a solid response spectrum. Accordingly, these results must not be advertised as discoveries of the current session. The useful present synthesis is the explicit original-vertex even-path capacity formula with a positive combinatorial witness, its quantitative bound and mincut birth interpretation, and the finite-word subclass certificate. This is an assessment of the inspected canon, **not a literature-priority claim**.

## 5. Sharp semidegree/connectivity threshold

For k>=1 and n>=2, define lambda(T) as the minimum directed cut size and delta0 as minimum in/outdegree. If n<=4k+1, then lambda>=k iff delta0>=k.

Singleton and complementary singleton cuts give necessity. For sufficiency, suppose a cut X,Y has r<k arcs Y->X and delta0>=k. Put a=|X|. Summing indegrees inside X gives

    ak <= a(a-1)/2 + r.

For integers 1<=a<=2k, the concave expression a(2k+1-a)/2 has minimum k at the endpoints. Thus a>=2k+1. Summing outdegrees in Y gives |Y|>=2k+1 as well. Hence n>=4k+2, contradiction. This verifies the argument without relying on any finite census.

At order 4k+2 take two cyclic regular tournaments of order 2k+1, orient all cross arcs A->B, and reverse one cross arc. For k>=2, both blocks are strong, both crossing directions exist, and the result is strong. Its semidegree is exactly k, but the reverse cut has one arc, so lambda=1<k. For k=1, the ordinal sum of two triangles without reversal is the sharp nonstrong obstruction at six vertices; no strong obstruction can exist because lambda>=1 is exactly strong connectivity. These constructions are sound. No literature novelty is asserted.

The related f_k(n)>max S(n-1) implication does not need the order threshold: by definition any non-first-creator has lambda<k, and reversing the fewer-than-k backward arcs of a minimizing cut gives an ordinal sum. The order threshold merely allows replacement of cut constraints by semidegree constraints.

## 6. Independent executable checks and limits

`check.py` is standalone standard-library Python. It uses explicit recursive enumeration of simple paths and direct permutation counts for extended tournaments, rather than importing the producer's subset-DP implementation. It independently enumerates odd directed cycles with a canonical starting vertex and all disjoint cycle families, checking the full OCF identity and exact support-family derivative equality. It checks all labeled bases of orders 1..5, every attachment, every coordinate-pair curvature face, neutrality, the strict lower bound, all marked support derivatives, padding, and every applicable mincut birth certificate. Separate exhaustive cut tests cover every labeled tournament of orders 2..6. Sharp examples cover k=1..4, while the proof covers every k. The 31 residual words are fully enumerated.

Run `python3 check.py > RESULTS.json` in this directory. Counts and outputs are recorded there. All checks use explicit `require`/raise guards and remain active under Python's `-O` flag. Normal and optimized runs were verified byte-identical to the frozen `RESULTS.json`. The large order-eleven census is outside this audit. No bounded enumeration here is offered as a substitute for the all-order arguments, and all-odd realization remains **OPEN**.

## Recovered prior connectivity credit

Concurrent S17 audit identified THM-4811(B) as the prior source of the same degree/cut threshold. The independent proof check remains valid, but this argument is inherited rather than new.
