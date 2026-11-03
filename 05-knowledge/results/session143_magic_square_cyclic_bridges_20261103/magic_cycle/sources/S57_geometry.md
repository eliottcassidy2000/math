# Square-root varieties for bare square-sum bipartite graphs

Scope: exact characteristic-zero geometry, a rational-point reformulation, and verification of a cited, previously known nondegenerate K4,4. No magic square of squares, new integer prize solution, or resolution classification is established. The bare-variety arguments in Sections 1–3 also apply geometrically over any field of characteristic different from 2. The arithmetic equivalence and source comparisons are stated over Q.

## 1. Complete intersection and dualizing sheaf

Fix r,s >= 2, N=rs, n=r+s-1, and m=(r-1)(s-1)=N-n. In P^(N-1) let X_(r,s) be defined by

f_ij = q_ij^2 + q_11^2 - q_i1^2 - q_1j^2 = 0,  i,j >= 2.

The boundary variables q_i1 and q_1j (counting q_11 once) form a polynomial ring B in n variables. Order the m interior variables arbitrarily. Each successive equation is monic quadratic in its own interior variable, and earlier equations do not involve that variable. A monic polynomial is a non-zero-divisor over any coefficient ring: the highest nonzero coefficient of a putative annihilator remains the highest coefficient of the product. Thus the f_ij form a regular sequence. The quotient is free over B with basis the squarefree monomials in the m interior variables. It has dimension n, and Proj has dimension n-1=r+s-2 and codimension m. In particular this is a scheme-theoretic complete intersection, not merely an expected dimension count. Its degree is 2^m, and its Hilbert series is (1-t^2)^m/(1-t)^N.

The complete-intersection dualizing-sheaf formula gives

omega_X = O_X(2m-N) = O_X(rs-2r-2s+2).

For completeness this formula follows by taking the Koszul resolution of O_X in P^(N-1), applying Hom(-,omega_P), and taking its unique nonzero sheaf Ext in degree m. Its final twist is omega_P(2m)|_X, with omega_P=O_P(-N). It does not require smoothness.

Cases: X_33 has dimension 4, degree 16, and omega=O(-1); X_34 has dimension 5, degree 64, and omega=O(0); X_44 has dimension 6, degree 512, and omega=O(2). These are dualizing-sheaf statements about the singular models, not assertions that their resolutions are smooth Fano, Calabi–Yau, or general type.

## 2. Geometric integrality

Let L be the n-dimensional vector space of additive matrices z_ij=a_i+b_j; fix b_1=0 to remove the one-dimensional gauge redundancy. Its coordinate ring R is polynomial in a_1,...,a_r,b_2,...,b_s. The coordinate linear forms ell_ij=a_i+b_j are nonzero, pairwise nonassociate irreducibles in R, including ell_i1=a_i. The affine cone of X has the presentation

A = R[q_ij for all i,j]/(q_ij^2-ell_ij for all i,j).

It is free of rank 2^N over R. If a product of a nonempty subset of the ell_ij were a square in Frac(R), its valuation at any selected irreducible ell_ij would be odd, a contradiction. The N square classes are therefore independent. The standard elementary multiquadratic lemma implies that adjoining their square roots has degree 2^N. One proof of the lemma is induction: in a multiquadratic extension the sign-change characters have one-dimensional eigenspaces spanned by the products of the previously adjoined roots; a square root of a base-field element is a character eigenvector, so its square class is in the span of those previously adjoined classes.

Consequently A tensor_R Frac(R) is a field. Freeness embeds A in this field, so A is a domain. The argument works after extension to an algebraic closure, proving geometric integrality.

The finite squaring map X -> P(L) has generic degree 2^(N-1): the affine degree 2^N loses the simultaneous global root sign under projectivization. The earlier degree 2^m is instead the degree of X in its q-coordinate projective embedding; these are different degrees and must not be conflated.

## 3. Singular-locus audit

Choose a full-rank m-by-N matrix C defining the additive linear space L as ker C. The Jacobian of the defining quadrics is 2 C diag(q). At a point, its rank drops below m exactly when there is a nonzero linear relation among the forms ell_ij supported entirely on the zero coordinates q_ij.

To see the graph criterion, regard ell_ij=a_i+b_j as edge incidence vectors of K_(r,s), changing the sign of all column potentials. A set of these vectors is independent exactly when its edge graph is a forest: a forest relation vanishes by removing leaves, while every cycle supplies an alternating relation. Hence a point is singular exactly when its zero-entry bipartite graph contains a cycle.

At a point of L, a zero edge says a_i=-b_j. Thus each nonempty connected component of the zero-edge graph is a complete bipartite block: all its row potentials are equal and all its column potentials are their negatives. A cycle therefore exists exactly when some 2-by-2 block of entries vanishes.

For any fixed two rows and two columns the four zero conditions have rank 3 on L. Their locus in P(L) has dimension n-4=r+s-5, unless r=s=2, in which case it is empty. Because the squaring map is finite and surjective geometrically, for every other r,s>=2 the singular locus has precisely dimension r+s-5 and codimension 3 in X. X_22 is the smooth quadric surface. The other X_(r,s) under discussion are singular. Since complete intersections are Cohen–Macaulay (hence S2) and X is regular in codimension 1, Serre's criterion also proves normality. No discrepancy computation or resolution-type conclusion is being supplied.

## 4. Exact rational-point equivalence for a bare graph

A bare K_(r,s) means r+s distinct positive integer vertices, partitioned into rows and columns, with every cross-part sum a nonzero integer square. It does NOT require that the subgraph be induced; extra within-part square-sum edges are permitted.

Let U be the open subset of X where every q_ij is nonzero, the values q_i1^2 are pairwise distinct in i, and the values q_1j^2 are pairwise distinct in j. Then

U(Q) is nonempty if and only if a bare positive-integer square-sum K_(r,s) exists.

Forward proof: choose rational homogeneous coordinates q. Put a_i=q_i1^2 and b_j=q_1j^2-q_11^2, so q_ij^2=a_i+b_j. All matrix entries are strictly positive. For rational c in the interval

-min_j b_j < c < min_i a_i,

the gauge a_i'=a_i-c, b_j'=b_j+c makes every vertex positive. This interval has width min_i a_i+min_j b_j=min_(i,j) q_ij^2>0. Within-part distinctness is unchanged. A cross-part collision a_i'=b_j' excludes just the one value c=(a_i-b_j)/2, so finitely many forbidden rational choices can be avoided. Choose an integer D divisible by all denominators of the resulting vertices and roots. Then D^2 a_i', D^2 b_j' are distinct positive integers and their cross sums are (D q_ij)^2. The reverse implication simply chooses integer square roots of the cross sums.

This is an exact existence equivalence, not a parametrization producing points of U. Neither geometric dimension nor the canonical twist alone proves U_44(Q) nonempty; the independently verified known example below does. Bare K4,4 existence is established in the literature, and must not be presented as an open problem. It does not settle the stronger decorated magic-square problem.

Concrete gauge check (only K2,2): q=[[1,4],[7,8]], a=(1,49), b=(0,15), c=1/3, D=3 gives row vertices (6,438), column vertices (3,138), and roots [[3,12],[21,24]].

### Cited known point: U_44(Q) is nonempty

Ajai Choudhry, “Two sets of integers such that all elements of the sumset of the two sets are perfect squares,” arXiv:2508.07806v1, 11 August 2025, Introduction immediately before equation (1.1), https://arxiv.org/html/2508.07806v1, records a numerical example from Guy, alongside Jean Lagrange’s parametric family. Take

A=(324,54756,119716,264196), B=(0,79200,227205,1258560).

The exact positive square-root matrix of A_i+B_j is

[18, 282, 477, 1122]
[234, 366, 531, 1146]
[346, 446, 589, 1174]
[514, 586, 701, 1234].

The gauge A'=A-1, B'=B+1 gives positive disjoint vertex sets

A'=(323,54755,119715,264195), B'=(1,79201,227206,1258561).

Every one of the 16 integer square identities has been independently verified, proving U_44(Q) nonempty. This is verification of an existing example, not a new construction. The same source gives known and new parametric families.

This specific point is not on the decorated locus, even after arbitrary independent row and column permutations, common nonzero scaling, and additive gauge. A decorated four-element potential set is a translate of {0,x,y,x+y}, so its sorted outer pair must sum to its sorted inner pair. For A those sums are 264520 and 174472; for B they are 1258560 and 306405. Neither is a parallelogram. Moreover a reflected matching A_i+B_sigma(i)=constant is impossible: the ranges of A and B are 263872 and 1258560, respectively. Exhaustive checking of all 24 permutations confirms zero parallelogram orderings in either part and zero constant reflected matchings. These failures persist under gauge and common nonzero scaling. Thus this certificate does not provide a nondegenerate MSOS point.

## 5. The decorated magic-square locus inside X_44

Order S=(0,x,y,x+y), use the same ordered S for row and column indexing, and impose z_(s,t)=C+t-s. The affine additive-matrix base has dimension 7; this decorated linear subspace has exactly dimension 3 (C,x,y are recovered from z_00,z_0x,z_0y). Its projectivization is P2 inside P6, hence codimension 4. Pullback by the finite squaring map is pure dimension 2, still codimension 4 in X_44. Equivalently it is a complete intersection of 13 quadrics in P15, though this whole pullback is reducible. Indeed, coordinate squaring makes k[q] a free module over k[z], so pulling back the regular sequence of 13 independent linear equations for this three-dimensional vector subspace preserves regularity.

The 16 entries C+t-s comprise exactly nine distinct linear forms

C+epsilon*x+delta*y, epsilon,delta in {-1,0,1}.

Their multiplicities are 4 for C, 2 for each of C±x and C±y, and 1 for each of C±x±y. The nine forms are pairwise nonassociate, so their square classes are independent. For each repeated occurrence, its root equals plus or minus the chosen representative root. There are 16-9=7 independent sign choices, hence exactly 128 geometrically integral components, all defined over Q. Each component is linearly embedded in P8, where its nine root coordinates square to the nine forms above. It is a surface complete intersection of six quadrics, with dualizing sheaf O(3). Distinct components meet on coordinate-zero loci. The full pullback is reduced: its coordinate ring is free over Q[C,x,y], and its generic fiber is a product of 128 separable multiquadratic fields; torsion-freeness excludes nilpotents.

Arrange the nine values as

[C+x,   C-x-y, C+y]
[C-x+y, C,     C+x-y]
[C-y,   C+x+y, C-x].

All rows, columns, and both diagonals sum to 3C. Conversely every rational 3-by-3 magic square has this form after a harmless choice of parameters/order. Thus each sign component is the magic-square-of-squares surface, not the six-dimensional bare X_44. The requirement that all nine square values be nonzero and pairwise distinct removes a finite arrangement of parameter hyperplanes. That requirement is stronger than merely having distinct vertices in the bare graph: one must not drop the decoration or the nine-value distinctness condition when transferring an MSOS conclusion.

## 6. Verified primary-source relationship and limits

Asher Auel and Benjamin Singer, “The algebraic geometry of 3-by-3 magic squares of squares,” arXiv:2609.09351v1, 8 September 2026, https://arxiv.org/html/2609.09351v1 (Introduction and Section 1), study this six-quadric surface in P8. They report 256 A1 singularities and general type for its minimal resolution. Those statements concern the magic-square surface, not X_44. Our component identification above is an explicit coordinate calculation; the bare-model proofs do not rely on importing their classification.

Onno Cain, “Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares,” arXiv:1908.03236v2, 12 August 2019, https://arxiv.org/html/1908.03236v2 (Sections 2–4), studies quartic/factorization formulations and a Gaussian-integer specialization. Those constructions and finite-ring tests do not assert that every bare K4,4 has the decoration above, nor supply a nondegenerate rational point by dimension counting. No use is made here of any stronger unstated converse about Cain's restricted specialization.

All three versioned primary sources were checked on 25 October 2026. The finite and symbolic controls accompanying this note are checks, not a substitute for the proofs, and do not search for or establish a prize solution.
