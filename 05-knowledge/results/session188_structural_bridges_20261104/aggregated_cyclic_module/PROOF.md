# Aggregating the cyclic-module quartic removes endpoint balance

## 1. Physical grouped observer

Use the fixed physical conventions of `../ordered_cyclic_module/PROOF.md`: A={p,q,r} is a homogeneous directed triangle p->q->r->p in T. In R=T-A there is a module M=X union Y, with nonempty X->Y and cyclic incidences A->X->Y->A. The exterior E=R-M is arbitrary, as are its uniform incidences toward M and A. The sides X and Y may now be arbitrary nonempty tournaments. No individual endpoint-balance or vertex-transitivity assumption is imposed.

Let T^A reverse only A. Set all port bits to zero except

    a_p=s, a_i=u for every i in X,
    c_q=t, c_j=v for every j in Y.

The actual ports satisfy x->y, with a_w=[w->x] and c_w=[w->y]. Define the grouped quartic contrast

    Omega_XY=[s u t v](F_T-F_(T^A)) under these substitutions.    (1)

**Theorem.** Omega_XY is a strictly positive integer. Its sixteen actual finite-difference queries include a separator. More precisely,

    Omega_XY >= 2 H(X) H(Y) H(E+A_singleton)>0.                 (2)

Here the extra singleton has exactly A's exterior incidences. All statements concern physical grouped queries, not a relabeling of old names or a truncation of the original response.

## 2. Exact alias and run-profile proof

The complete two-port polynomial has bidegree at most (2,2). A monomial contributing to s u t v must contain a_p, exactly one a_i with i in X, c_q, and exactly one c_j with j in Y. No additional a or c factor is possible, and no smaller monomial can acquire four distinct variables through identification. Therefore

    Omega_XY = sum_(i in X,j in Y)
                 [a_p a_i c_q c_j](F_T-F_(T^A)).              (3)

This is an exact alias identity for the full response. In particular it does not discard cubic or quartic terms before substitution.

For h ordered nonempty spanning paths in X, let E_i^X(h) and S_i^X(h) count covers whose one prescribed slot respectively ends or starts at physical i. Every cover has exactly one such endpoint in that slot, hence

    sum_(i in X) E_i^X(h)=sum_(i in X) S_i^X(h)=q_h(X),        (4)

and the same holds for Y. This equality holds without any balance assumption on an individual i. In the exact disjoint-kernel expansion every side contains exactly one physical mark. Summing these marks replaces its complete marked profile by q_h, regardless of whether a term requires a start or an end. There is no extra factor h: the marked slot is prescribed.

Consequently the entire coefficientwise positive generating function K of `../ordered_cyclic_module/PROOF.md`, equation (2), applies to (3) with the functional

    x^s y^t -> q_s(X) q_t(Y),  q_0(X)=q_0(Y)=0,

followed by full E-support extraction. The identity, all eight interfaces, and every return through arbitrary exterior labels are unchanged. This is a use of the established ordered-module positive kernel, not a new claim of an identity beyond it. Its xy term gives (2), since q_1=H and every tournament has a Hamilton path.

The result also covers singleton sides and empty E. The transitive attachment X->A->Y is outside the cyclic hypothesis and is not included.

## 3. Quotient converse with arbitrary children

Suppose A,X,Y are children of a fixed quotient Q with the preceding cyclic ordered-module witness after deleting A. Let T=Q[A,D_1,...] and U=Q^op[A,D_1,...], holding all child interiors and physical names fixed. The same grouped quartic is positive for F_T-F_U.

To see this, whole-converse U. The degree-four part of the complete response is invariant under whole converse. Every child except A is homogeneous to both actual ports in each of the sixteen grouped queries, including X and Y. Its entire ordered-cover profile therefore suffices for exact substitution. Converse preserves q, so all those child interiors can be restored, leaving exactly the local A-converse contrast in (1). Equivalently one can apply (4) in the inherited quartic restoration proof.

This restoration continues to hold if the target non-A children have the same complete q profiles as the source, with A restored to its original orientation after quotient converse. No physical one-marked invariance of an individual child is needed. For the sixteen grouped queries, the q restoration preserves the entire actual count on each side of the whole-converse comparison, not merely one selected coefficient.

These quotient comparisons are whole-core statements. No arbitrary strong outer-context lift follows from positivity alone. The local theorem still allows any exterior satisfying its two module hypotheses.

## 4. All original root-containing families

Let Q be prime of order m>=3, and let z,u,v be distinct quotient names with z->u->v->z and {u,v} a module of Q-z. Inflate z by C3 and every other quotient name by an arbitrary nonempty tournament on a fixed physical child set. Then every family of original modules containing the whole root oddly changes the full named response. A fixed battery of at most m+20 actual queries, depending only on the source, detects all such targets.

Proper original modules lie inside individual children. Under simultaneous XOR reversal of a family, the quotient is conversed. The inherited original-module profile-monotonicity theorem gives H(target)>=H(source). If strict, query (0,0) separates; it is already among the sixteen grouped quartic queries.

Assume scalar equality. The same theorem's hereditary equality consequence gives equality of every corresponding child q profile. If the triangle child is selected, its orientation is restored after the root converse. Section 3 restores all other children at the grouped queries, so the strictly positive grouped quartic separates.

If the triangle child is unselected, its orientation is conversed. Split its three names into singletons and keep every other child one uniform block. Equality of child q profiles restores the target to the whole-converse response of the source. The exposed-singleton anchor theorem in `../returns_next/UNIFORM.md` then gives an actual separator. There are m+2 compressed blocks and at most m+4 anchor queries. Taking their union with the sixteen grouped queries gives m+20. This includes crossing and nested original-module families inside arbitrary children. It does not assume those modules remain modules during a sequential execution.

This strengthens the corresponding root-family theorem in the ordered-cyclic-module packet by removing vertex transitivity of the two witness children. The older physical-mark theorem remains useful when one wants fixed individual marked names rather than aggregate port blocks.

## 5. A whole-child orientation cube

In the setting of Section 4, allow the quotient orientation Q or Q^op and independently allow each child to retain its interior or take its whole converse. Let k be the number of children of order at least two. These 2^(k+1) physically distinct tournaments have pairwise distinct full named responses.

Compare two states. If their quotient orientations agree, their differences are whole converses of disjoint children; the established disjoint-module isolation theorem separates them. If the quotient orientations differ, relative to the first state their difference is a root-containing family of its whole children. The triangle remains C3 in either state. When the quotient is reversed, exchange the witness roles u and v to recover z->u->v->z. Section 4 applies from that source state, with no physical response relabeling. This proves pairwise distinction.

If every child is nontrivial, k=m and the cube has 2^(m+1) states. Singleton choices do not create distinct states and are not counted twice. The cube conclusion concerns root and whole-child converses. For arbitrary children it is NOT a theorem of injectivity for every possible internal original-module XOR state. Nor is m+20 claimed to be a common reconstruction battery across this entire cube.

Every concrete critical T,U,W quotient has the required deleted-root cyclic pair, as established in the ordered-cyclic-module packet. Thus the theorem applies to those quotients when the child at the designated cyclic-witness root z is C3 and the other children are arbitrary, at every order. For the displayed W witness z is the hub; for U it is an A_i. A C3 child at an arbitrary different quotient name is not asserted to supply this witness. It is not arbitrary-prime closure or unrestricted reconstruction.

## 6. Reproduction and limits

The standalone checker retains exact grouped finite-difference values, root-containing families including a crossing TT3-child control, a complete whole-child orientation cube, and a complete named polynomial bank with cubic/quartic coordinates. The reused run identity is checked over integers. Direct literal Hamilton DP checks the physical augmented graphs independently of the run functional.

The all-order conclusions are proved above. The finite data are controls, not a proof by enumeration at all orders. No numerical sign transfer to an arbitrary strong ambient context is asserted.
