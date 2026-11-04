# Ordered-cover quartics transfer the cyclic-module separator through every ambient

## 1. A full-cover converse identity

For a named tournament T and actual ports x->y, let P_h^T(a,c) be the number of ordered covers of the augmented tournament by h nonempty directed paths. As usual a_v=[v->x] and c_v=[v->y]. This is the complete Boolean polynomial in all physical attachment bits.

**Lemma 1.** Each P_h has bidegree at most (2,2). Its entire bidegree-(2,2) part is invariant under whole converse of the old tournament:

    [a_i a_j c_k c_l] P_h^T = [a_i a_j c_k c_l] P_h^(T^op)

for distinct i,j and distinct k,l, allowing overlap between the two pairs.

Proof. A cover uses at most one incoming and one outgoing arc at each port. Thus it has at most two a factors and two c factors. A degree-(2,2) contribution requires both ports to be internal in their respective cover paths, with two old-name contacts each. The ports cannot be adjacent: using their edge leaves at most one old-name contact at each port.

For a fixed ordered cover word skeleton, the two contact factors at x are a_i(1-a_j), whose degree-two part is -a_i a_j; likewise the y contact contributes -c_k c_l. Their product has sign +1. Reverse every path, keeping its ordered slot and the physical names of both ports. Every old arc is now directed in T^op. The port contact directions are exchanged, but their top-degree factors are the same. The edge x-y is not used, so its opposite direction after a literal global converse is irrelevant. This is an involutive bijection of the degree-(2,2) contributing skeletons and proves equality for every h. In particular no Hamilton-only assumption is used.

**Lemma 2.** Fix a coefficient using two a variables and two c variables on selected old names, possibly after identifying a specified group of names with one bit. If it already needs two original a incidences and two original c incidences, then fixing any other attachment bits cannot change it.

Indeed an extra factor from an unselected name would exceed one of the bidegree bounds. Boolean identification does not create a route around this: in the grouped observer below a_p is distinct from a_X, so both already require separate original a factors; c_q and c_Y likewise. This explains why nonzero exterior lift bits do not contaminate the quartic.

## 2. Fixed ambient and exact quotient restoration

**General transfer identity.** Let T=Q[A,D_1,...,D_k] with an arbitrary nonempty distinguished child A, and let U=Q^op[A,B_1,...,B_k] on the same physical child sets. Suppose q_h(B_i)=q_h(D_i) for every i,h. Let W converse only A in T. Take any saturated bidegree-(2,2) coefficient observer whose attachment bits are uniform on each non-A child, while the attachment bits inside A may remain individually named. Embed the whole core in any fixed rooted ambient and give both ports the core's exterior incidences. Then the source-minus-U observer equals the source-minus-W observer exactly in that ambient.

Indeed, for every ordered-cover coordinate h, Lemma 1 permits whole-conversing U in that coefficient. Its quotient becomes Q, A becomes A^op, and every other child becomes B_i^op. Those children remain homogeneous to both ports throughout the observer's actual Boolean evaluations. Their complete q profiles equal those of D_i, so exact run-cover substitution restores them in every h coordinate. Thus the coefficient agrees with W for every h. The common ambient run expansion is a linear combination of those coordinates, proving the assertion. Lemma 2 accounts for the fixed nonzero exterior bits. No triangle, ordered-pair, cyclic-incidence, positivity, or parity assumption is needed for this identity itself.

We now apply this master identity to the positive cyclic-module observer. Its structural hypotheses enter only the strict lower bound, not the restoration identity.

Let the core T=Q[A,D_1,...] have A=C3 and the cyclic ordered-module witness from `../aggregated_cyclic_module/PROOF.md`: A->X->Y->A, and X union Y is a module after deleting A. The other children and the remaining core exterior are arbitrary. Let U=Q^op[A,D_1,...], with all child interiors and physical names held fixed. Put W=T with only A conversed.

Let K be any fixed rooted ambient quotient. Replace its distinguished root by the core N. Thus K denotes the collapsed quotient, not the expanded tournament. The exterior names, arcs and incidences to N are identical in the source and target. For each query give BOTH ports N's exterior incidences, and use the grouped internal bits

    a_p=s, a_X=u, c_q=t, c_Y=v,

with all other core bits zero. These are sixteen actual queries on fixed physical names. Write L for the coefficient [s u t v].

The augmented core N union {x,y} remains a module. The exact ambient run expansion is therefore

    F_ambient = sum_(h>=1) d_h P_h^core,  d_h>=0,

with the same fixed coefficients d_h for T,U,W. They depend only on the unchanged collapsed ambient. No congruence or sign transfer is used here.

For each h, Lemma 1 replaces the quartic part of P_h^U by that of P_h^(U^op). In U^op the quotient is Q and all children are conversed, including A. Every other child is uniform to both ports at all sixteen grouped queries. Its complete q profile is converse-invariant, so restoring it leaves every ordered-cover count unchanged. The only remaining change is A converse. Thus

    L(P_h^T-P_h^U)=L(P_h^T-P_h^W) for every h,

and consequently

    L(F_ambient(T)-F_ambient(U))
       = L(F_ambient(T)-F_ambient(W)).                      (1)

The same restoration applies if target non-A children have the same complete q profiles as their source counterparts and A is restored after quotient converse. The equal profiles, rather than literal equality of child interiors, are what is required.

Now apply the arbitrary-exterior local positivity theorem to the entire ambient tournament. A remains a homogeneous triangle, and X union Y remains a module after deleting A: every added exterior name sees all of N uniformly. Lemma 2 shows that the nonzero prescribed exterior port bits do not alter this saturated quartic. If E_total consists of every ambient old name outside A union X union Y, equation (1) has the strict lower bound

    L(F_ambient(T)-F_ambient(U))
       >= 2 H(X) H(Y) H(E_total+A_singleton)>0.             (2)

**Theorem.** The grouped cyclic-module quartic separates the quotient-converse comparison in every fixed ambient tournament, including an arbitrary strongly connected outer context. At least one of the same sixteen lifted actual queries separates.

This is a specific exact quartic identity. It does not assert that an arbitrary positive core query gap, or an arbitrary positive core quartic, keeps its sign under ambient substitution. The hypotheses making the local kernel positive remain essential.

## 3. Physical balanced-mark version

The same ambient conclusion holds for the individual physical coefficient [a_p a_i c_q c_j] in `../ordered_cyclic_module/PROOF.md` when i in X and j in Y have complete balanced one-marked profiles with positive Hamilton entries. The only change in Section 2 is restoration: each marked child has one prescribed marked endpoint slot among its run slots, so every P_h uses its complete E_i or S_i profile. Converse exchanges these equal profiles. Unmarked children still use q. Lemma 1 and the common ambient expansion then reduce to the local arbitrary-exterior physical-mark theorem exactly as above.

This does not replace individual endpoint balance by an unsupported factorization of joint marks in one child. There is exactly one mark in each of the two exterior marked children; the two triangle marks are handled by the exact C3 kernel.

## 4. Root families and orientation cubes in arbitrary fixed ambients

Use the prime quotient and designated cyclic-witness root z hypotheses of Section 4 of `../aggregated_cyclic_module/PROOF.md`. The z child is C3; all others are arbitrary nonempty tournaments. All selected original-module operations are inside the core. The exterior is fixed.

Every root-containing original-module XOR family is detected by at most m+20 fixed lifted actual queries, exactly the grouped-quartic/anchor union from that packet. If the core scalar Hamilton count is unchanged, hereditary equality supplies the child q restoration. A selected triangle child is handled by Section 2; an unselected one by the inherited exposed-anchor separator and its fixed-ambient modulo-four lift.

For completeness the scalar-strict branch needs no added query. The internal all-zero member of the grouped battery makes the augmented core the ordered sum x->y->core. Every q coordinate of this sum is a nonnegative linear combination of the core q coordinates, and its Hamilton count is exactly H(core). Original-module q monotonicity therefore gives coordinatewise monotonicity of the augmented core profiles, with a strict first coordinate when H(core) increases. In the fixed ambient run expansion d_1=H(K)>0. Thus this same lifted all-zero internal query strictly increases. It need not be the globally all-zero attachment mask; its exterior bits remain the core incidences.

The root-plus-whole-child orientation cube consequently has 2^(k+1) distinct full named responses in every fixed ambient, where k is the number of nontrivial children. If the quotients agree, use disjoint-module isolation; if they differ, use the root-family result from the chosen source state. This is not the full internal original-module orbit of arbitrary children. The m+20 detection battery is not asserted to be one common reconstruction battery for the cube.

All concrete critical T/U/W quotients have the needed witness. The C3 child must be at a valid cyclic-witness root, for example the W hub or an A_i in the displayed U convention; an arbitrary C3 location is not enough for this theorem.

## 5. Scope and reproduction

Lemma 1 is an all-order ordered-cover identity. The subsequent strict ambient conclusion uses the special cyclic ordered-module kernel. It strengthens the earlier packets' deliberately limited whole-core scope, without modifying those frozen proofs. It does not give arbitrary-prime closure or unrestricted reconstruction.

The verifier retains complete ordered-cover coefficient controls, exact sixteen-query source/quotient/local contrasts in strong ambients, physical balanced-mark controls, complete root-family comparisons including scalar-strict crossing operations, and an entire small orientation cube. All integers and named masks are retained. Finite controls complement the proofs, not replace them.
