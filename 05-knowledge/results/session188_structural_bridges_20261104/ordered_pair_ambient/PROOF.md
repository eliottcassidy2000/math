# Ordered-module pair-child separators in every fixed ambient

## 1. Exact ambient theorem

Let `T=Q[D_i]` be a tournament substitution on fixed physical child sets. A distinguished child is exactly `D_z={p,q}`, with `p->q`. Let `U=Q^op[D_i']` have the opposite quotient, the same distinguished pair orientation, and `q_h(D_i')=q_h(D_i)` for every other child and every ordered-cover coordinate h.

After deleting the pair from T, suppose there is an ordered module `M=X->Y`, where X and Y are nonempty disjoint unions of other complete child blocks and `|M|>=3`. All other child interiors and sizes are unrestricted.

Embed the whole core T or U at the root of any fixed ambient tournament K. Outside names, arcs and incidences to the core do not change. Attach physical ports `x->y`. For each exterior name, give both ports exactly the original core-root incidence. Inside the core use four independent Boolean bits

    a_q=s, a_X=u, c_p=t, c_Y=v,

and set every other core bit to zero, with `a_w=[w->x]`, `c_w=[w->y]`. Let L denote the coefficient `[s u t v]` of this restricted physical response.

Write `R_total` for the complete ambient old tournament after deleting p and q. Then

    L(F_ambient(T)-F_ambient(U))
      = C_YX(R_total)-C_XY(R_total)
      >= H(R_total/M)
          [q_2(X)H(Y)+H(X)q_2(Y)+q_2(X)q_2(Y)] > 0.     (1)

Here the two ordered cover slots are nonempty, the first ends in the first indicated set, and the second begins in the second indicated set. At least one of the same sixteen actual finite-difference queries therefore separates in every fixed ambient, including arbitrary strong outer contexts.

## 2. Exact transfer, not a presumed sign lift

For each h, let `P_h^V(a,c)` count all ordered h-path covers of the port-augmented core V. The general master identity in `../quartic_cover_transfer/PROOF.md`, Sections 1--2, applies to the distinguished pair and the saturated grouped observer L. It gives

    L(P_h^T-P_h^U)=L(P_h^T-P_h^W),                      (2)

where W reverses only the pair in T.

For clarity, its mechanism is exact. Reversing every path in a degree-(2,2) contributing cover proves invariance of that entire degree under old-tournament converse. Both ports are internal and nonadjacent in such a contribution, so their fixed mutual edge is unused. After whole-conversing U, its quotient is Q and its pair is reversed. Every other child is uniform to both ports; its complete ordered-cover profile restores it in each h-coordinate. No endpoint-profile assumption is needed on these unmarked children.

The grouped monomial already requires one a incidence at q and a distinct one in X, and one c incidence at p and a distinct one in Y. Hence its original bidegree is saturated. No exterior port bit or additional name inside an aliased group can contaminate it. In particular the prescribed nonzero outside bits do not change the coefficient.

The common augmented core is a module in the ambient. Its exact run substitution is a linear combination of all `P_h` with the same exterior coefficients on both sides. Applying it to (2) identifies the ambient quotient-converse coefficient with the ambient local-pair-flip coefficient.

The pair remains a module globally. So does M after deleting the pair: every added outside name sees the entire core uniformly, and hence all of M uniformly. The local marked-pair identity gives the first equality in (1), now on the full exterior `R_total`. The support-preserving injection and three disjoint missing-cover classes in `../ordered_module_skew/PROOF.md` give its lower bound. Strictness uses `|M|>=3`; no strictness is claimed for two singleton sides.

This argument does not claim that an arbitrary positive core quartic retains its sign in an ambient. It uses the exact all-cover restoration to a particular local pair flip and that local theorem's arbitrary-exterior hypothesis.

## 3. Root-containing original-module families

Assume additionally that Q is prime of order m>=3. Start from T and apply any simultaneous XOR family of its original modules, including the whole core oddly, with every operation contained inside that core. All exterior vertices remain fixed. The original pair and ordered-module hypotheses above remain the source hypotheses.

Every proper original module lies within one original child. The hereditary monotonicity theorem gives coordinatewise nondecrease of the core's ordered-cover profile. If its scalar Hamilton count increases, the internal-zero member of the sixteen-query battery already separates in the ambient. Indeed the augmented core is `x->y->core`, whose profile depends nonnegatively on the original core profile and whose Hamilton count is exactly `H(core)`. The ambient run coefficient of its first coordinate is `H(K)>0`, so the increase remains strict. This is not necessarily the globally all-zero physical mask: the exterior bits retain the core-root incidence.

If the core Hamilton count is equal, hereditary equality restores all child q profiles. If the pair child is selected as well as the root, its orientation is restored. Equation (1) applies. If the pair is unselected, its orientation is conversed; queries uniform on all other children restore the whole-core-converse comparison in every cover coordinate. The compressed exposed-anchor battery at p from `../returns_next/UNIFORM.md` supplies a core contrast `2 modulo4`. Its standard fixed-ambient transfer multiplies by the odd `H(K)`, so it remains nonzero.

Thus the existing fixed quartic/anchor union from `../prime_pair_root_family/PROOF.md`, with the specified outside bits added, separates every such root-containing target in every fixed ambient. It uses at most `m+19` actual queries, improved to `m+18` for odd m, without any new query for the scalar branch. Its choice depends only on the source partition, pair and ordered module, not the selected family. This is source-relative detection, not a common reconstruction battery for all targets.

## 4. Concrete consequences and boundaries

The earlier elementary four-vertex classification gives an eligible ordered module for every five-name quotient with nontrivial children and a pair child. Therefore its quotient-converse theorem now holds in every fixed ambient. For prime five-name quotients, the root-family conclusion holds as well.

The cyclic, staircase and alternating-hub quotient families, and hence all critical prime quotients covered by `../critical_quotient_synthesis/PROOF.md`, satisfy the same ordered-module criterion under their previously stated size hypotheses. Their pair-child quotient-converse and source-relative root-family results likewise hold in every fixed ambient. The sharper local requirement remains `|X|+|Y|>=3`; no arbitrary replacement for that requirement is made.

There is no claim for a quotient without the specified ordered module after pair deletion, for a distinguished child larger than a pair, for operations involving outside names, or for arbitrary-prime reconstruction. This is a separate strengthened corollary; the older packets' deliberately narrower frozen statements are unchanged.

## 5. Verification

The verifier retains full ordered-cover coefficient controls for small cores, sixteen physical source/quotient/local-pair comparisons in strong outer contexts, and complete root-containing child-family controls including crossing intervals of a transitive child. The literal Hamilton engine counts every whole augmented physical tournament, rather than truncating its response. The all-order proof is Sections 1--3; finite controls are not extrapolated to arbitrary ambient size.
