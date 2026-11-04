# A sign-definite grouped-skew polynomial criterion

## 1. Contribution and precise scope

This packet gives an all-order sufficient criterion for a named mixed-quartic battery to separate an unchanged odd transitive child under a source-relative quotient converse, in every fixed ambient. The criterion is coefficientwise sign of a finite polynomial on the physical contracted core. It does not require an ordered witness pair, equal profiles between two other children, or log-concavity.

Four exact noncirculant prime-nine examples satisfy the criterion when every other child is TT2. This gives unbounded odd-transitive-child families and arbitrary fixed ambients by proof, not finite extrapolation. The same four quotient examples have subsequently received broader provisional arbitrary-child affine certificates in the parallel prime-nine work. This packet therefore claims a distinct structural criterion and its exact certificates, not exclusive credit for closing those quotients or unrestricted prime-nine closure.

Generic local module-converse detectability is already covered by incoming E4. The contribution here is an explicit restricted mixed-quartic certificate with source-relative root and fixed-ambient scope. The companion kernel-span packet concerns a different, size-dependent local restricted battery; it is not needed for this theorem.

## 2. The physical contracted-core polynomial

Let the source core be C=Q[A,D_1,...,D_s], where A is transitive of odd order m>=3 and all children are nonempty. Let G be one complete child D_j, or a specified union of complete children. Contract A to a singleton z without contracting any other physical name. Write the resulting tournament as V=Q[z,D_1,...,D_s], of order N. Every exterior name of a fixed ambient sees this whole core uniformly.

In the squarefree algebra on the physical names of V, let M(t) be the simple-path matrix with edge weight [u->v]+t. Singleton paths are included. Write C_i(t) for its column sum and R_i(t) for its row sum, and sum these over i in G to obtain C_G,R_G. Define the full polynomial

    k_G(t)=[V]{C_z(t)R_G(t)-C_G(t)R_z(t)}.              (1)

The two factors in each term have disjoint support. Thus the coefficient of t^r is a signed joint two-cover count, not a product of separately evaluated endpoint marginals. Its degree is at most N-2. In particular,

    k_G(0)=sum_(i in G)(C_zi(V)-C_iz(V)),              (2)

where C_uv(V) is the ordered spanning two-path-cover capacity ending at u and starting at v, with the gap arc untested. All names of every uncontracted child remain present in (1). Replacing them by unweighted quotient vertices would generally change the polynomial.

## 3. Exact transfer to every fixed ambient

Let E be any fixed tournament exterior to V, with each of its vertices uniformly adjacent to V. Let alpha,beta,tau be its usual prefix, suffix, and excursion series:

    alpha=1+sum_(v->V) E_v(E),
    beta =1+sum_(V->u) S_u(E),
    tau  =sum_(V->u,v->V) G_uv(E).

Here all products remain squarefree on E. Eliminating maximal exterior runs in a path gives full endpoints C_i(full)=alpha C_i(tau) and R_i(full)=beta R_i(tau), for i in V. Consequently the actual grouped skew in V+E is

    K_G(V+E)=[E]alpha beta k_G(tau)
            =sum_(r>=0) theta_r [t^r]k_G(t),
    theta_r=[E]alpha beta tau^r >=0.                    (3)

Moreover theta_0=H(E+v)>0, where v is one virtual vertex with the core's incidences. Indeed, a Hamilton path of E+v is uniquely a permitted prefix, v, and permitted suffix. The equality remains valid for empty E, when theta_0=1.

**Theorem 1 (sign-definite transfer).** If k_G(0) is nonzero and all coefficients of k_G have its sign, then K_G(V+E) has that same strict sign in every fixed ambient E.

**Proof.** Equation (3) is a nonnegative linear combination, with a strictly positive weight on the nonzero constant coefficient. QED.

Zeros among the coefficients are permitted. The hypothesis is sufficient, not claimed necessary; a mixed-sign polynomial does not establish cancellation, and failure of the criterion is not a response collision.

## 4. The actual mixed-quartic battery

Use the reviewed odd-transitive identity in ../session189-odd-transitive/PROOF.md, Theorem 1. Number the physical A order as 0->1->...->m-1, put m=2k+1, and define

    B_G=Delta_(01;m-1,G)
          +sum_(r=2)^k (-1)^r binom(m-2,r-1) Delta_(01;r,G),
    Delta_(01;r,G)=sum_(i in G)[a_0 a_1 c_r c_i](F_T-F_(local A converse)).

That exact identity yields

    B_G=2 K_G(V+E).                                    (4)

Each term is implemented by making c uniform on G and taking the coefficient of a_0 a_1 c_r c_G. Bidegree at most(2,2) makes this coefficient exactly the sum of the indicated complete named quartics: with two a and two c variables already used, no higher alias can occur. Each quartic has a legitimate sixteen-query Boolean finite-difference extraction. The binomial battery therefore uses at most8(m-1) actual queries. The physical group must not overlap A.

Combining (3)--(4), the coefficientwise sign hypothesis gives a nonzero binomial combination for every odd m>=3 and every fixed ambient, independent of the size of that ambient. The size of A does not occur in k_G.

## 5. Quotient-only and original-module root-family scope

First consider C'=Q^op[A,D'_1,...,D'_s], with A unchanged in its physical order and each D'_i having exactly the same complete ordered-cover profile as its source D_i. The target is a quotient-only converse with profile-preserving interior changes; it is not the local A converse.

For every query used above, both ports are uniform on every D_i. The reviewed all-cover uniform-child restoration therefore replaces D'_i by D_i^op without changing any augmented-cover coordinate. The restored target is the whole-core converse of Q[A^op,D_1,...,D_s]. The inherited top-quartic whole-core-converse invariance then identifies the quotient-target quartic contrast with the local A-converse contrast in every augmented-cover coordinate. This is the same precise restoration used in the reviewed odd-transitive root theorem, but its nonvanishing input is now (3), not an ordered-witness skew inequality.

In the actual fixed ambient, both ports imitate the core at every exterior name e: a_e=c_e=[e->C]. The all-cover restoration lifts through that ambient. Equations (3)--(4) therefore prove the quotient-target separation with these physical masks; they do not silently converse any ambient arc.

**Theorem 2 (source-relative root-family criterion).** Let Q be prime, let A=TT_m with odd m>=5, and suppose some complete-child group G satisfies Theorem 1 on V. Then every simultaneous XOR family of original source modules selecting the whole core oddly is separated from the source by a fixed finite two-port battery in every fixed ambient.

**Proof.** Use the inherited scalar/profile monotonicity reduction: a strict scalar increase is detected by the core-zero query with the ports imitating the core outside. In the remaining branch all corresponding child profiles agree. The reviewed internal transitive-distance quartics force the physical target A to be A or A^op. The latter branch is separated by the general conversed-child compressed-anchor lemma. In the unchanged branch, the preceding all-cover restoration applies and (3)--(4) give a strict binomial mixed contrast. The union of these batteries depends only on the source and ambient, not on the chosen original-module family. QED.

One conservative query bound is

    1 +16*m*binom(m-1,2)+8*(m-1)+(m+|Q|+1).

The four terms are the scalar query, internal-order tomography, one grouped binomial battery, and the compressed-anchor battery. This is source-relative detection, not pairwise injectivity among arbitrary outputs. The m>=5 condition is used for the inherited internal-order reconstruction; the unchanged-A mixed identity itself holds for odd m>=3.

## 6. Four exact prime-nine certificates

Use z=0 and the following outgoing-row quotients, with every other physical child TT2:

    Q860 =(192,289,99,71,15,153,306,94,189), G=child1;
    Q2005=(320,193,83,71,43,15,432,61,190), G=child4;
    Q4067=(448,209,195,71,45,15,432,312,62), G=child4;
    Q6501=(384,177,147,135,105,77,15,368,126), G=child1.

The numeric suffixes record their incoming census node IDs, but no completeness assertion about that census is needed here: the displayed rows are self-contained certificates. All four are prime by direct examination of every proper nonsingleton subset. Their degree sequences are nonregular, hence they are not circulant tournaments. Direct exact endpoint and full repeated-run calculations give p=sigma=0 at all nine vertices. There is no ordered witness pair in Q-z of the type required by the earlier odd-transitive ordered-witness theorem.

The contracted physical core V has17 names, not9. Each selected k_G has every coefficient nonnegative and has constant coefficient, respectively,

    175478400, 158246148, 248375296, 50196404.             (5)

All16 coefficients of every one of the nine grouped skew polynomials are retained in selected.out, including zeros. The four selected complete coefficient lists are also printed in receipt.json. They are independently reproduced by two different exact methods:

1. weighted.cpp keeps polynomial endpoint counts on every physical subset of the17-name V, and convolves complementary subsets to compute (1);
2. runword.cpp compresses only the eight TT2 modules. A run word visits each such child once or twice, with ordered-cover weights q_1(t)=1+2t and q_2(t)=2. Its quotient transitions have weight [u->v]+t. Perturb the z->G edge by +epsilon and the G->z edge by -epsilon; the coefficient of epsilon in the Hamilton path polynomial is exactly k_G(t). This independent derivative recurrence agrees coefficient by coefficient.

Both computations are exact integers. weighted.cpp uses 64-bit endpoint arrays, justified by the bound16!*max_r binom(16,r)<2^63 at17 names, and128-bit products/accumulators. The compressed engines use128-bit integers. No floating-point sign test is used.

For each displayed quotient, Theorem 2 now covers Q[TT_m,TT2,...,TT2] for every odd m>=5, every source order labeling of those transitive children, every original-module root-containing family, and every fixed ambient. This conclusion is proved from a finite polynomial certificate; it is not inferred from testing several m values.

## 7. Physical and complete named controls

The208 retained actual queries cover the four m=3 unchanged-A quotient comparisons, one m=5 comparison, and a genuinely cyclic two-name ambient for the first quotient. Every query records both physical masks, the complete physical augmented tournament, and the compressed TT2 representation. Source and quotient-target counts remain distinct.

The bare binomial contrasts are350956800,316492296,496750592,100392808. For the first quotient the m=5 contrast remains350956800. The cyclic ambient has theta_0=3, theta_1=1 and gives10037256456=2(3k_0+k_1). The nonzero outside mask explicitly sets both port incidences at the incoming exterior name. These controls use the full physical graph, not a singleton replacement for a TT2 child.

count_physical.cpp evaluates actual Hamilton counts by the exact run-word decomposition, with q_1=1,q_2=2 for each TT2. An additional100 independently recounted physical inflated tournaments check this compression directly by a predecessor subset recurrence.

Nine smaller source, quotient-target and local-A-converse graphs retain15,219 complete named response coefficients, including every cubic, quartic and zero. Their actual values are independently recounted by predecessor.cpp, then Mobius-inverted again. They verify the ambient polynomial transfer and all-cover restoration at the response level. The large prime-nine examples retain all grouped excursion-polynomial coefficients and the specified actual query masks/counts, not an uncomputed full named response bank; no claim of such a bank is made.

## 8. Sources and remaining boundaries

Dependencies are the reviewed Session189 transitive-distance theorem, general conversed-child anchor integration, and odd-transitive exact mixed identity/restoration. Their precise statements and inherited S188 source links are in the sibling proof packets; this note does not strengthen them by assumption. Incoming generic local-converse E4 is pinned at commit76b3b1a77d124100f03026420f6b821533cdec15, path05-knowledge/results/two_port_parity_separation_20261104/PROOF.md.

No full named-response collision, unrestricted reconstruction theorem, arbitrary-child prime-nine theorem, or new global frontier is asserted here. If a k_G polynomial has coefficients of both signs, this sufficient criterion simply does not decide that group. The maximum-matching local theorem in the separate kernel-span packet must not be substituted for (3) without its physical-exterior and admissible-mark restrictions.
