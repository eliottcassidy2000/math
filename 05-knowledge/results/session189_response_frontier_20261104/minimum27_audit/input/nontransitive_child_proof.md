# Positive mixed quartics for crossed TT4 root targets

## 1. Statement and boundary

Let C=Q[A,D_1,...,D_s], where s>=1, A has physical transitive order0->1->2->3, and every D_i is nonempty. Consider the physical target

    C'=Q^op[A',D'_1,...,D'_s],

where corresponding D_i,D'_i have equal complete ordered-cover profiles. Suppose the two ends of A' are exactly the source's two middle names1,2. There are precisely four such orders:

    (1,0,3,2), (1,3,0,2), (2,0,3,1), (2,3,0,1).          (1)

**Theorem.** A fixed battery of named mixed-quartic queries separates C from every target in(1), in every fixed ambient. Both ports imitate C at every ambient name. The battery is uniform on every other child, and can use one common bit on their entire union. The result requires no restriction on Q beyond having another child, and no log-concavity or equality between different source children.

The theorem closes the crossed-TT4 branch left by the incoming common-end result. It does not use that incoming result in its proof. A statement excluding **all** transitive children combines this crossed-branch theorem with the separately reviewed common-end theorem, plus the established source-family reductions. That integration awaits independent review of the present crossed-branch packet.

The targets in(1) are physical orders, not arbitrary relabelings of the whole core. The source, the quotient target, the restored whole-core converse, and the local replacement below remain distinct.

## 2. General two-kernel replacement identity

This section applies to arbitrary tournaments A and B on the same physical name set, not only transitive ones. Let A be a module in T with exterior R. Replace only its interior by B, preserving every exterior arc and every name, to obtain T_B. Fix distinct p,q,r in A and i in R.

Use the squarefree path algebra of the reviewed odd-transitive packet. On R set

    alpha=1+sum_(v->A) E_v,  beta=1+sum_(A->u) S_u,
    tau=sum_(A->u,v->A) G_uv,
    ell_i=sum_(v->A) G_iv,  rho_i=sum_(A->u) G_ui.

For each internal tournament Z=A or B, let M_Z(t) be its simple-path matrix with edge weights [u->v in Z]+t. Column sums are C_j, row sums are R_j, and H is their total. All extraction below is joint squarefree extraction on A's physical set. Define

    L_Z(t)=[A](C_p M_qr+C_q M_pr),
    R_Z(t)=[A](M_rp R_q+M_rq R_p),
    P(t)=L_A(t)-L_B(t),  Q(t)=R_A(t)-R_B(t).           (2)

Here R_Z(t) is a scalar polynomial, distinct from the exterior set R. Then the complete named quartic difference is

    [a_p a_q c_r c_i](F_T-F_(T_B))
       =[R]{alpha S_i P(tau)+beta E_i Q(tau)
                 +alpha beta(ell_i P'(tau)+rho_i Q'(tau))}.    (3)

**Proof.** The full eight-term quartic expansion into ordered three-path covers, after exterior elimination, is

    alpha S_i L_Z +beta E_i R_Z
      +alpha beta ell_i [A]{H(C_p M_qr+C_q M_pr)
                            +C_r(C_p R_q+C_q R_p)}
      +alpha beta rho_i [A]{H(M_rp R_q+M_rq R_p)
                            +R_r(C_p R_q+C_q R_p)},           (4)

with each internal expression evaluated at t=tau. This is the same complete expansion used in the reviewed two-channel converse identity, before imposing any converse relation.

Differentiating the weighted path resolvent gives M'_uv=R_u C_v, C'_j=H C_j, and R'_j=H R_j. Therefore the two bracketed expressions in(4) are L'_Z and R'_Z respectively. Subtract the formulas for A and B. This proves(3). No equality of internal endpoint marginals or even of internal profiles was used. Each product keeps all shared-name support restrictions. QED.

The older local-converse identity is the special case Q=-P. Arbitrary physical replacement needs both independent kernels; replacing Q by -P here would lose the positive channel that proves this theorem.

## 3. Four complete TT4 kernel certificates

For each quotient target order A' in(1), put B=(A')^op, the reversed target order on the same four names. For the moment compare the source A=(0,1,2,3) with this local B. Write

    D_(pq;r,i)=[a_p a_q c_r c_i](F_T-F_(T_B)).

The following combinations have equal constant kernels P=Q=gamma, and have P'=Q'=0:

- A'=(1,0,3,2): D_(01;3,i)-D_(01;2,i), gamma=2.
- A'=(1,3,0,2): 2D_(01;3,i)-D_(01;2,i), gamma=2.
- A'=(2,0,3,1): 2D_(01;3,i)+3D_(01;2,i), gamma=-2.
- A'=(2,3,0,1): D_(02;3,i)-D_(02;1,i), gamma=1.

These are exact degree-two polynomial identities. A standalone finite proof is to enumerate all24 permutations of the four physical names, cut each at each of its three internal positions, and expand its two remaining tested edges with weights [u->v]+t. L counts covers whose first path ends at p and whose second starts at q and ends at r, plus the p/q exchange. R counts covers whose first starts at r and ends at p and whose second starts at q, plus the p/q exchange. These are precisely(2). inner_kernels.json retains all12 marked triples for each of the five internal orders, both L and R, and every coefficient including zero. receipt.json additionally lists every difference kernel in the four certificates.

Inserting each constant certificate into(3) gives exactly

    gamma [R](alpha S_i+beta E_i).                     (5)

No excursion term remains, and there is no asymptotic approximation or parity-only claim.

## 4. Positive symmetrized capacity and arbitrary ambient

Contract A to a singleton z, leaving every other physical name present. For any exterior name i, the two ordered spanning two-cover capacities satisfy

    C_zi=[R]alpha S_i,   C_iz=[R]beta E_i.

Thus(5) is gamma(C_zi+C_iz), a symmetrized capacity rather than a skew. Let G be the union of all other core children. Summing over i in G is implemented by giving c one uniform bit on G; all other children therefore remain uniform in every query. Bidegree at most(2,2) proves that the coefficient a_p a_q c_r c_G is exactly the sum of the named quartics. Higher aliases cannot contribute to a saturated quartic.

Let V=Q[z,D_1,...,D_s] be the physical contracted core, with at least two vertices. Define

    s_G(t)=[V]{C_z(t)R_G(t)+C_G(t)R_z(t)},             (6)

again with every internal edge weighted [u->v]+t. Every coefficient of this polynomial is nonnegative, since each is an ordinary weighted two-cover count with a plus sign.

Its constant coefficient is strictly positive. Take any Hamilton path of V. Because V has at least two names, some arc of that path is incident with z; its other endpoint lies in G. Delete that arc. The resulting ordered two-path cover contributes to C_zi if the arc is z->i, or to C_iz if it is i->z. More quantitatively,

    s_G(0) >= 2H(V)-S_z(V)-E_z(V) >= H(V)>0,           (7)

because a Hamilton path cannot both start and end at z when |V|>=2.

Now put the core in any fixed ambient E. Its maximal exterior-run decomposition gives

    C_zG(V+E)+C_Gz(V+E)
       =sum_(r>=0) theta_r [t^r]s_G(t),
    theta_r=[E]alpha_E beta_E tau_E^r>=0,
    theta_0=H(E+v)>0.                                 (8)

This is the exact positive excursion transfer, not a product of separately evaluated marginals. Equations(6)--(8) make the grouped symmetrized capacity strictly positive in every fixed ambient. They do not require that each individual name i be a separator. Using all of G guarantees positivity without any claim about a prescribed individual endpoint.

Consequently each crossed target has a nonzero grouped certificate with its displayed sign gamma in every ambient. Each certificate uses two grouped quartics, so at most32 actual queries. Their union uses the four marked triples(01;2),(01;3),(02;1),(02;3), hence at most64 queries before deduplication, independent of the target order or ambient size.

## 5. Actual quotient-target restoration

The local replacement B above is **not** the actual target A'. For each actual query, both ports are uniform on every other child D_i. By equal complete ordered-cover profiles, the inherited all-cover uniform-child restoration replaces D'_i by D_i^op in every augmented-cover coordinate. The resulting quotient target is

    Q^op[A',D_1^op,...,D_s^op]
       = (Q[B,D_1,...,D_s])^op,   B=(A')^op.           (9)

The inherited top-quartic whole-core-converse invariance identifies the quartic of the right side with the quartic of the local-replacement core Q[B,D_1,...,D_s]. This invariance holds in every augmented-cover coordinate, exactly as in the reviewed odd-transitive restoration. Thus the grouped source-minus-quotient-target quartic equals the local source-minus-B quartic used in(3)--(5).

At each fixed ambient vertex e, set a_e=c_e=[e->C]. These are physical bits and are fixed throughout the finite differences. The all-cover restoration in(9) lifts with these imitating bits. It never reverses or relabels an ambient arc. Equations(8) and(9) together establish the theorem in every fixed ambient.

The restoration does not assert equality of individual quartics whose exterior mark singles out one name inside a multi-name D_i. The proof and controls compare the **whole uniform group** G at that step. The general local-replacement identity(3), by contrast, is valid individually and is checked individually.

## 6. Original-module consequence and incoming dependence

Consider a simultaneous XOR family of original source modules at a highest changed prime node having a TT4 child A. In the scalar-increase branch the inherited scalar battery already separates. Otherwise hereditary equality gives equal child profiles and keeps the target of A transitive. If its two target ends are exactly the source's two middle names, the preceding theorem separates it. This crossed-branch consequence uses the established highest-node/profile reduction and ambient restoration; it does not use incoming Theorem B or E4.

The remaining20 of the24 target transitive orders share at least one source end as a target end. The independent incoming-anchor audit approves the common-end theorem with the necessary quotient-order condition m>=2; it separates precisely this complementary case. Combining that corrected version with an independently approved crossed theorem would therefore remove the TT4 exception. Together with incoming TT2/TT3 separation and reviewed distance reconstruction for order>=5, it would exclude all transitive children at a surviving highest changed prime node. This is a clearly conditional integration statement until the crossed-branch review is complete, not an additional unconditional theorem claimed by this packet.

Incoming E4 already handles generic local **module converse** detectability. It is neither needed nor misused for the arbitrary replacements B in(2), and it does not replace the quotient/ambient restoration in(9). The new contribution is the explicit positive mixed-quartic source-root certificate for the crossed physical TT4 targets.

## 7. Complete controls and reproduction

The packet retains36 complete named response polynomials: four contexts, each with its source, four physical quotient targets, and four local replacements. The contexts are C3[TT4,1,1], C3[TT4,TT2,TT2], the latter core in a cyclic two-name ambient, and the prime-five carousel with a TT4 child and four singleton children. All old orders are at most10.

Every cubic and quartic coefficient, including zero, is kept. An independent predecessor recurrence recounts all actual small-support values, and the verifier redoes the full Mobius inversion. The general replacement identity(3) is checked for all768 individual mixed quartic differences, including marks in the fixed ambient; quotient restoration is checked only on the appropriate uniform core group.

The complete contracted capacity matrices have positive grouped symmetrized capacities4,24,216,36 respectively. Thus the four signed gaps in each context are gamma times that value. The cyclic ambient is nontrivial and retains its nonzero exterior port masks.

Finally,544 distinct actual grouped finite-difference queries are independently counted on their full physical augmented tournaments, with source and target rows/masks recorded. These are literal Hamilton counts, not values inferred from the kernel certificate. Normal and optimized verifiers must reproduce identical receipts. No finite experiment is used to extrapolate the arbitrary-ambient theorem: its proof is(3), the four complete TT4 polynomial identities, positivity(6)--(8), and restoration(9).

## 8. Dependency pins

The complete eight-term quartic expansion and path-resolvent derivative identities are the reviewed Session189 odd-transitive packet, Sections1--3; its Section7 and the accompanying independent audits pin the inherited S188 all-cover restoration. The positive ambient transfer is also the elementary prefix/suffix/excursion decomposition used in the separately frozen skew-polynomial packet. The general source-family scalar/profile and highest-node reductions are inherited S188 results explicitly cited in those reviewed packets.

The newer incoming common-end proof is at commit97632252db49c2d3dce55cb9393a1cb0725c263b, path05-knowledge/results/reversal_cube_mod4_anchor_20261104/PROOF.md, Git blob51390c8fe6423e7bd5a7326495ae81e01121c3aa. It is mentioned only for the conditional complementary integration above. Its preregistration or claimed status is not treated as independent approval.

The independent session189-incoming-anchor-audit/REVIEW.md approves the incoming common-end proof with the explicit quotient-order correction m>=2, and records its exact source and physical controls. This packet assumes at least one other core child, so its own setting satisfies that condition. The broader audited common-end lemma also allows arbitrary nontransitive distinguished interiors sharing a physical source or sink; that strengthening is not needed for the TT4 integration.
