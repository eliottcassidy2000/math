# A transitive-triple circulation equals the positive cyclic-child quartic

## 1. Physical statement

Let A={0,1,2} be a module of a named tournament T, with physical order 0->1->2 and 0->2. Let W converse only A. The actual ports are x->y, with a_v=[v->x] and c_v=[v->y]. For distinct exterior physical names i,j, put

    Omega_pq = [a_p a_i c_q c_j](F_T-F_W).

Let V replace only the internal order of A by the cycle 0->2->1->0. Write V* for V with that triangle conversed, and define

    Omega_cycle = [a_0 a_i c_2 c_j](F_V-F_V*).

All exterior names, arcs and incidences are unchanged. Then the exact integer identity is

    Omega_01 + Omega_12 - Omega_02 = Omega_cycle.       (1)

This is a relation between complete disjoint quartic coefficients. The virtual cycle is a proof device, not a permitted relabeling or a replacement of the actual source. No assertion about lower-degree coefficients is hidden in (1).

## 2. Rank-one path proof

Work over the commutative squarefree algebra on exterior labels R=T-A. Use the notation G,E,S,I,O,ell,rho,alpha,beta,tau of ../cyclic_child_disjoint/PROOF.md. Explicitly G_uv sums all directed u-to-v paths, E_i sums paths ending at i, S_i sums paths starting at i, I consists of exterior names beating A, and O of names beaten by A. Define

    ell_i=sum_(v in I) G_iv, rho_i=sum_(v in O) G_vi,
    alpha=1+sum_(v in I) E_v, beta=1+sum_(v in O) S_v,
    tau=sum_(u in O,v in I) G_uv.

Products discard intersecting physical supports. Eliminate exterior excursions: the effective A path matrix is

    M=(I-Z_A(B_A+tau J))^(-1) Z_A.

Let C_p and R_p be its column and row sums, and H its total sum. Exterior elimination is exact because each successive A visit has an independently marked, disjoint exterior excursion. Converse of A transposes M, interchanges C and R, and fixes H.

For the fixed exterior marks define

    K_S=alpha(ell_i S_j+ell_j S_i)
         -beta(E_i rho_j+E_j rho_i),
    K_D=alpha beta(G_ij-G_ji),
    B_0=E_i S_j-E_j S_i,
    B_1=alpha(rho_i S_j-rho_j S_i)
         +beta(E_i ell_j-E_j ell_i).

Write S_pq=C_p C_q-R_p R_q, D_pq=C_p R_q-C_q R_p and E_pq=M_pq-M_qp. The complete eight-term disjoint three-cover expression, reproduced in the cyclic-child packet, gives after subtraction and extraction of z_A:

    Omega_pq=[R][z_A]{S_pq K_S+D_pq K_D
                          +E_pq B_0+H E_pq B_1}.     (2)

To check truncation, S_pq and D_pq have A-degree at least three, while E_pq has degree at least two. Before extraction the omitted terms have the form H S_pq, H D_pq or H^2 E_pq, hence A-degree at least four and vanish. Thus (2) is exact for three old names, not an asymptotic truncation of the actual observer.

Direct path enumeration on the three effective names gives the following full-support coefficients, as polynomials in tau:

    order / (p,q)     [z_A]S  [z_A]D  [z_A]E       [z_A]HE
    TT3 / (0,1)         -2       0       0              1
    TT3 / (1,2)          2       0       0              1
    TT3 / (0,2)          0      -2       1+2tau         1
    cycle / (0,2)        0       2      -1-2tau         1

For example M_02-M_20=z_0 z_2+(1+2tau)z_A in the transitive case, whereas the cycle gives z_0 z_2-(1+2tau)z_A. The table yields

    Omega_01=[R](-2K_S+B_1),
    Omega_12=[R]( 2K_S+B_1),
    Omega_02=[R](-2K_D+(1+2tau)B_0+B_1),
    Omega_cycle=[R](2K_D-(1+2tau)B_0+B_1),

and proves (1). No sign assumption on K_S or B_1 is needed. In fact B_1 can have either sign in finite arbitrary-exterior controls.

## 3. Grouped positivity and a 24-query battery

Suppose R contains two nonempty disjoint sets X,Y such that X union Y is a module of R, X->Y, and A->X while Y->A. The remaining exterior is arbitrary. Identify a_i=a_X for every i in X and c_j=c_Y for every j in Y, keeping all other exterior attachment bits fixed arbitrarily. In every Omega use the grouped quartic [a_p a_X c_q c_Y]. This equals the sum over physical i in X,j in Y: bidegree at most (2,2) leaves no higher-degree alias contribution once each port already has its distinct A and exterior mark.

The aggregated cyclic-child positivity theorem, ../aggregated_cyclic_module/PROOF.md together with ../ordered_cyclic_module/PROOF.md, applies to the virtual cycle. Therefore

    Omega_01+Omega_12-Omega_02
      >= 2 H(X) H(Y) H(E+A_singleton)>0,              (3)

where E=R-(X union Y), and the singleton has A's incidences toward E. H denotes the Hamilton-path count. This is all-order strict positivity of the signed circulation. Individual TT3 quartics need not be positive.

The three fourth differences combine without loss into 24 actual queries. With all other A bits zero, use the following six (a restricted to A,c restricted to A) masks and weights:

    (1,2): +1; (2,4): +1; (1,4): -1;
    (2,0): -1; (0,2): -1; (0,0): +1.

For each, toggle a_X and c_Y independently and multiply by the two-variable finite-difference sign (-1)^(u+v). These 24 masks are distinct because A,X,Y are disjoint and nonempty. Their signed source-minus-target sum is exactly (3), so at least one separates.

## 4. Quotient converse and every fixed ambient

Let T=Q[A,D_1,...,D_k], with A=TT3 in the physical order above. Let U=Q^op[A,B_1,...,B_k] on the same child sets, where every B_t has the same complete ordered nonempty path-cover profile as D_t. Assume X,Y are two child sets with the cyclic ordered-module witness from Section 3 after deleting A.

For each of the three saturated grouped quartics, the general master identity of ../quartic_cover_transfer/PROOF.md replaces quotient converse by local A converse. Its hypotheses require no triangle: uniform attachment to each non-A child permits exact q-profile restoration in every ordered-cover coordinate. Hence the same 24-query circulation separates T and U.

Embed the whole core in any fixed rooted ambient tournament, retaining every physical exterior name. Give both actual ports the core's exterior incidences at every query. The augmented core remains a module, so the common outer run-cover substitution and the all-cover top-quartic reversal lemma still yield that master identity. Apply Section 3 to the entire exterior: the witness X union Y remains a module after deleting A because every ambient name sees the core uniformly. Thus (3), with E enlarged by all ambient names, remains strictly positive. This covers arbitrary strongly connected outer contexts, without transferring the sign of an arbitrary core coefficient through a congruence argument.

The designated root restriction is real. The theorem applies, for example, at a cyclic ordered-module witness in the critical T/U/W quotients, with TT3 at that witness root and arbitrary nonempty remaining children. It does not assert that every root of every prime quotient has such a witness.

## 5. Scope

The all-order conclusions are the exact physical circulation, its grouped strictness, the 24-query quotient-converse separator with equal non-A q profiles, and its arbitrary fixed-ambient lift. This packet does not claim separation for every original-module family inside the distinguished TT3. Its proper interval modules can permute the physical transitive order, which is not covered by the restoration assumption that the target A has the original interior. It also does not claim general prime-quotient closure or unrestricted named reconstruction.

Finite controls retain exact physical query values and a complete named polynomial bank of bidegree at most (2,2), including every cubic and quartic coefficient and every zero. They complement rather than replace the formal path proof. The path elimination and top-degree response identities are inherited from the cited frozen packets; the three-edge transitive circulation is the specialization derived here.
