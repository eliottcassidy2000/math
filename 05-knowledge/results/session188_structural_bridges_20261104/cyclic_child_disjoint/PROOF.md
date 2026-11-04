# Disjoint quartics of a cyclic child and pair-free alternating-hub separation

## 1. Fixed physical observer and exterior path algebra

Let T have a module A={p,q,r} oriented p->q->r->p, and let W converse only A. Put R=T-A. For distinct physical i,j in R, consider

    Omega_ij = [a_p a_i c_q c_j](F_T-F_W).

Ports are the actual x->y, with a_v=[v->x], c_v=[v->y]. All four old names are distinct. This is the exact disjoint quartic of the complete response, not a cubic approximation.

Work over the commutative squarefree algebra on variables z_v for v in R, so z_v^2=0. Let G_uv be the sum of all directed paths in R from u to v, weighted by the product of their vertex variables. A singleton path is included. Define

    E_i=sum_u G_ui,       S_i=sum_v G_iv,
    I={v in R:v->A},      O={v in R:A->v},
    ell_i=sum_(v in I) G_iv,
    rho_i=sum_(v in O) G_vi,
    alpha=1+sum_(v in I) E_v,
    beta =1+sum_(v in O) S_v,
    tau=sum_(u in O,v in I) G_uv.

Every product discards intersecting vertex supports. Write [R] for the coefficient of the product of all z_v. In particular these products are exact disjoint path-cover sums, not independent endpoint marginals.

**Disjoint-kernel theorem.**

    Omega_ij = [R] {
        2 alpha beta (G_ij-G_ji)
        -(1+2 tau)(E_i S_j-E_j S_i)
        +alpha(rho_i S_j-rho_j S_i)
        +beta(E_i ell_j-E_j ell_i) }.                   (1)

This holds for arbitrary R and arbitrary homogeneous incidences I,O.

## 2. Rank-one elimination proof

Give p,q,r their own squarefree variables. An excursion from one A-name through R to another has generating series tau independent of its two A-endpoints, by modularity. Thus eliminate R and use the effective A adjacency matrix B_A+tau J, where J is all ones. Its path matrix M is the finite squarefree resolvent

    M=(I-Z_A(B_A+tau J))^(-1) Z_A.

Repeated names vanish automatically. Products of exterior excursions likewise retain disjoint support. Let C_u and R_u denote the column and row sums of M, and H=sum_(u,v in A) M_uv. The full old path matrix has

    G^T_up=ell_u C_p,    G^T_pu=R_p rho_u,
    G^T_uv=G_uv+ell_u H rho_v,
    E^T_p=alpha C_p,     S^T_p=beta R_p,
    E^T_u=E_u+alpha H rho_u,
    S^T_u=S_u+beta ell_u H.

Here u,v are exterior names; G^T_pq=M_pq. These identities follow by splitting a path at its successive visits to A. Conversely, concatenation recovers exactly one such path. Empty exterior prefix/suffix are the constants in alpha,beta.

Converse of A transposes M, exchanges C_u and R_u, and fixes H. Let z_A=z_p z_q z_r. Directly from the three possible A names,

    M_pq-M_qp=z_p z_q(1-z_r(1+2 tau)),
    [z_A](C_p C_q-R_p R_q)=0,
    [z_A](C_p R_q-R_p C_q)=2,
    [z_A]H(M_pq-M_qp)=1.                               (2)

The differences of the C/R products begin in A-degree three, so multiplication by H kills their z_A coefficient. Also H^2(M_pq-M_qp) has degree at least four and cannot contribute. For example C_p=z_p+tau z_p z_q+(1+tau)z_p z_r+z_A(1+2tau+2tau^2), while R_p interchanges the two degree-two coefficients; the other rows/columns follow cyclically. These formulas verify every identity in (2).

The complete disjoint quartic normal form is the sum of the following eight disjoint three-cover products:

    E_p G_iq S_j + E_q G_jp S_i
    +E_p G_ij S_q + E_j G_qp S_i
    +E_i G_pq S_j + E_q G_ji S_p
    +E_i G_pj S_q + E_j G_qi S_p.

Substitute the eliminated expressions, subtract their M-transposes, and extract z_A. The first and last paired terms vanish by the second identity in (2); the two cross products give 2 alpha beta(G_ij-G_ji). The M_pq terms give the remaining three terms of (1). This proves the full integer coefficient identity.

A related general identity requires no three-name assumption: for any module A, distinct p,q in A and distinct i,j outside, the local A-converse contrast satisfies

    [a_p a_q c_i c_j](F_T-F_(A conversed))=0.            (3)

The same elimination works with arbitrary A. Its eight-term expression factors through C_p R_q+C_q R_p, which is invariant under transpose, while H and every exterior factor are fixed. Exchanging ports gives the corresponding c_p c_q identity. This is an exact same-port support restriction, not a statement that all disjoint quartics vanish.

## 3. Applying the kernel to quotient converse

Let T=Q[A,D_1,...,D_k], with the distinguished child A a directed triangle, and let U=Q^op[A,D_1,...,D_k], all interiors and names unchanged. Suppose i and j lie in two different exterior children, and each is endpoint-balanced in its own child: its ordered-cover endpoint counts E_i(h) and S_i(h), in one prescribed cover slot, agree for every h, and likewise for j.

Whole converse preserves the entire quartic part of the response. After whole-conversing U, the quotient is Q and every child interior is conversed. In each eight-term three-cover expression, an exterior child containing exactly one physical mark contributes only one marked endpoint profile E_i or S_i. Its converse exchanges these, so the balanced hypothesis permits restoring it. Every unmarked child contributes its ordinary q profile, which is always converse-invariant. Thus

    [a_p a_i c_q c_j](F_T-F_U)=Omega_ij.                 (4)

Only the C3 child remains conversed. No joint same-child marked profile is discarded: there are exactly two marks in A and one in each of the two other marked children. Vertex-transitive children satisfy the balanced condition, with E_i(h)=S_i(h)=q_h/|D_i|.

Consequently (1) is an exact 2+1+1 disjoint-quartic kernel for arbitrary quotients with these marked children. Formula (3) gives zero when both marks in A occur at the same port under the same exterior balance hypotheses.

## 4. Alternating-hub positivity without any pair child

Let k=2r>=2. Take the alternating-hub quotient W_(2r+1): chain indices 0->1->...->k-1, with a hub beating even indices and losing to odd indices. Inflate the hub by A=C3 and each chain index t by an arbitrary nonempty tournament D_t. Choose i in D_0 and j in D_(k-1). Assume

    E_i(h)=S_i(h)=b_0(h),  E_j(h)=S_j(h)=b_(k-1)(h),
    b_0(1)>0, b_(k-1)(1)>0.                            (5)

For example each extreme child may be any vertex-transitive tournament, including a singleton or a directed triangle. Other chain children have no endpoint-balance requirement.

**Theorem.** The coefficient (4) is strictly positive. Therefore one of its sixteen actual queries separates the substitution from quotient converse with interiors fixed, for every r and every allowed child size. In particular all-C3 substitutions of every W quotient form an unbounded pair-free family satisfying this theorem.

### Positive polynomial reduction

Compress exterior runs to chain colors z_0,...,z_(k-1). A chain path has

    G_uv=z_u z_v product_(u<t<v)(1+z_t), u<v,
    G_uu=z_u,             G_uv=0, u>v,
    E_v=z_v product_(t<v)(1+z_t),
    S_u=z_u product_(t>u)(1+z_t).

Use these in (1), with I the odd indices and O the even indices. The final linear functional L sends z_t^h to q_h(D_t) at unmarked chain children, and to b_t(h) at the two marked extremes. It factors across colors. This is exact run-cover compression; powers count distinct runs, not repeated physical names. All weights are nonnegative and their order-zero values are zero.

Put g=G_(0,k-1), e=E_(k-1), s=S_0. Since the extremes are respectively outgoing/incoming neighbors of the hub,

    rho_0=z_0,  ell_(k-1)=z_(k-1),
    S_u g=G_(u,k-1)s for every even u,
    E_v g=e G_(0,v) for every odd v.

Thus the first term of (1) and the negative pieces of its last two terms cancel exactly to (alpha+beta)g. The entire color expression becomes

    (alpha+beta)g + (1+2tau) e s
        +(alpha+beta-1-2tau) z_0 z_(k-1).              (6)

The first two summands have nonnegative coefficients. The first has strictly positive L-value, at least

    2 b_0(1)b_(k-1)(1) product_(0<t<k-1) H(D_t).       (7)

The last summand also has nonnegative L-value. If k>=4, covering all intermediate colors forces the alpha term to end at k-1, the beta term to start at 0, and the tau term to start at 0 and end at k-1. Direct multiplication gives its L-value as

    [b_0(1)b_(k-1)(2)+b_0(2)b_(k-1)(1)]
          product_(0<t<k-1) H(D_t) >=0.                (8)

For k=2, add b_0(1)b_1(1) to this expression. This accounts for the otherwise absent constant term. Equations (6)–(8) prove strict positivity without a bounded inference or a sign assumption about a general exterior kernel.

The primary structural proof that W is strong and prime for r>=2 is retained in ../alternating_prime_quotient/PROOF.md. At r=1 it is C3; the same positivity proof and the explicit k=2 boundary apply.

## 5. Root-containing families when every child is C3

For W_(2r+1)[C3,...,C3], every proper nontrivial original module is exactly one of the triangle children, by primeness of the quotient and of C3. Hence every root-containing original-module XOR family consists of the root together with some of these children.

If the hub child is selected, its internal orientation is restored. Reversals of other children cannot affect the quartic in (4), because its marked extreme children are endpoint-balanced and all remaining children are unmarked. The positive coefficient therefore separates.

If the hub child is not selected, its three names are exposed anchors: queries uniform on every other child make the target exactly the whole-converse response. The compressed uniform-anchor theorem in ../returns_next/UNIFORM.md gives a residue-two separator. This branch uses actual queries and does not discard quartics of an arbitrary partial contrast.

Thus every root-containing family is visible, including the fully covered family selecting the root and all children. Families omitting the root are disjoint child reversals and are separated by module isolation. The reasoning applies from every quotient/child-cycle-orientation state: reversing the chain labels identifies W^op with the same alternating-hub convention, and every triangle's two orientations remain endpoint-balanced. Therefore the 2^(2r+2) tournaments obtained by independently choosing W or W^op and the orientations of all 2r+1 triangle children have distinct full named responses.

This is a structurally proved cube distinction, not an exhaustive response-bank claim. A fixed detection battery against one source root-containing family can combine the sixteen quartic queries with the compressed anchor battery at one hub name. There are m+2 compressed blocks when m=2r+1; the inherited bound gives at most m+20 queries in their union. This is not claimed as a common reconstruction battery for all cube states.

## 6. Scope

The all-order conclusions are the exact disjoint kernel, the same-port zero identity, alternating-hub positivity, and the all-C3 root-family/cube consequence. General prime quotients outside this family remain open here. The selected hub child must be C3 in (1), not an arbitrary larger child. Arbitrary strong outer ambient lifting is NOT asserted: these positive integer quartics need not have a nonzero residue in a transfer theorem's modulus. The whole displayed core and separately justified ordered SCC padding are the scope.

Dependencies are the complete quartic normal form in ../response/PROOF.md, exact module run-cover substitution, and the separately audited uniform-anchor and alternating-hub results. No new historical priority is claimed.
