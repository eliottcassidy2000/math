# Binomial mixed-quartic tomography for every odd transitive child

**New exact all-order identities and a scoped root-family consequence; independent review pending.** A binomial combination of named three-inside/one-outside quartics eliminates every internal excursion term for a transitive module of any odd order. It recovers exactly twice a contracted exterior two-cover skew. Combined with internal-order tomography and the inherited anchor and ordered-witness theorems, this extends the earlier transitive-triple root-family result to odd transitive children of every order at least five. No arbitrary-prime closure, unrestricted reconstruction, or full-response collision is claimed.

## 1. Physical setup

Let A be a module of an arbitrary named tournament T. Write R=T-A. Fix three distinct names p,q,r in A and one name i in R. Let T' reverse only A, fixing every exterior name and arc. The ports have x->y and a_v=[v->x], c_v=[v->y]. Define

    Delta_(pq;r,i)=[a_p a_q c_r c_i](F_T-F_T').

This is a complete named quartic coefficient. Its four Boolean variables are distinct even when other formulas use a physical name in both ports. Here p,q,r,i themselves are all distinct.

Use the squarefree algebra on all physical old names, with repeated support zero. On R let G be the simple-directed-path matrix, including singleton paths, S_i its row sum, and E_i its column sum. Define

    alpha=1+sum_(v->A)E_v,   beta=1+sum_(A->u)S_u,
    tau=sum_(A->u,v->A)G_uv,
    ell_i=sum_(v->A)G_iv,    rho_i=sum_(A->u)G_ui.

Every support extraction below is joint in this algebra. The exterior may be arbitrary and strongly connected; no factorization into independent endpoint marginals is made.

## 2. A general two-channel local-converse identity

On A introduce a formal commuting variable t and the simple-path matrix M(t) whose edge weight from u to v is [u->v in T[A]]+t. Its entries also contain the squarefree variables of every visited A name. Let

    C_j=sum_u M_uj,    R_j=sum_v M_jv,    H=sum_j C_j.

(The letter R_j in this display is an internal row sum, not the exterior vertex set.) Put

    X_(pq;r)(t)=[A](C_p M_qr+C_q M_pr-M_rp R_q-M_rq R_p).

**Lemma 1.** For every module A and every four specified names,

    Delta_(pq;r,i)=[R]{(alpha S_i-beta E_i) X_(pq;r)(tau)
                  +alpha beta(ell_i-rho_i) X'_(pq;r)(tau)}.       (1)

The derivative is the ordinary formal derivative in t.

**Proof.** The complete eight-term disjoint-pair quartic normal form is the sum of ordered three-path-cover products E_u G_vw S_z over both orientations and orders of the two marked pairs. Eliminate the exterior R by its maximal runs. The full path matrix restricted to A is M(tau). Full endpoints in A have end sum alpha C_j and start sum beta R_j. The cross entries are G_(q,i)=R_q rho_i and G_(i,q)=ell_i C_q; the full exterior start and end sums are S_i+beta H ell_i and E_i+alpha H rho_i.

Consequently the quartic before taking its converse difference is

    alpha(C_p M_qr+C_q M_pr)(S_i+beta H ell_i)
      +alpha beta(C_p R_q+C_q R_p)R_r rho_i
      +alpha beta C_r(C_p R_q+C_q R_p)ell_i
      +beta(E_i+alpha H rho_i)(M_rp R_q+M_rq R_p).

All products retain disjoint support. Under converse of A, M is transposed, C and R are exchanged, and H is unchanged. Subtraction gives the first exterior channel times X, plus alpha beta(ell_i-rho_i) times

    [A]{H(C_p M_qr+C_q M_pr-M_rp R_q-M_rq R_p)
                 +(C_p R_q+C_q R_p)(C_r-R_r)}.

Differentiating the path resolvent gives M'=M J M, where J is the all-ones matrix. Hence M'_uv=R_u C_v, C'_j=H C_j, and R'_j=H R_j. The last display is exactly X'(t). Substitution t=tau proves (1). No endpoint marginal is substituted for a joint cover. QED.

## 3. Exact transitive kernel

Now let A have physical transitive order 0->1->...->m-1, with m>=3. Fix p=0, q=1, and r in {2,...,m-1}. Put L=m-2 and u=r-1. Define the ordinary polynomial

    P_u(s)=(s+1)^u s^(L-u)-s^u(s+1)^(L-u).

For an ordinary polynomial P define its finite Eulerian transform

    E(P)(t)=sum_(j>=0) (Delta^j P)(0) t^j,                 (2)

where Delta P(s)=P(s+1)-P(s). The sum is finite. Then

    X_(01;r)(t)=2 E(P_(r-1))(t).                          (3)

In particular X has degree at most m-3. This is a symbolic identity, not an interpolation from the controls.

**Proof.** In the transitive simple-path matrix G_0 on A,

    (G_0)_ij=z_i z_j product_(i<h<j)(1+z_h)   for i<j,
    (G_0)_ii=z_i,
    C_j^0=z_j product_(h<j)(1+z_h),
    R_j^0=z_j product_(h>j)(1+z_h),
    H_0=product_h(1+z_h)-1.

Rank-one expansion gives

    M=G_0+t R^0 C^0/(1-tH_0),
    C=C^0/(1-tH_0),   R=R^0/(1-tH_0).

Let B be the m-3 unmarked names and put P_left=product_(2<=h<r)(1+z_h), P_right=product_(h>r)(1+z_h). Their cardinalities are a=r-2 and b=m-1-r. After extracting the three already present marked variables, the definition of X gives

    X/2=[B]{P_left/D
            +t P_left P_right(P_left-P_right)/D^2},
    D=1-t(P_left P_right-1).                              (4)

Expand D^-1 and D^-2 geometrically in t/(1+t). Using [B]P_left^v P_right^w=v^a w^b and shifting the second sum yields

    X/2=(1/(1+t)) sum_(s>=0)(t/(1+t))^s
       [(s+1)^(a+1)s^b-s^(a+1)(s+1)^b].                 (5)

Newton's finite-difference expansion P(s)=sum_j Delta^jP(0) binom(s,j), together with the geometric series for binom(s,j), identifies the right side with (2). This proves (3). The leading powers in P_u cancel, so its degree and that of X are at most m-3. QED.

## 4. Odd order cancels every excursion term

Suppose m=2k+1 is odd. For each exterior name i define a physical linear combination of quartic contrasts

    B_i=Delta_(01;m-1,i)
           +sum_(r=2)^k (-1)^r binom(m-2,r-1) Delta_(01;r,i).    (6)

For m=3 the sum is empty. Let z be a virtual singleton with A's original exterior incidences in R+z, and let C_uv denote the ordered spanning two-path-cover capacity ending at u and starting at v, with the gap arc untested.

**Theorem 1.** For every odd m>=3 and every fixed exterior R,

    B_i=2 K_i,
    K_i=C_zi(R+z)-C_iz(R+z).                              (7)

For m=5 the combination is Delta_(01;4,i)+3 Delta_(01;2,i). For m=7 it is Delta_(01;6,i)+5 Delta_(01;2,i)-10 Delta_(01;3,i).

**Proof.** Since L=m-2 is odd, pairing opposite terms in the binomial identity ((s+1)-s)^L=1 gives

    P_L(s)+sum_(u=1)^((L-1)/2)
                   (-1)^(u+1) binom(L,u) P_u(s)=1.       (8)

Apply E and (3). The corresponding combination of X polynomials is exactly2, so the combination of their derivatives is exactly0. Lemma1 therefore makes (6) equal to2[R](alpha S_i-beta E_i).

A cover counted by C_zi has first path ending at z and second starting at i. Deleting z gives an alpha prefix and an exterior path starting at i, with joint full support; therefore C_zi=[R]alpha S_i. Similarly C_iz=[R]beta E_i. This proves (7). QED.

The identity cancels every term involving tau, including all possible exterior reentries, rather than assuming they vanish modulo a power of two. Its integer coefficients may be negative; each named quartic is recovered by a legitimate sixteen-query Boolean finite difference. There are k=(m-1)/2 quartic terms, so at most8(m-1) actual queries suffice per exterior mark.

## 5. The even-order boundary is explicit

The same calculation does not give a constant kernel for even m. For L=m-2 even, pairing the expansion of ((s+1)-s)^(L-1)((s+1)+s)=2s+1 gives a kernel combination2(1+2t), with derivative4. The coefficients are1 at r=m-1 and

    (-1)^r [binom(m-3,r-1)-binom(m-3,r-2)]

at r=2,...,(m-2)/2. Thus the resulting contrast is

    2[R]{(alpha S_i-beta E_i)(1+2tau)
                         +2 alpha beta(ell_i-rho_i)}.    (9)

For m=4 this is already the single r=3 contrast. No general nonvanishing of (9) is asserted. No extension of the odd theorem to all even transitive children is hidden in the statement.

## 6. Odd transitive children with an ordered exterior witness

Use the reviewed S188 exterior-skew theorem in [tt3_degree_tomography/PROOF.md, Section3](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/tt3_degree_tomography/PROOF.md), with its [final independent scope review](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/tt3_degree_tomography_audit/REVIEW.md). That skew theorem itself depends only on the virtual root z and its exterior, not on the size of the distinguished child.

Specifically, suppose R has an ordered module X->Y with A->X and Y->A. Suppose X,Y have the same complete ordered-cover profile q, and q is log-concave. The cited theorem proves

    K_X-K_Y >= 2 q_1^2 H(E+z_M)>0,
    K_X=sum_(i in X)K_i,  K_Y=sum_(i in Y)K_i,
    E=R-(X union Y).                                    (10)

Group the c_i bit uniformly over X or over Y in each quartic in (6). Because the monomial already uses two a variables and two c variables, no higher alias route exists; its coefficient is exactly the sum over the physical names of that group. Therefore the difference of the two grouped binomial combinations is

    2(K_X-K_Y) >= 4 q_1^2 H(E+z_M)>0.                    (11)

This is valid for every odd transitive A and every fixed ambient included in R.

## 7. Root-containing original-module family extension

Let Q be prime with k_Q children, and let its child z be an odd transitive tournament A=TT_m with m>=5. Suppose distinct children u,v satisfy z->u->v->z and {u,v} is a module of Q-z. Let their physical child tournaments X,Y have identical log-concave complete ordered-cover profiles. Every other child is arbitrary nonempty. Put the core into any fixed ambient. Only operations inside the core are allowed in this root-family statement.

**Theorem 2.** Every simultaneous XOR family of original source modules selecting the whole core oddly is detected from the source by a fixed finite two-port battery in every fixed ambient.

**Proof.** Include the core-internal zero query, with both ports imitating the core at all exterior names. If the scalar core Hamilton count strictly increases, hereditary original-module profile monotonicity and positive ambient substitution detect it there. Otherwise all corresponding child profiles agree and the target quotient is Q^op.

Include all internal shared-name quartic queries from the independently reviewed Session189 transitive-distance theorem. If they do not separate, its all-cover restoration and reconstruction show that the target physical A is either A or A^op.

If it is A^op, the general conversed-child anchor lemma, explicitly inherited from S188 UNIFORM U2--U4, separates. Its battery leaves A's names as singletons and compresses every other child.

If it is A, use the two grouped mixed binomial batteries of Section6. Top-degree invariance under whole old-core converse and equal-profile restoration of all other uniform children identify the quotient-converse quartic comparison with the local A-converse comparison, in every augmented-cover coordinate. This is exactly the already reviewed all-cover quartic restoration mechanism. Both ports imitate the core on the exterior, so ambient substitution preserves that identity. The source ordered-witness hypotheses still hold after adjoining the ambient names, since they see the whole core uniformly. Equation(11) gives a strictly positive actual quartic combination, hence one actual query in its finite-difference battery separates. QED.

One conservative fixed query bound, without removing duplicates, is

    1 + 16*m*binom(m-1,2) + 16*(m-1) + (m+k_Q+1).

These are respectively the scalar query, internal-order tomography, the two grouped binomial contrasts, and the compressed anchor battery. The entire battery depends on the source and the specified ambient incidences, not on which original modules are selected. It is source-relative detection, not all-pairs reconstruction among arbitrary outputs.

The m=3 theorem was already proved separately in S188 by degree tomography. The new argument extends its ordered-witness source-family conclusion to every odd m>=5. It does not remove equal-profile or log-concavity hypotheses on X,Y. With common odd transitive children, the inherited critical T/U/W quotient witnesses supply concrete unbounded prime-quotient families; the present packet does not reprove their classification.

## 8. Reproduction and boundaries

The accompanying verifier retains every inner kernel coefficient, checks the formal derivative channel, compares literal ordered two- and three-cover enumeration with (3), and verifies the odd/even binomial combinations. The physical controls include complete named source/local-converse response vectors, all cubic/quartic zero coordinates, complete physical masks, and exact contracted skews. Finite examples support the implementation; the all-order assertions follow from Sections2--7.

The fixed-ambient mixed identity is unconditional. Its nonzero root-family consequence needs the separately stated ordered-witness hypotheses. No prime-eight general closure or full named-response collision is claimed.
