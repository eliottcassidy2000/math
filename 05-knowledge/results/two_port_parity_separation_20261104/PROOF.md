# The cover/current parity certificate

**PROVED E1-E4; same-session proof review and independent exact algorithms.**
Local theorem names E1-E4 below; cite this file with the identifier.
This file does not claim unrestricted two-port reconstruction of arbitrary
different tournaments. The comparison is a fixed named core with its
whole converse, and subsequently a single module reversal in a fixed quotient.

## Inherited definitions and mechanisms

All vertices have fixed physical names. Add x->y with a_u=[u->x] and
c_u=[u->y]; F_T(a,c) is the resulting Hamilton-path count. Set

    Delta_T=F_T-F_(T^op),      delta_T=Delta_T/2 modulo2.

Division is over integers before reduction, in values and Boolean
coefficients. Every tournament has odd H, so it is well-defined.
All responses are multilinear Boolean polynomials. Whole converse has the
full cubic bound of **T4**, [twelve principles](../two_port_twelve_principles_20261103/PROOF.md).
No quartic cancellation is asserted for an arbitrary internal reversal.

S_u,E_u count Hamilton paths starting/ending at u. C_uv counts ordered
spanning two-path covers ending at u then starting at v. Define

    b_u=(S_u+E_u) modulo2,
    D_uv=(C_uv-C_vu)/2 modulo2,      D_uu=0,
    h_uv=[a_u a_v c_v]delta_T,       u!=v.

The half-cover skew is integral by **P3**,
[squarefree parity](../two_port_beyond_nine_20261103/audit/PARITY.md).
**T2/T3**, [twelve principles](../two_port_twelve_principles_20261103/PROOF.md),
give

    D_uv=[a_u c_v]delta_T
        = parity of odd-support simple paths between u,v either way, (1)

and, at even n,

    h_uv=sum_(A subset V-v, |A| even, 2<=|A|<=n-4) S_u(T[A]) mod2. (2)

Singleton paths are included in the underlying path matrix. The cubic
part creates no quadratic terms by Boolean reduction: a repeated port
index in the same port factor repeats an old name in two disjoint paths
and is zero. Thus (1) is the actual mixed quadratic coefficient of T2.

**D1**, [exact deletion currents](../prime_twelve_current_20261104/PROOF.md),
also identifies

    h_uv=-g_u(T-v)/2 modulo2,
    g_u=2(E_u-S_u)+sum_w(C_uw-C_wu).                         (3)

The sign disappears in characteristic two. S187's
[complete support cutoff](../session187_named_response_20261104/prime12/PROOF.md)
retains every mixed cubic and supplies a second derivation of (2).
The increment below is a compatibility relation between (1) and (2),
not a replacement of the other named cubic coefficients.

## E1. At every vertex, the two parity obstructions cannot both vanish

For every tournament of order n>=3 and every v,

    sum_(u!=v) D_uv + sum_(v->u) h_uv = 1 modulo2.           (4)

**First proof, using the short-support cutoff at even n>=4.** Work in the squarefree algebra over F2, with path matrix G,
start vector s=G1, end vector e=G^t1, and U=product_w(1+z_w).
Inherited P2 says s+e=(U+1)s. By (1),

    sum_(u!=v) D_uv = L_odd(s_v+e_v).

The singleton term is counted twice on the right and cancels, consistently
with the zero diagonal of G+G^t. For a monomial z_A in s_v, its weight in
L_odd((U+1)s_v) is the parity of the number of odd supersets of A,
plus the indicator that |A| is odd. If at least two names are unused,
the first number is even. If exactly one name is unused, even n makes
|A| odd and the two contributions cancel. With no unused name both
contributions vanish. Hence

    sum_(u!=v) D_uv
       =sum_(A contains v, |A| odd, |A|<=n-3) S_v(T[A]) mod2. (5)

Since n>=4, the singleton {v} occurs in (5), with contribution1.
Every other path counted there starts v,u with v->u. Remove its first
vertex v. This is a bijection onto paths starting at u, avoiding v, with
even support size2,...,n-4. By (2), their total parity is sum_(v->u)h_uv.
Thus (5) is1+sum_(v->u)h_uv, proving (4). QED.

This proof uses a linear coefficient functional; it never substitutes z=1
inside a squarefree product. The singleton survives the cutoff precisely
from n=4. Formula(4) is not asserted at n=2.

**Second proof, by path partition, at every n>=3.** The general shared-label
identity T2(3), equivalently (3) above, says

    h_uv=sum_(A subset V-v, |A| even) (S_u(T[A])+E_u(T[A])) mod2.

Absent endpoints and the empty set contribute zero. In sum_(v->u)h_uv,
prepend v to the paths starting at u. These are exactly the odd-support
paths of size at least three starting at v. For the paths ending at u,
the condition v->u is complementary to u->v. The latter condition
appends v and counts exactly the odd-support paths ending at v. The sum
of both endpoint conditions, over each fixed even nonempty A in V-v,
counts H(T[A]), which is odd. Thus

    sum_(v->u)h_uv
       =sum_(u!=v)D_uv
          + number of nonempty even subsets of V-v (mod2)
       =sum_(u!=v)D_uv + (2^(n-2)-1) (mod2).

For n>=3 the last term is1. This proves (4) without the squarefree
superset cancellation, and shows why its constant is unavoidable. For
n=2 it is0 instead, explaining the exact failure boundary. QED.

In current notation, (4) reads

    sum_(u!=v)(C_uv-C_vu)/2
       +sum_(v->u)g_u(T-v)/2 = 1 modulo2.                   (6)

Consequently cover-skew symmetry modulo4 and vanishing of all shared-label
cubics modulo4 are incompatible at every n>=3. Endpoint balance and
primality are not needed for this incompatibility.

## E2. Twelve actual queries exclude every endpoint-balanced twelve-core

The general even-order statement is stronger. Fix v in a tournament of
even order n>=4; let d be its outdegree. Define n literal queries, all
with c=e_v:

    a=e_u              for each u->v;
    a=e_u+e_v          for each v->u;
    a=e_v if d is odd, and a=0 if d is even (one final query).

Call their set Q_v. The n attachment masks are distinct. Then

    sum_((a,c) in Q_v) delta_T(a,c) = 1+b_v modulo2.        (7)

In particular, if b_v=0, at least one of these n named queries has
Delta_T(a,c)=2 modulo4. For endpoint-balanced twelve-cores this holds
at every one of the twelve possible anchors, without any additional
assumption about primeness, strongness or the other cubic coordinates.

**Proof.** T2 and P2 give the singleton-a identity

    delta_T(a,0)=sum_u b_u a_u.                            (8)

Indeed its coefficient is L_all(s_u+e_u). In (U+1)s_u every proper
support has odd enlargement weight and the full support has weight zero;
P2 identifies this with S_u(T)+E_u(T). No pure-a quadratic remains.
Also delta_T(0,0)=0 and sum_u b_u=0, since sum S=sum E=H.

Put f(a)=delta_T(a,e_v). Boolean finite differences give

    D_uv=f(e_u)+f(0)+b_u,
    h_uv=f(e_u+e_v)+f(e_u)+f(e_v)+f(0).                    (9)

The second equation includes all lower-degree terms: those with c=0
cancel by (8), and each f value contributes exactly the indicated
two-a finite difference. Substitute (9) into (4). The f(e_u) terms
cancel for out-neighbors and remain for in-neighbors. The coefficient
of f(e_v) is d modulo2, and that of f(0) is n-1+d=1+d modulo2.
The remaining b sum is sum_(u!=v)b_u=b_v. This is (7). QED.

Equivalently the actual integer contrasts satisfy

    sum_((a,c) in Q_v) Delta_T(a,c) = 2(1+b_v) modulo4.    (10)

If b_v=1, the single query(e_v,0) already separates by (8). Thus n+1
specified queries at any chosen anchor certify whole-converse separation
for every even n>=4. For a known endpoint-balanced core only the n
queries in Q_v are needed. This is a detector for the stated converse
comparison, not a reconstruction theorem from these n values.

## E3. Every nontrivial tournament differs from its named converse modulo4

**Theorem.** For every named tournament T on n>=2 vertices, some named
two-port query has F_T-F_(T^op)=2 modulo4.

Even n>=4 is E1 or E2. At n=2, name the actual arc0->1. The query
(a,c)=(e_0,0) gives F_T=3 and F_(T^op)=5, so their difference is2 modulo4.

For odd n>=3, sum(1) over unordered endpoint pairs. Every odd-support
directed path on at least three vertices occurs exactly once. Oddness of H
on each induced subset gives

    sum_(u<v)D_uv=sum_(k odd, k>=3)binom(n,k)
                  =2^(n-1)-n=1 modulo2.

Therefore a mixed quadratic coefficient of delta_T is nonzero. The
integer Boolean finite-difference transform is invertible modulo4, so
some actual query on its two involved attachment bits separates.
This odd-order mechanism is inherited from the earlier odd-cover theorem;
the even-order incompatibility is the new step. QED.

At orders0 and1 converse changes no arc, so no such separator exists.
Whole-converse isomorphism does not evade E3: the vertex names and query
masks are fixed, with no simultaneous relabeling of attachments.

## E4. No single nontrivial module converse is hidden in any ambient quotient

Let Q have distinguished root z, let M have at least two fixed names, and
form T=Q[M], T'=Q[M^op] with exactly the same exterior names and arcs.
By E3 choose a core query with Delta_M=2 modulo4. Give each new port the
root z's incidences to all old exterior names, keeping that core query
inside M. The extended core remains a module in each actual extension.
Inherited **T1**, [ambient residue transfer](../two_port_twelve_principles_20261103/PROOF.md),
then gives exactly

    F_T(lifted a,lifted c)-F_T'(lifted a,lifted c)
       = H(Q) Delta_M(a,c) = 2 modulo4,

because H(Q) is odd. Thus the ambient full named responses differ.

This closes in particular the endpoint-balanced strong prime twelve-core
branch in every ambient order. It does not cancel a strong quotient over
the integers, claim the ambient internal-reversal contrast is cubic, or
exclude collisions between arbitrary unrelated tournaments or simultaneous
changes of several modules. Those broader injectivity questions require
their own reduction and are not consequences asserted here.
