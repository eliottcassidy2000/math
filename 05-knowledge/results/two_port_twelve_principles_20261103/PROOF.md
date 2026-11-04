# Ambient residue transfer and the cubic converse observer

**PROVED: T1-T4 below.** Controls are **FINITE-EXACT**, separately scoped.
Universal named two-port reconstruction, exact even endpoint imbalance, and
exclusion of every twelve-vertex module remain **OPEN**. No full-response
collision or literature-priority claim is made.

All tournaments have physical, fixed vertex names. A bit a_u=1 means u->x;
c_u=1 means u->y; the exterior arc is x->y. F_T(a,c) is the Hamilton-path
count of this actual extension. Empty induced sets have H=1. A module has
uniform incidences to every exterior vertex. Its substitution at z in Q is
Q[M]. The binary relation is the actual arc; there are no ties or name gauges.

## Inheritance and dependencies

The closest mechanisms are B3-B4 in
[the ambient descent proof](../two_port_seven_core_20261103/algebra/PROOF.md),
P1-P3 in [squarefree parity](../two_port_beyond_nine_20261103/audit/PARITY.md),
and A4 in [the full three-path formula](../two_port_beyond_nine_20261103/algebra/PROOF.md).
The hostile is the eight-vertex endpoint-parity-balanced nonselfconverse
core (4,17,42,131,45,75,31,119). The corrected near miss is primitive
response normalization: ordered triangle factors can have arbitrarily
large odd content. The least-used sidecar here is the parity of the size
of each path's support, summed over all supports rather than only the full
one. The [board](BOARD.md) records the portfolio and losses.

T2-T3 use the inherited odd-cycle congruence

    H(T) = 1 + 2 o(T) (mod 4),

where o counts directed simple odd cycles, with rotations identified, on
all subsets. This is THM-466,
[the 2-adic digits of H](../../../01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md),
from THM-002 [OCF](../../../01-canon/theorems/THM-002-ocf.md). We import the
proved congruence, not its older speculative interpretation of digit data.
Our checker also compares literal cycle enumeration with independent path DP.

## T1. Full response modulo four descends through every ambient quotient

**Theorem T1.** For any two nonempty same-size tournaments M,N and any
quotient Q with distinguished z,

    H(Q[M])-H(Q[N]) = H(Q)(H(M)-H(N)) (mod 4).       (1)

Consequently, if Q[M] and Q[N] have identical full m-port responses, then
M and N have identical full m-port responses modulo four. This holds for
every m>=0 and each fixed tournament orientation among the ports.

**Proof.** The uniform-run expansion, retained with integer skeleton counts,
is H(Q[M])=sum_k c_k q_k(M), where q_k counts ordered k-path covers and
c_1=H(Q). For k>=2, use unordered set partitions Pi of the common named set:

    q_k(M)-q_k(N)
      = k! sum_Pi (product_(B in Pi) H(M[B])
                    - product_(B in Pi) H(N[B])).

Each parenthesized difference is even, because every nonempty induced
tournament has odd H. Thus the difference is divisible by 2k!, in particular
by four. Only k=1 survives modulo four, proving (1).

For any requested attachments inside M and N, give every new port the same
incidences to the old exterior vertices as z. The old module together with
the ports is then a module in the actual ambient extension. Apply (1) to
these two larger modules. H(Q) is odd and invertible modulo four. Equality
of the ambient query therefore implies equality of that core query modulo
four. All attachment masks and all exterior port tournaments are allowed.
QED.

This is a necessary descent test, not an exact cancellation theorem. It
does not divide the full integer response by H(Q), discard path-cover runs,
or require response content one. It compares arbitrary same-size module
pairs; converse is only one application. Equal modules are positive controls.

## T2. The full modular converse contrast has an explicit cubic formula

Work in R=F2[z_v]/(z_v^2). Let G enumerate simple directed paths including
singletons, with their complete squarefree vertex support, and put

    s=G1, e=G^t1, U=product_v(1+z_v),
    p=a^t s, q=c^t s, X=a^t G c.

For f in R, L_even(f) and L_odd(f) sum its coefficients on even and odd
supports respectively. These are linear functionals, not ring evaluations.
In particular replacing z by 1 inside products would be invalid.

Define the Boolean polynomial

    delta_T(a,c) = (F_T(a,c)-F_(T^op)(a,c))/2 mod 2.

The division is over the integers, before reduction: both counts are odd.
Its unique Boolean coefficients are well-defined by integer finite differences.

**Theorem T2.** Exactly

    delta_T = L_even((e+s)^t(a+c))
       + L_odd((e+s)^t a + pq + (U+1)(p+q)X).       (2)

In particular it has degree at most three. For u!=v the coefficient of
a_u a_v c_v (also of a_v c_u c_v) is

    h_uv(T) = sum_(W subset V-v, |W| even)
                        (S_u(T[W])+E_u(T[W])) mod 2. (3)

Terms not containing u are zero. At even n>=2 this simplifies to

    h_uv(T) = sum_(W subset V-v, |W| even, |W|<=n-4)
                                      S_u(T[W]) mod 2. (4)

For n=12 the support sizes are precisely 2,4,6,8; the empty term is zero.
Every hidden ambient module converse must have all h_uv=0 by T1.
Their logical independence from the old combined sieve is not asserted.

**Proof of (2).** By the odd-cycle congruence, delta is the parity difference
of odd cycles in the two extended tournaments. Cycles using no port cancel
under path reversal. A cycle using only x closes an even-support old path;
its attachment monomial is a_end(1-a_start). Reversing that path changes
this, in characteristic two, by a_end+a_start. The two one-port terms
therefore give L_even((e+s)^t(a+c)).

An odd cycle using x->y has an odd-support old path from after y to before x.
Reversing it gives the difference (e+s)^t a+a^t(G+G^t)c. Inherited P1 is
G+G^t=ss^t, so the latter term is pq.

Finally an odd cycle through both ports, with neither adjacent to the other,
is uniquely a path P after x and a disjoint nonempty path Q after y. Its
contribution is L_odd(((1+a)^t G c)((1+c)^t G a)). No factorial or independent
endpoint-marginal product is introduced. Under converse replace G by G^t.
Inherited P2 gives e=Us and U^2=1; P1 gives c^tGa=X+pq. Every positive-degree
element of R squares to zero, hence p^2=q^2=0. Expanding the two products,

    (Uq+X)(Up+X+pq) + (q+X+pq)(p+X)
       = (U+1)(p+q)X.

This retains the paired paths and proves (2), including singleton paths.

**Proof of (3).** Only the last, cubic term can contribute. For a_u a_v c_v,
disjoint support forces the middle path in X to be the singleton v and
the endpoint in p to be u. The other apparent choice repeats u or v in
two paths and vanishes. Its coefficient is
L_odd((U+1)z_v s_u). Removing v changes odd to even support and removes
v from every old path. P2 on T-v turns (U_(V-v)+1)s_u into e_u+s_u, proving
(3). The other monomial follows by the symmetry in a,c of the cubic term.

**Proof of (4).** Put D=V-v, so |D|=n-1 is odd. For a monomial z_W in s_u,
the number of even supersets in D is even if |D-W|>=2, is one if |D-W|=1,
and is zero if W=D. Thus L_even(U_D s_u) is the sum of the size n-2
coefficients of s_u. In L_even((U_D+1)s_u), those terms cancel, leaving
exactly the even support sizes at most n-4. This also handles n=2,4,
where the sum is empty. QED.

The full three-path tensor has not been replaced by H,S,E,C. Equation (2)
is a specific converse projection modulo four, with a separate paired-path
term. Zero of this polynomial is necessary for a hidden module, not sufficient
for an integer full-response collision.

## T3. A parity-only path interpretation of the old cover-skew obstruction

**Theorem T3.** For distinct u,v,

    (C_uv(T)-C_vu(T))/2
       = number of simple directed u-to-v or v-to-u paths
         having an odd number of vertices                 (mod 2). (5)

**Proof.** C_uv-C_vu is even by inherited P3. Flipping the mutual edge changes
H by C_reverse-C_forward, since the Hamilton polynomial is affine in the
two mutually exclusive arc variables; each gradient entry is independent
of that mutual edge. Modulo four, the odd-cycle congruence makes this change
twice the parity of the odd cycles destroyed or created. Removing the flipped
edge is a bijection from those cycles to directed paths between u,v with an
even number of edges, hence an odd number of vertices. Neither path uses the
mutual edge. Signs disappear after reduction modulo two. QED.

Thus C=C^t modulo four can be tested using only path parities and odd support
size. This is an equivalent implementation of B4's intrinsic condition, not
an additional theorem excluding its residual universe. The even-support
part is missing from this observer; P3 uses all supports for C modulo two.

## T4. Exact converse comparison needs only a cubic query battery

**Theorem T4.** Over the integers, F_T-F_(T^op) has total degree at most three
and bidegree at most (2,2). Pure a and pure c quadratic coefficients vanish.
Consequently the following at most 1+2n+n^3 attachment assignments suffice
to decide exact equality of these two *core* responses: zero; each singleton
variable; each pair a_u,c_v; and each triple a_u,a_v,c_w or a_w,c_u,c_v with
u!=v. This is not a battery for comparing arbitrary tournament pairs.

**Proof.** Reversing the entire extended Hamilton word and exchanging the
two port names gives

    F_(T^op)(a,c)=F_T(1-c,1-a).

Each response has degree at most two in either port's bits. Its degree-four
part only counts words in which both ports have two old neighbours. Such
words never use the port-port arc. Exchanging the port names therefore makes
this homogeneous part symmetric in a,c, even while the fixed exterior arc
is x->y. The degree-four parts in the displayed identity cancel exactly.

With a=0, x is a universal source. Every Hamilton path starts at x, and
the remaining count is the one-port response at y. Its quadratic coefficients
are -(C_uv+C_vu), unchanged by converse. For c=0, y is *not* a universal
source because x->y. Instead specialize the inherited integer squarefree
formula A4, writing P=1^tG1:

    F_T(a,0)=[V]{1+P+2(e^t a)+2(e^t a)P-(e^t a)(a^t s)}.

Its quadratic coefficients are again -(C_uv+C_vu). Thus both pure-quadratic
parts of the difference vanish. The remaining possible monomials are
the constant, 2n singletons, n^2 mixed pairs, and 2n binom(n,2) mixed triples.
Boolean finite differences recover their coefficients from the stated
assignments; omitted pure-pair values are the sums of their singleton values.
The count is 1+2n+n^3. QED.

At twelve vertices this is 1753 assignments (one Hamilton count for each
of the two converse cores per assignment), versus 6241 for the complete
weight-at-most-two-per-port interpolation battery. This saves comparison
queries; it supplies no smaller complete universe of tournaments to examine.

## Validation and scope boundary

Run `python check.py --output check.json` and the corresponding `python -O`
command. The standard-library checker pins its source, uses explicit errors,
compares complete query tables through order four, literal Hamilton words,
literal odd cycles, integer endpoint/cover DP, and independent squarefree
path multiplication. Universes and query counts are in the receipts.
The full small truth tables also test the exact degree-three bound and
vanishing pure-quadratic difference, independently of residue reduction.

The inherited residual remains even, nonselfconverse, size at least twelve.
Equations (2)-(5) add exact representations and tests, not an assertion that
the residual is empty. The incoming paired-endpoint forest grammar is
compatible with these path products; its arbitrary-depth state cannot be
replaced by this one modular, two-port projection.
