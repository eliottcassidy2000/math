# Deletion currents in the exact named cubic

**PROVED D1-D3; FINITE-EXACT D4 and the declared controls.**
Universal prime twelve-core exclusion and exact full-response collision
remain **OPEN**.

## Inheritance and definitions

Incoming **C2/C4**,
[two_port_eight_core_20261103/algebra/PROOF.md](../two_port_eight_core_20261103/algebra/PROOF.md),
exclude homogeneous pairs and leave a strong, prime, nonselfconverse
twelve-core. Incoming **T1/T4**,
[two_port_twelve_principles_20261103/PROOF.md](../two_port_twelve_principles_20261103/PROOF.md),
give arbitrary-ambient response descent modulo4 and the full integer cubic
for whole-core converse. The nearest mechanism for the following identity
is **R1**, [deletion contraction](../two_port_unit_distance_20261103/reconstruction/PROOF.md),
together with Session27's exact [one-port current, section6](../session27_repunits_lattices_20261022/session27-response/proof.md).
D1 is a coefficient-level refinement of that inherited deletion observer.

All old names remain fixed. F_T(a,c) counts Hamilton paths after adding
x->y with a_i=[i->x],c_i=[i->y]. Let Delta_T=F_T-F_(T^op).
S,E are start/end Hamilton counts; C_uv counts ordered spanning covers
by two nonempty directed paths, first ending at u and second starting at v.
Define the exact one-port current

    g_u(T)=2(E_u-S_u)+sum_w(C_uw-C_wu).                  (1)

For the one-port function f_T(b), with b_i=[i->z], Session27 proves

    f_T(b)-f_(T^op)(b)=sum_u g_u(T)b_u.                 (2)

This is a complete criterion for whole-converse one-port equality. It is
not a criterion for full two-port equality or for arbitrary ambient module
invisibility.

## D1. A shared named cubic is exactly a deletion current

For distinct u,v,

    [a_u a_v c_v] Delta_T
       =[a_v c_u c_v] Delta_T
       =-g_u(T-v).                                    (3)

**Proof of the first identity.** In a contributing Hamilton word, the two
port incidences at v force either x,v,y or y,v,x. To supply a_u also,
the first order extends to u,x,v,y. Its signed coefficient is minus the
number of spanning covers of T-v whose first path ends at u and whose
second path may be empty: -E_u(T-v)-sum_w C_uw(T-v).

For the second order, supplying a_u requires y,v,x,u. A preceding old
vertex at y would introduce a different c-bit, so y must be first.
This contributes +S_u(T-v). These are the only possible words for the
three named variables. Thus the coefficient in F_T itself is
S_u-E_u-sum_w C_uw, computed on T-v.

Converse exchanges S,E and transposes C, giving (3) by subtraction.
For the second monomial family, the whole-converse identity
Delta_T(a,c)=-Delta_T(1-c,1-a), together with T4's cubic degree bound,
makes its homogeneous cubic part symmetric in a,c. Alternatively reverse
the same two word families. The old name v is shared, but a_v,c_v are
distinct Boolean variables. Singleton middle paths are essential. QED.

Equation(3) evaluates only one family of cubic coordinates. The remaining
mixed cubics and all lower terms remain in the retained full polynomial.

### Exact versus ambient scope

Exact equality F_T=F_(T^op) forces g(T-v)=0 for every v.
For a hidden module M in an arbitrary quotient, T1 instead gives
Delta_M=0 modulo4. Its necessary consequence is

    g_u(M-v)=0 modulo4 for every u!=v.                 (4)

This weaker congruence is the one usable directly in the unknown ambient
problem. Do not substitute exact zero for (4). No quartic cancellation is
asserted for a partial module reversal in an asymmetric ambient tournament.

The symmetric positive-cut weights from
[Session176 circuits](../session176_tournament_information_20261104/circuits/PROOF.md)
give an additional interpretation. If
w_uv=(C_uv+C_vu)/2, then

    g_u=2(sum_(v->u) w_uv - sum_(u->v) w_uv).           (5)

Thus the exact shared cubic observes weighted divergence on the deletion
card. These weights are realized Hamilton-cover counts; arbitrary formal
balanced flows are not interchangeable with them.

## D2. An endpoint-balanced deletion card excludes its parent core

Let T have at least two vertices. Put
b_u(T)=S_u(T)+E_u(T) modulo2. For every v,

    b_v(T)=1+sum_(v->u) b_u(T-v) modulo2.              (6)

**Proof.** A Hamilton word beginning at v next visits some out-neighbor u,
so S_v(T)=sum_(v->u) S_u(T-v).
Likewise E_v(T)=sum_(u->v) E_u(T-v).
Add these and use sum_u E_u(T-v)=H(T-v)=1 modulo2. QED.

Therefore b(T-v)=0 implies b_v(T)=1. A hidden module core has b(M)=0
by inherited **B3**, [ambient descent](../two_port_seven_core_20261103/algebra/PROOF.md);
consequently **none of its vertex-deletion cards can be
endpoint-parity balanced**. This is an all-order necessary obstruction,
not a classification of the remaining cards.

## D3. Every odd vertex-transitive card is excluded, with arbitrary birth edges

Let Q be vertex-transitive of odd order m, with h=H(Q). Transitivity gives
S_u(Q)=E_u(Q)=h/m; this integer is odd. Add a new vertex v with exactly k
out-neighbors in Q. Then

    S_v=k*h/m,    E_v=(m-k)*h/m,
    S_v-E_v=(2k-m)*h/m,                                (7)

which is nonzero and odd. Therefore every such extension fails the intrinsic
endpoint-parity requirement, irrespective of its other properties. It cannot
be a hidden module core in any ambient order.

This explains the complete rejection of all births over cyclic11 and
Paley11 in earlier scans, and covers all other odd vertex-transitive cards
without scanning their attachment masks. It does not cover every eleven-card
of a prime twelve-core. Ordinary regularity is insufficient: the inherited
regular seven-vertex H171 example has unequal endpoint data and a visible
one-port current.

## D4. Complete finite exclusion in one prime twelve-core switching family

Take the eleven-name tournament whose residues j-i modulo11 are
{1,3,4,5,9}, and adjoin source11. Its outgoing rows are

    (570,1140,233,466,932,1864,1681,1315,583,1166,285,2047).

For each X subset {0,...,10}, reverse exactly the edges crossing X,X^c.
These are2^11=2048 distinct named tournaments: fixing the source's gauge
removes the global sign duplication. The complete integer subset-DP
certificate gives S_u+E_u odd at **every** name of **every** member.

Exactly2024 members are strong and prime. The remaining24 have a universal
source or sink. Primality was checked by all4082 proper nontrivial subsets
in the author algorithm and by an independent pair-generated module-closure
algorithm; strongness is checked by transitive closure.

| H | Named presentations | Strong and prime |
|---:|---:|---:|
|95095|24|0|
|196815|264|264|
|236447|1320|1320|
|256023|440|440|

By the inherited endpoint-parity obstruction, none supports a hidden
module converse in any ambient quotient. This is a finite exact theorem
about the stated entire orbit; it does not exclude arbitrary prime12 cores.

## Retained full response and hostile controls

The three prime switches at masks1,3,11 have full whole-converse cubic
contrasts with respectively1582,1574,1578 nonzero cubic coefficients.
The separately inherited prime core

    (3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002)

has H300453, balanced endpoint parity, and1584 nonzero cubic coefficients.
It preserves the branch which the conference family never enters; it
fails other modular conditions and is not a collision.

The complete bidegree(2,2) responses for all four cores and their converses
are saved, with physical masks and integer coefficients, in
[response_tables.json](response_tables.json). At twelve this is6241
coefficients for each full response. The contrast is exactly cubic and
uses the1753 allowed supports of T4. No sorted profile, quotient by labels,
or scalar normalization replaces it. The corresponding actual mod4 queries
are in [response.json](response.json).

D1 is independently checked against1056 entries from those full tables.
Both scalar/path and spectral hostiles are retained in [BRIDGE.md](BRIDGE.md).

## Other finite observations and the stopping reason

The authenticated complete bank of198947 isomorphism classes in orders2..9
contains respectively0,1,0,1,0,2,0,3 classes with exact g=0.
All seven are strong and regular, but this finite fact does not prove
regularity or odd order for every exact one-port-blind core.
None has g(T-v)=0 on every deletion card.

The SAT scout encoded necessary endpoint parity and shared-cubic parity
under a valid Hamilton-path labelling gauge. It imposed neither the full
ordered-cover condition nor all full-response coefficients. The runs at8
and12 returned UNKNOWN after their requested conflict budgets; the12-run
reported40377 conflicts and no model. No UNSAT certificate, exact collision,
or prime12 residual candidate resulted. Its records are exploratory stopping
evidence only. Extending that search without a stronger reduction is not a
proof of the remaining branch.

