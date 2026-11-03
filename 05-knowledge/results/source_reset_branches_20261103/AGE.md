# A complete module-age divisor law and a linear native reset

**PROVED all sizes in the operations below.** The formulas for the four
rooted coordinates are inherited from [BLOCK_SWITCH B1](../adaptive_state_induction_20261103/BLOCK_SWITCH.md).
The continuation-divisor interpretation is inherited from
[ROOTED_RESETS R1/R3](../two_port_return_induction_20261103/ROOTED_RESETS.md).
The new deductions are the age classification, stable rotation, clock
selector, and their precise observation/implementation boundaries.

## A1. Objects and the oldest/non-oldest classification

Let(T,r) be a nonempty rooted tournament. Use

    z=(h,s,q,F), a=h-s, d=F-q,

where h=H(T), s counts paths starting at r, q=H(T-r), and F counts spanning
two-slot path covers whose nonempty first path ends at r; the second may
be empty. Empty tournaments have H=1. Add k>=2 successive vertices at r,
each with r->x->all other old vertices. The new vertices form a transitive
module. A marked new vertex w has A newer vertices before it and B older
vertices after it in that module's directed order; A+B+1=k.

The inherited profile of w is

    h'=a+q+(2^k-1)d,
    s'=2^B d                         (A>0),
        a+(2^B-1)d                   (A=0),
    q'=a+q+(2^(k-1)-1)d,
    F'=2^A a+(2^(B+1)3^A-2^A)d+2^(A+1)q
        + 1_(B=0)(s-q).                                      (1)

Forming a root-child interface at w gives exact all-future two-port
divisor D_w=gcd(h',s',q',F'). This divides every continuation count in
the inherited six-command language, and is their gcd, not just a bound.

**Theorem A1.** Put g0=gcd(a,q,d) and g=gcd(d,a+q). Then

    D_w = g0                         if B>0,
    D_w = gcd(g,s+(2^(k-1)-1)q)       if B=0.                 (2)

Thus every non-oldest position has the same divisor. The exceptional
oldest position can decrease, preserve, or increase it; there is no
universal monotonicity claim across a reset.

**Proof.** Every tournament H count is odd, including q and q'. Since
h'-q'=2^(k-1)d, any common divisor of h',q' divides d (the divisor is odd),
and then divides a+q. If A,B>0, s'=2^B d and F' reduces modulo g to
2^A q. Hence the gcd is gcd(g,q)=g0. For A=0, q'-s'=q and
F'=h'+q, again giving exactly gcd(q,a,d). At B=0, s'=d, so the first
three coordinates have gcd g. Reducing F' modulo g, with a=-q and d=0,
gives s+(2^(k-1)-1)q. Conversely each displayed divisor divides all four
coordinates. This proves equality in all cases. The integer g is positive
and odd: a+q>0, and h' odd prevents a+q,d from both being even. QED.

At fixed z,k a single oldest/non-oldest bit therefore suffices for this
divisor observation. When the two values differ, zero bits cannot suffice.
This is not a one-bit carrier for later native moves or full H responses.

**Theorem A2.** If T is strong of order at least3, the k marked positions
have pairwise different single-bridge future-H observers on their converses.

Indeed d>=a>0 by the inherited positive rooted defect and strongness.
For A>0, s'=2^B d distinguishes the k-1 positions. For the newest position,
s'=a+(2^(k-1)-1)d exceeds 2^(k-2)d, the largest of the others (also at k=2).
On the converse graph the golden coordinates are

    v=(h',(F'-s')/2,(F'+s')/2).                            (3)

They retain s', and the inherited three-probe theorem recovers v from
future H counts. Consequently any exact observer for all such contexts
must distinguish at least k states among these actual presentations.
This is a target-relative lower bound, not a theorem that a fixed number
of integer coordinates cannot encode the ages.
The halves in(3) are integers: the exact adaptive birth identity
H(B_w(U))=H(U)+F_w-s_w and oddness of both tournament counts give
F_w-s_w even, without restricting U to a previously generated word bank.

## A3. Stable rotation improves the implementation clock

Keep the newest physical name z fixed as the port. Move it from the first
to the last position of its transitive module, leaving every other leaf
in its previous relative order. Reverse just the k-1 arcs between z and
the other leaves. All exterior arcs and all other old arcs are unchanged.
The endpoint is isomorphic to the original unmarked tournament, so H is
unchanged, and z is now the oldest mark. This is a specified graph operation,
not an assignment of a new numerical state.

The inherited [single-arc compiler, METAGRAPH §5](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md)
implements every reversed arc using one temporary vertex, two actual
Hamilton-path complements, and deletion of that temporary name. Therefore
the exact clocks for this stable rotation are

    temporary births = temporary deletions = k-1,
    native complements = 2(k-1),   peak order = |T|+k+1.       (4)

This replaces the whole-module reversal in
[NATIVE_RESET_INDUCTION N7](../two_port_return_induction_20261103/NATIVE_RESET_INDUCTION.md),
whose sufficient clock was k(k-1) native complements. Both produce an
oldest physical mark; they differ on the order of the unmarked leaves.
The older divisor-breaking family is inherited, not claimed anew.

For that family, original profile(5,2,3,6), fixed original root3, and k>=2,
the old newest-mark interface has divisor3; the new oldest-mark interface
has divisor1. Its next oldest-root birth still has the inherited count
6*2^k+2*3^k-1. The stable rotation has an explicit marked reentry: create
the fixed name z as the first leaf, keep ports(3,z), then create the other
leaves in their original chronological order, always discarding the new
name. Every endpoint is in the actual two-port language. No intermediate
two-port membership, strongness, or H neutrality is asserted.

**Theorem A3 (source-retaining lower bound).** Any implementation by L
native complements which retains all N=|T|+k original vertices throughout,
while allowing arbitrary temporary births/deletions, satisfies

    L >= (k-1)/2             if L is even,
    L >= (N-1)/2             if L is odd.                    (5)

This assumes births retain old-old arcs and deletions remove only temporary
vertices; it does not apply to replacing or deleting original vertices.
On the persistent name set the final change graph is

    (L mod2) K_N XOR E_1 XOR ... XOR E_L,

where E_i is the set of pairs of persistent names adjacent in the ith
Hamilton spine. Every E_i has maximum degree at most2, even if a temporary
vertex occurs between two persistent names. The target change graph is a
star of degree k-1 at z, with degree zero at every exterior vertex. Compare
degrees at z for even L and at an exterior vertex for odd L. This proves(5).
Thus the new clock is optimal in order of growth in k within this declared
source-retaining language. Exact optimal native lengths remain open.

## A4. Complete arithmetic selector for the oldest mark

For fixed actual z, define

    D(k)=gcd(g,s+(2^(k-1)-1)q),  g=gcd(F-q,h-s+q).          (6)

This function is periodic in k with a period dividing ord_g(2); for g=1
use period1. This gives a terminating exact classifier for every k>=2:
enumerate the doubling orbit modulo g and retain its gcd values. Its least
period may be smaller than ord_g(2). No search over tournament orders is
needed once the actual rooted source is supplied.

For each prime power p^e dividing g, the full precision guard is exactly

    p^e | D(k)  iff  q*2^(k-1) = q-s (mod p^e).            (7)

If q is a unit this is a membership test for one point in the doubling
orbit; otherwise keep the displayed linear congruence rather than divide
illegally. The [older order-tower theorem](../fermat_wieferich_order_20261014e.md)
gives ord_(p^e)(2)=d_p p^max(0,e-s_p), where d_p=ord_p(2) and
s_p=v_p(2^d_p-1). Thus Wieferich delay changes the selector clock when
such a prime actually divides this source's g. No claim that1093 or3511
occurs in our finite source census is made.

Two adjacent test ages have the exact common divisor

    gcd(D(2),D(3))=gcd(h,s,q,F).                           (8)

Proof: subtract (s+3q)-(s+q)=2q and use g odd, then expand a,d. This does
not imply that either D(2) or D(3) equals that gcd: different primes can
block the two ages. Nor does a primitive original profile by itself prove
that some D(k)=1; (6) gives the exact finite orbit question for that source.

The actual rooted source with rows

    (74,36,17,6,11,29,190,63), root2,
    adaptive roots(0,2,1,0,6), profile(179,10,41,566)

has a=169,d=525,g=105. As k=2,...,13 its D values are

    3,7,3,5,21,1,3,35,3,1,21,5.

This is an exact period12; the first primitive age is k=7. The three local
forbidden conditions on t=k-1 are t odd modulo2 (prime3), t=0 modulo4
(prime5), and t=2 modulo3 (prime7). Their union explains the initial delay;
the first good exponent is t=6. This is joint incidence of clocks, not an
independent-prime probability argument.

The complete H<=295 adaptive census contains25198 rooted presentations,
22856 with primitive original profile. Every one has a primitive oldest
reset, with first k<=7; seven fail at both k=2 and3. These are finite-exact
facts in that source universe, not a proof for all rooted tournaments.

## A5. What the arithmetic quotient forgets

The 3^A term disappears from(2), but remains in F' and therefore in every
complete golden observer. For the(5,2,3,6) family, k=2 and4 both give D=1
yet give different F49 coordinates. At k=6 the oldest profile is
(195,3,99,1649), with v=(195,823,826), field coordinate0 and charge6.
It is primitive. A nonzero-phase-only field encoding would lose this case.

Modulo7, the full family(1) can be computed from k and A using actual powers
of2 and3 and the boundary indicators; it is not the field Frobenius map.
If a subsequent q2 attachment is used, that separate operation multiplies
alpha by theta, theta^2=4theta+2, and retains charge. Frobenius is still
alpha->alpha^7 with theta^7=4-theta. The exact source profile, positional
boundary and operation label prevent the three clocks from being conflated.
