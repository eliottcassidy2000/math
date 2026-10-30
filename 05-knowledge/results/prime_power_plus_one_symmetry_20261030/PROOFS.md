# Prime-power-plus-one growth: three different symmetry laws

2026-10-30. **PROVED scoped deductions; CITED classical inputs;
FINITE-EXACT controls; same-agent review.** No literature-priority claim.
All groups preserve literal tournament arcs unless explicitly called
converse maps. Counts of specified families are not A000568 totals.

**Concurrent credit.** Incoming commit51cd593a3e already proves the
general carousel family and scalar/regular-orbit recurrences, including
the unique ternary fixed palette and the1554/87460002 counts; see
[T1--T2, GROWTH.md](../prime_power_plus_one_20261029/GROWTH.md).
G1 recovers that shared result and retains its full order polynomial.
G3 extends the preceding bounded-observable limits to nonuniform
palettes and every fixed positive moment. G5 uses a root-of-unity
carrier distinct from the incoming norm-one conic and paired lattice.

## G1. An exact child inventory at every p^k+1

For odd prime p let R_p be the carousel: vertices Z/p, with x->y when
1<=y-x<=(p-1)/2 modulo p. Put T_(p,0)=one vertex and
T_(p,k+1)=R_p[T_(p,k),...,T_(p,k)]. Thus |T_(p,k)|=p^k.

R_p has full group C_p. Indeed translations act; an automorphism fixing
0 fixes its out-neighborhood and in-neighborhood, each a transitive
tournament with a unique linear order, hence fixes everything. It is a
prime tournament: its canonical maximal-module partition, if nontrivial,
would be a nontrivial block system of this transitive prime-degree action.
The canonical decomposition therefore gives the FULL iterated wreath group

    W_(p,k)=Aut(T_(p,k)),  |W_(p,k)|=p^g_k,
    g_k=(p^k-1)/(p-1).

The canonical-module mechanism is inherited from Session110's
[transitive-action carrier](../session110_overnight_20261030/CHECKPOINT_02.md).
For a palette of Q named colors define P_(p,k)(Q,z) as the sum over color
orbits of z^(log_p stabilizer order). Then

    P_(p,0)=Q,
    P_(p,k+1)=[P_(p,k)(z)^p-P_(p,k)(z^p)]/p
               +z P_(p,k)(z^p).                         (1)

Nonconstant p-tuples of child COLOR-ORBIT TYPES have a free root C_p
action. Constant tuples have stabilizer H wr C_p. This proves (1),
extending the credited ternary [S1 recurrence](../stabilizer_necklace_kernels_20261028/STABILIZERS.md).
It is the classical wreath-product distinguishing-color mechanism; see
[Chan](https://arxiv.org/abs/math/0601414). The polynomial records group
orders, not full abstract group types or future joint-port responses.

With ONE added vertex use Q=2. The regular-parent theorem
[R1--R2](../regular_parent_symmetry_20261029/PROOFS.md) applies: the new
vertex is intrinsically fixed, and the only duplicate root choices are
vertex clones. Since the parent is transitive, the exact UNMARKED child
inventory at order p^k+1 is

    E_(p,k)(z)=P_(p,k)(2,z)-z^(g_k-k).                   (2)

The subtractand is the point stabilizer, of order p^(g_k-k). For k=0
the same formula gives one two-vertex child. If F_k=P_(p,k)(2,1), then

    F_0=2,  F_(k+1)=(F_k^p+(p-1)F_k)/p,
    child classes=F_k-1.                                (3)

This removes child graph-isomorphism work within the stated regular-parent
sector once its subset-orbit inventory is known. It is not a benchmarked
speedup of a full tournament enumerator.

### Why the ternary branch is exceptional

Let A_k be the constant coefficient of P_(p,k)(2,z). Then

    A_0=2,  A_(k+1)=(A_k^p-A_k)/p.                      (4)

At p=3, A_k=2 forever; hence there are exactly two rigid unmarked children
for k>=2. For EVERY prime p>=5, A_1=(2^p-2)/p>2 and A_k strictly grows.
For k>=2 the rigid unmarked count equals A_k; for k=1 it equals A_1-1.
In particular A_(5,2)=1554 and A_(7,2)=87460002. The ternary fixed point
is an exact primitive-necklace recurrence fact, not a pattern inferred
from the first few exponents.

The labelled mass identity is an independent normalization:

    |W_(p,k)| P_(p,k)(Q,1/p)=Q^(p^k),
    |W_(p,k)| E_(p,k)(1/p)=2^(p^k)-p^k.                 (5)

Rigid classes proliferating does not mean a random binary attachment is
usually rigid. Its exact probability is p^g_k A_k/2^(p^k), and it is at
most (1-2^(1-p))^(p^(k-1)): every bottom block must be nonconstant.
This tends to zero for every fixed p.

## G2. Same size, different tree: flat and hierarchical growth

The flat carousel R_n, n=p^k, has full group C_n by the same fixed-vertex
argument, but T_(p,k) has the wreath group in G1. Both are regular.
For 0<=j<k the number of binary C_n-orbits with stabilizer p^j is

    [2^(p^(k-j))-2^(p^(k-j-1))]/p^(k-j).                (6)

There are two constant orbits with stabilizer p^k. This is the ordinary
least-period necklace count: a word of length p^a is nonprimitive iff
its period divides p^(a-1). A vertex stabilizer in C_n is trivial, so
subtract one from the rigid coefficient to obtain UNMARKED children.
For example, at order10 the flat parent R9 gives59 children,55 rigid;
the hierarchical parent C3[C3] gives23 children,2 rigid.
At order26 the flat R25 gives1342183 children,1342175 rigid; the
hierarchical R5[R5] gives6559 children,1554 rigid.

These are disjoint child families because their regular parents are
nonisomorphic (R2). At size7 there is already another base choice:
Paley7 has group21 rather than R7's group7. Its towers are not covered
by the p-group polynomial (1). Size p^k+1 does not determine the
automorphism action, its decomposition tree or the attachment law.

## G3. The previously missing weighted-tail bound

This addresses the explicit full-moment and nonuniform-palette targets
in [the regular-parent handoff, items3/5](../regular_parent_symmetry_20261029/HANDOFF.md).
Fix a faithful finite base group G on b points, |G|=h>1, and let
W_k=G wr ... wr G act on b^k leaves in its imprimitive action. Write
g_m=h^((m-1)/(b-1)) for the group order on a subtree of m leaves.
Let d=min_(g!=1)(b-number of cycles of g)>0.

Color leaves independently with an arbitrary finite probability palette
pi=(pi_a), and let a=max pi_a. Colors of zero probability may be removed.
For any fixed real s>=1 set M_(k,s)=E |Stab(coloring)|^s. Then

    M_(1,s)^(b^(k-1)) <= M_(k,s)
      <= M_(1,s)^(b^(k-1)) exp(epsilon_s),              (7)

    epsilon_s=(h^s-1)(h-1)
      sum_(j=2..k) b^(k-j) (g_(b^(j-1))^s a^(b^(j-1)))^d.  (8)

These inequalities hold at every finite k. They are useful when a is
small; there is no hidden assumption of a uniform or freely sampled
arithmetic palette.

**Proof.** For a child color orbit i, write H_i for its stabilizer order,
p_i its probability. Every literal coloring in this orbit has the same
probability w_i because leaf colors are IID. Thus p_i=g_m w_i/H_i.
Tilt the child-orbit measure by H_i^s, obtaining

    rho_i=p_i H_i^s/M_(j-1,s) <= g_m^s a^m.

Here w_i<=a^m, H_i<=g_m, s>=1 and M_(j-1,s)>=1. On a tuple of child
orbits the full stabilizer order is product H_i times the order L of
the base stabilizer of the orbit-type tuple. This follows by lifting
every permitted block permutation and multiplying the kernel sizes;
it does not assume a canonical section or identify joint port laws.
Consequently

    M_(j,s)/M_(j-1,s)^b = E_(rho IID) |L|^s.

A nonidentity base permutation fixes such a tuple with probability
product_over_cycles sum_i rho_i^length <= (max rho_i)^d.
The union bound over h-1 elements, and |L|^s-1<=h^s-1, bound the ratio
between1 and1+(h^s-1)(h-1)(g_m^s a^m)^d. Iterate and use log(1+x)<=x.
This proves (7)--(8). QED.

For fixed b,h,d,s and a->0, the sum is O(n a^(bd)), n=b^k: write
g_m^(sd)a^(md)=h^(-sd/(b-1)) [h^(sd/(b-1))a^d]^m and bound the sparse
sum by the geometric tail over all m>=b. Constants may depend on the
fixed action and s, never on k or the varying palette.

### Cyclic prime towers: the correct entropy is a collision moment

For G=C_p let c_p=sum_a pi_a^p. A bottom block has stabilizer C_p iff
all its colors agree, probability c_p. If n=p^k->infinity and

    n c_p -> lambda < infinity,                         (9)

then a^p<=c_p forces a->0, and n a^(p(p-1))<=n c_p^(p-1)->0.
Thus (8) tends to zero. The bottom moment is exactly
M_(1,s)=1+(p^s-1)c_p, so

    E |Stab|^s -> exp((p^s-1)lambda/p),                 (10)
    Stab converges in law to C_p^Z,
    Z distributed Poisson(lambda/p).                   (11)

The distribution assertion follows from the same higher-node collision
union bound without moment tilting: outside an event of probability o(1)
the full group is the product of the independent bottom stabilizers.
Equation (10) controls ALL fixed positive moments, not merely bounded
observables. For 0<s<1 it follows from convergence in distribution and
uniform integrability supplied by a larger moment. The convergence of
any fixed higher moment also proves uniform integrability of |Stab|^s.

For a uniform Q-palette, c_p=Q^(1-p), recovering the known critical
scale n/Q^(p-1). For nonuniform palettes, cardinality Q and squareclass
rank alone are insufficient. A color of mass1/2 keeps c_p>=2^-p even
if arbitrarily many tiny colors are added. The corresponding symmetry
does not disappear merely because the displayed palette grew.

### Paley-seven towers: a different collision statistic

For the actual Paley7 group of order21, the seven C3 subgroups each have
cycle partition(1,3,3), and distinct ones generate the whole group.
Put c_j=sum pi_a^j. The exact bottom law is

    P(H=C3)=7(c_3^2-c_7),  P(H=G)=c_7,
    M_(1,s)=1+7(3^s-1)c_3^2+(21^s-7*3^s+6)c_7.        (12)

If n=7^k and n c_3^2->lambda<infinity, then c_7<=c_3^(7/3) gives
n c_7->0. Also a^3<=c_3 and d=4 give n a^28->0, so (8) vanishes.
It follows that Stab converges to C3^Poisson(lambda), with

    E |Stab|^s -> exp((3^s-1)lambda)                    (13)

for every fixed s>0. Thus R7 towers are controlled by c_7, while Paley7
towers are controlled by c_3^2. This is a literal permutation-action
distinction at the SAME sizes, not a numerical analogy between3 and7.
All multiport statements fix the ports by name. They do not assert that
the enlarged unmarked graph has this entire automorphism group.

### All prime cyclic-complement affine bases

The same proof applies to the specified affine action
G=C_p semidirect C_l, where p,l are odd primes and l divides p-1.
Put t=(p-1)/l and d=t(l-1). Its p order-l subgroups have cycle type
(1,l,...,l); distinct ones generate G. Consequently, for any IID palette,

    P(H=C_l)=p(c_l^t-c_p),  P(H=G)=c_p,
    M_(1,s)=1+p(l^s-1)c_l^t+[(pl)^s-p l^s+p-1]c_p.   (13a)

If n=p^k and n c_l^t->lambda<infinity, then
c_p<=c_l^(p/l), so n c_p->0 because p/l>t.
Also a^l<=c_l and p d/l>t give n a^(p d)->0.
Equation(8) therefore yields, for EVERY fixed s>0,

    Stab -> C_l^Poisson(lambda),
    E |Stab|^s -> exp((l^s-1)lambda).                 (13b)

This concerns this FULL wreath action when used as a tournament law.
The incoming [local-base theorem D3](../prime_power_plus_one_20261029/GROWTH.md)
supplies literal full groups at P7, P11 and local13; its uniform-palette
distribution result is inherited, while(13a--b) add the nonuniform law
and moment control. Their critical statistics are c3^2, c5^2 and c3^4,
respectively. Equal palette size need not equal any of these statistics.

### Exact uniform-palette moment recurrences

With Q uniform, let a_m=g_m/Q^m. For cyclic prime towers,

    M_(k+1,s)=M_(k,s)^p
       +(p^s-1) a_m^(p-1) M_(k,p(s-1)+1), m=p^k.       (14)

For Paley7 towers,

    M_(k+1,s)=M_(k,s)^7
      +7(3^s-1)a_m^4 M_(k,s)M_(k,3s-2)^2
      +(21^s-7*3^s+6)a_m^6 M_(k,7s-6), m=7^k.          (15)

M_(0,s)=1. Sum over child color orbits, whose probabilities are a_m/H_i.
The constant-tuple and (1,3,3)-tuple contributions give these formulas.
Depth decreases in every right-hand moment, so they are finite recursive
evaluators even though moment order increases. A finite list of marginal
moments is not declared a closed joint-attachment carrier.

For the general affine base just specified, set m=p^k,
a_m=(pl)^((m-1)/(p-1))/Q^m. The same orbit calculation gives

    M_(k+1,s)=M_(k,s)^p
      +p(l^s-1)a_m^d M_(k,s) M_(k,l(s-1)+1)^t
      +[(pl)^s-p l^s+p-1]a_m^(p-1) M_(k,p(s-1)+1).   (15a)

This is an exact recursion, not an independence approximation. The last
coefficient corrects the shared constant tuple, which otherwise belongs
to every one of the p order-l fixed sets.

## G4. What the factors of p^k+1 actually prescribe

For any odd prime p and k>=1,

    v2(p^k+1)=v2(p+1) if k odd, and1 if k even.          (16)

For every odd prime ell dividing p^k+1, write a=v2(k). Then

    ord_ell(p)=2^(a+1) u,  u odd,  u divides k/2^a.      (17)

Indeed p^k=-1 makes the order divide2k but not k, with precisely this
2-part; equivalently its half divides k with odd quotient. Conversely
this order condition implies p^k=-1. The same argument works modulo
ell^e when ell^e divides p^k+1, since the cyclic unit group has the unique
element of order2 equal to-1.

If k is even, every odd factor ell is1 mod4. This agrees with the
sum-of-two-squares form (p^(k/2))^2+1, and excludes those factors as
vertices of the prime3-mod4 Legendre tournament. If k is odd, factors
of both residue types may occur; the order's 2-part remains exact.

**CITED existence.** Zsigmondy's theorem applied to (p,1,2k) provides
a prime with ord_ell(p)=2k for every odd p and k>1. At k=1 the only
exception is p+1 a power of2. The (2,1,6) exception is outside odd p.
The theorem and its exceptions are recorded in Avci,
[*Large Zsigmondy Primes*, introduction](https://arxiv.org/pdf/2011.06136).
This applies the classical theorem, without claiming any new prime-density
law. The repository's [older Zsigmondy triad](../collatz_mod6_20260917_zsigmondy_triad.md)
is retained for the fixed-base/exception distinction.

## G5. An all-depth signed Frobenius bank from those factors

Fix odd primes p,ell, ell!=p, with ord_ell(p)=2^(a+1)u_1, u_1 odd.
Set delta=v_ell(p^ord_ell(p)-1). The standard order-lifting law gives

    ord_(ell^e)(p)=2^(a+1)u_e,
    u_e=u_1 ell^max(0,e-delta).                         (18)

It follows either by the binomial/LTE valuation step for odd ell, or by
the inherited [all-depth order tower](../fermat_wieferich_order_20261014e.md).
Let r=p^(2^a) modulo ell^e and H_e=<r^2>, of odd order u_e.
Then r belongs to -H_e, and r^u_e=-1. On the additive vertices Z/ell^e,
choose a connection set S with S disjoint from -S and S union -S all
nonzero elements. Define x->y iff y-x belongs to S.

The following are equivalent:

    (i) H_e preserves S;
    (ii) multiplication by r reverses every tournament arc.

For (i)->(ii), r belongs to -H_e. For (ii)->(i), square the action.
Thus these are actual signed Cayley tournament actions. Multiplication
by p itself reverses such a tournament exactly when a=0 is possible;
when a>0, some even power of a proposed reverser would be negation, a
contradiction. Multiplication by p cannot preserve any such tournament,
since its cyclic group contains negation. Removing the2-part of its
order is essential, not optional.

Every nonzero H_e orbit is disjoint from its negative, because H_e has
odd order. On the valuation layer of additive order ell^j its orbit
length is u_j. Therefore the COMPLETE labelled bank has2^b_e states,
with exactly

    b_e=sum_(j=1..e) phi(ell^j)/(2u_j)
       =(ell^min(e,delta)-1)/(2u_1)
          +max(0,e-delta) ell^(delta-1)(ell-1)/(2u_1).  (19)

After the lifting delay delta, each additional level adds a constant
number of independent orientation bits. These are connection-set counts,
not graph-isomorphism counts. Translations and H_e give an automorphism
subgroup of order ell^e u_e; extra full graph symmetries are not excluded.

**Literal Frobenius map.** Choose a primitive ell^e-th root zeta in
F_(p^ord_(ell^e)(p)). Map exponent x to zeta^x. The field Frobenius
z->z^p becomes exponent multiplication by p, inversion becomes x->-x,
and differences become ratios. The selected connection set becomes a
set of allowed ratios between root-of-unity vertices. Thus Frob^(2^a)
is a converse map and Frob^(2^(a+1)) preserves the chosen orientation.
The field and primitive-root gauge are retained; this is not the
unrelated claim that multiplication by p is Frobenius on F_(ell^e).

For factors initially supplied by p^k+1, a=v2(k). As e increases, the
same family persists, with p^(2^a u_e)+1 divisible by ell^e. Hence (19)
grows a whole infinite family, not just the factors of one finite term.

Examples, all with the exact gauge above:

| base p | factor ell | a | delta | b_e |
|---|---:|---:|---:|---|
|3|7|0|1|e|
|5|7|0|1|e|
|3|5|1|1|2e|
|7|5|1|2|2 at e=1; 10e-8 for e>=2|
|5|13|1|1|6e|
|7|43|0|1|7e|

The 7/5 delay is visible already in7^2+1=50 and
v5(7^4-1)=2. At modulus25, r=7^2=-1 and H is trivial: all4096
circulant orientations are allowed. At modulus125 the bank has22 bits.
At p=3,ell=7,e=1 the two bank states are the two Paley7 orientations.

The field carrier matters even when the characteristic is unchanged.
Frobenius on the additive Paley F27 preserves its square differences
([the inherited literal check](../charged_recurrence_memory_20261028/FROBENIUS.md)).
Here Frobenius on the seven roots of unity in F729 reverses the chosen
ratio relation. The source vertices, difference-to-ratio map and primitive
root distinguish these operations; the word "Frobenius" alone does not.

## G6. Composition boundary and the hostile control

At child order28, the regular parent C3^[3] has group3^13, so every
child group is a3-group. The prime7 dividing28 cannot appear as a child
automorphism prime in this family. At order126, R5^[3] children have
only5-power automorphisms, despite126 having factors3 and7. At order344,
R7^[3] children have7-power groups; choosing P7^[3] instead permits3
and7. The factor43 from7^3+1 has a real role in G5's different signed
Frobenius tournament, not as a newly created symmetry in these children.

Retain the actual twelve-vertex81/27 pair from
[regular-parent R5](../regular_parent_symmetry_20261029/PROOFS.md): their
binary cuts have1620 versus1728 trivial CORE stabilizers. Neither parent
is regular. Therefore (2) and its unmarked conclusion are not applied
to this pair. The same cyclic triples, the same size, or leading moment
scales cannot replace the actual child types and group action.

No statement here closes arbitrary Hamilton-path complement/re-marking
operations, counts all tournaments at p^k+1, or provides an IID law for
unconditioned arithmetic attachments. G3 imports a measured IID palette;
using it on arithmetic data requires a separate joint realization law,
such as the scoped [growing-bank bounds](../growing_attachment_families_20261029/PROOFS.md).

For clarity, a sufficient weighted transfer guard is pointwise relative
control on the WHOLE leaf coloring: if |mu(w)/nu(w)-1|<=eta on every
nu-positive state, with no extra mu support, then

    |E_mu f-E_nu f|<=eta E_nu f,  f>=0.

This follows term by term and applies to f=|Stab|^s. If a separate
exceptional component has mass beta, its worst-case contribution is
beta |W_k|^s. Small total variation without this weighted control does
not transfer G3. No existing fixed-pattern arithmetic bound is silently
asserted to meet this growing whole-state guard.
