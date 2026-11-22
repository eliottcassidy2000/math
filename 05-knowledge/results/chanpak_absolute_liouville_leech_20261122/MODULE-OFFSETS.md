# Unequal module widths at the full perfect-number target

**Status: PROVED for the inherited module families; FINITE-EXACT for
the bounded source census.** These are obstructions to specific
physical carriers, not to arbitrary unequal star arms or the
Chan–Pak spectrum. The literal strong straddling edge at exponent
31 is not produced here.

Put E_p=2^(p−1)(2^p−1), with odd p≥5. No assertion below needs
2^p−1 to be prime. The two specific prescribed targets at p=31
are 2305843008139952127 and 2305843008139952129.

## 1. The four existing two-port neighbor templates

Use the exact compressed-word profiles in
[TWO_PORT_NEIGHBOR_FAMILIES](../chanpak_moonshine_prime_blocks_20261118/TWO_PORT_NEIGHBOR_FAMILIES.md).
Both transitive module widths m,k can be any positive integers,
including highly unequal. Nonetheless **none of their four
templates gives either E_p−1 or E_p+1** at any odd p≥5.

**Family I.** Its two counts are U±1 for
U=8·2^(m+k)+2^(m+1)−2^(k+1). This center is divisible by four,
so hitting either odd neighbor requires U=E_p; U=E_p±2 is
impossible. If m=k, U is a pure power of two. If m<k,
v₂(U)=m+1, so m=p−2. But then
U/2^(p−1)=1+2^(k−p+2)(2^p−1)>2^p−1. If m>k, v₂(U)=k+1
forces k=p−2 and m≥p−1, and the leading term alone is
at least 2^(2p), with positive remaining terms.

**Family II.** Its two counts are U±1 for
U=3(2^m−1)(2^(k+1)−1)+7. Modulo three the lower output
could only hit E_p−1 and the upper only E_p+1, hence U
must equal E_p. If m≥3,k≥2 then U=2 mod8.
If m=1, U=3·2^(k+1)+4, divisible by eight only for k=1,
with U=16. If m=2 then U=9·2^(k+1)−2; the remaining case
k=1, m≥3 gives U=9·2^m−2. Neither is divisible by eight.

**Family III.** Its three counts are N−2,N,N+2, with

```text
N=3·2^m+2^k B−3,      B=8(2^(m+1)+3^(m+1)−2).
```

The coefficient B has v₂(B)=3. Modulo 16 the only way to get
E_p+1 is N−2 at m=1. There N−2=88·2^k+1, and its odd part
11 after subtracting one is not 2^p−1. The only way to get
E_p−1 is N+2 at m≥4. The required equation is then
E_p=3·2^m+2^k B. Modulo 3 this forces m odd and k even:
for m even B vanishes mod 3, and for m odd B is one mod 3.
Now m<k+3 would force m=p−1, even; m>k+3 would force
k=p−4, odd. If m=k+3, divide the right side by 2^m to get
2^(m+1)+3^(m+1)+1, which has 2-adic valuation one, so p=m+2.
The equality reduces to 3^m=2^(m+1)−1; it is false modulo
eight for odd m≥5.

**Family IV (a literal single-arc switch).** Its two counts are U±1
for U=2^m+2^k(8·2^m+12·3^m−8). The coefficient of 2^k has
valuation two. If m<k+2, equality U=E_p forces m=p−1 and
k≥p−2, in which case even 12·2^k·3^m exceeds E_p;
the rest of U is nonnegative. If m>k+2, it forces k=p−3 and
m≥p, and the leading dyadic term alone is at least 2^(2p).
If m=k+2, then

```text
U/2^m=2^(m+1)+3^(m+1)−1.
```

For even m this has valuation one, making v₂(U) odd instead
of p−1; for odd m it vanishes mod 3 whereas E_p is a unit.
Thus there is no straddling U=E_p. The crossed cases U=E_p±2
need m=1 by reduction mod 4; then U=2+44·2^k. The equations
44·2^k=E_p or 44·2^k=E_p−4 fail respectively by odd part 11
and by 2-adic valuation.

## 2. The two additional inherited full-height offsets

The ten-vertex marked cores in
[the perfect/Wieferich packet](../chanpak_perfect_wieferich_seams_20261122/PROOF.md)
give E_p−5 and E_p−9 at particular widths. Allowing their two
widths to vary arbitrarily still does not produce either full
neighbor at odd p≥5.

- The first count is 128·2^(m+k)−46·2^m+38·2^k−5,
  always 3 mod4. To get E_p−1, modulo eight forces m=1 or k=1.
  In the first case 294·2^k=E_p+96; for p≥7 its valuation
  forces k=4, then 4704<E_p+96. At p=5 the valuation instead
  forces k=3, which also fails. In the second case
  210·2^m=E_p−72, forcing m=2. This asks E_p=912, between
  E_5 and E_6.
- The second count is 256·2^(m+k)−144·2^m+14·2^k−9,
  also 3 mod4. To get E_p−1, reduction mod16 forces k=2,
  giving 880·2^m=E_p−48. For p≥7 the valuation forces m=0.
  At p=5 it instead forces m=2, but 3520≠448.

Thus neither the inherited low offset nor a different pair of
module widths on the same core corrects the full target.

## 3. The exact coefficient glue for a whole new ray

Suppose a source with two marked vertices has no compressed words
visiting either mark more than twice. Its nonnegative integer
profile is (c_11,c_12,c_21,c_22), where c_ij counts words with i
visits to the first, j to the second, and all other vertices once.
With transitive modules, compressed words have respectively
1 or 2^m−2 fillings of the first mark and 1 or 2^k−2 of the
second. The resulting count is

```text
A·2^(m+k)+B·2^m+C·2^k+D,
A=c_22, B=c_21−2A, C=c_12−2A,
D=c_11−2(c_12+c_21)+4A.
```

For constant integer width offsets m=p−a and k=p−b, equality
to E_p±1 at infinitely many p is an equality of quadratic
polynomials in 2^p. Without loss a≤b, and the exact fiber is

```text
A=2^(a+b−1),    2^(b−a)B+C=−2^(b−1),    D=±1.
```

Conversely a physical source in that fiber realizes the full
ray at every p with positive widths. For example a=2,b=1
(up to reversing the port names) reduces to finding the signatures
(c_11,c_12,c_21,c_22)=(15 or 17,6,10,4) in the unsorted port
order with widths (p−2,p−1).

This gives an integral local test for a physical change of core.
Keep the same two widths, write r=2^(b−a), and require a flip
between two **unmarked** core vertices. Its signed c-difference
must belong to the affine line

```text
Δc_22=0,   Δc_21=t,   Δc_12=−r t,
Δc_11=2(1−r)t±2,          t an integer.
```

These conditions are also sufficient at the profile level,
since they fix all height coefficients and change exactly the
constant by two. The four c-to-(A,B,C,D) coordinates are related
by an integer unimodular change. Thus their primitive transverse
condition, e.g. C≡−2^(b−1) mod r when r>1, is genuine exact
glue; agreeing only on A or on a near-perfect height discards
it. Membership in this arithmetic fiber is not itself a
realizability theorem for compressed words or arc flips.

**Finite-exact search.** In the complete stored isomorphism-class
lists of strong sources through order eight, with any two
marked vertices, no profile of the stated max-two-visit type
lies in either ray fiber. The program computes the full
compressed signature by DP on used unmarked labels, two visit
counters and last label, and tests the coefficient equations
exactly. For its offsets it is enough to test |a|,|b|≤30:
for sources of order at most ten, A≤12!/(2!2!)<2^27,
and |B|,|C|<2^28. If A is a positive power of two then
a+b≤27. If a≥31 then b≤−4; when C is nonzero, its
term in B·2^b+C·2^a dominates the other term and A;
if C vanishes, equality requires |B|=2^(a−1)≥2^30,
also impossible. Swap ports for the other side, noting
a+b≥1. The same search at distance at most two arc flips
from the two displayed ten-vertex cores finds only their
existing −5 and −9 ray offsets among |D|≤25. Neither search
excludes higher visit multiplicity, widths other than p minus
a fixed offset, or larger/different strong sources.

The prior [
six-label Leech audit](../chanpak_mrb_leech_clock_20261114/MRB_LEECH_SOURCES.md)
distinguishes the historical chain 13→2→5→3→7→4 from the
artificial identification 4=13. Their pure-A Coxeter labels do
give rank-24 root systems:
4 corresponds to A_3^8, 2 to A_1^24, 5 to A_4^6, 3 to A_2^12,
7 to A_6^4, and 13 to A_12^2. Rank 24 and sharing a list of
labels do not identify the lattices: these root-system determinants
are h^(24/(h−1)), in particular 4^8 versus 13² at the two
endpoints. The usable design lesson of the established Leech
connection is to keep exact integral transverse glue rather
than only ambient rank and an approximate norm. The rank-four
coefficient lattice above is a literal exact glue test in that
sense, **not** a Leech-lattice identification or a realization
from its deep holes. The pure-A correspondence is discussed in
[Borcherds's Leech paper](https://math.berkeley.edu/~reb/papers/leech/leech.pdf).

## 4. A shorter scalar factor witness, not a straddling edge

The scalar P31 upper neighbor factors into
29·71·137·1621·5042777503 as in the star-arms route. The last,
largest factor has a one-apex realization on 34 vertices, using
the 33-bit contact word

```text
101101010011101011010110111010100
```

or two-sided address (23197,108692) with a central 0. Its count
is the exact weighted inversion read: the two end contributions,
one per adjacent 01, and 2^(i−j−1) for each earlier marked
position j and later unmarked position i. The first bit 1 and
last bit 0 ensure strongness. This reduces that scalar factor's
earlier 35-vertex witness by one; ordinal joining the factor
realizers still does not give a single-arc straddling edge or
a general all-exponent correction.

## Reproduction

From repository root compile the accompanying C++17 checker, then run
with argument 8 (the complete stored source deck) and with argument 0
(the two distance-two neighborhoods). The first prints the deck sizes
to stderr and no qualifying rays to stdout. The second prints eight
labelled instances, all with the original −5 or −9 coefficients.
The source deck sizes are checked so a missing file is not interpreted
as a negative search.

## 5. Why a scalar ordinal product is not a unit seam

Let U⇒V mean all arcs from a nonempty block U to a nonempty
block V, and reverse the single cross arc a→b, with a in U
and b in V. Write e for full paths of U ending at a, s for
full paths of V starting at b, and define

```text
P_U(a) = Σ_(A⊆U,a∈A) (# paths of U[A] starting a) H(U−A),
Q_V(b) = Σ_(B⊆V,b∈B) (# paths of V[B] ending b)   H(V−B).
```

Use H(empty)=1. The exact ordinary count increase is
P_U(a)Q_V(b)−es. The destroyed paths were precisely the old
U-to-V concatenations using a→b. A new path using b→a has
block pattern U_pre, V_pre (ending at b), U_post (starting
at a), V_post; each of the unmarked outer pieces may be
empty. Choose the two subset paths containing the respective
marks and the two complementary paths independently. Safe
forward concatenations are automatic: a and b already occur
at the reversed seam and cannot cause a second forbidden
crossing. Paths avoiding the changed arc are unchanged.

Removing the last a of a full path shows e≤H(U−a)≤P_U(a).
If a has an outneighbor x in U, the additional subset {a,x}
contributes at least one, so e<P_U(a). If a is a sink then
there are no longer paths starting a, and equality holds.
Moreover P_U(a)−e is **even**: apply the exact increase formula
with V a singleton and use Rédei's oddness on both tournaments.
Dually Q_V(b)−s is a nonnegative even integer, zero exactly
when b is a source of V. Thus reversing an arc between **two
nontrivial strong blocks** increases the count by at least
1·2+1·2+2·2=8, since e,s≥1 from cuts of directed spanning
cycles. The bound is attained by two directed triangles,
for which P=Q=3 and e=s=1.

More generally let the ordinal components be C_1⇒…⇒C_q.
For a single reversal to make the whole tournament strong,
its original tail must be in C_1 and its original head in
C_q; otherwise an untouched end segment blocks reachability.
For q≥3, no original Hamiltonian path uses that cross arc,
since it cannot skip the middle components. Splitting at
any component boundary, the new path formula therefore reduces to

```text
increase = P_(C_1)(a) Q_(C_q)(b) · ∏_(1<i<q) Z_(C_i)(2),
```

where Z_C(2) is the two-slot profile. Indeed the selected
path and the complementary path concatenate independently
within each chronological component; in each interior
component there are two unrestricted labelled slots. All
factors Z_C(2) are at least 2^|C| by Rédei on each part of
each splitting. In particular, when q≥3 the only case of
increase two is three singleton components: two or more
interior vertices already force increase at least four,
and if there is just one interior singleton but an endpoint
component is nontrivial strong, its pointed factor is at
least three (singleton, length-two, and full orders).

Consequently an ordinal-product witness with at least two
nontrivial strong components cannot be upgraded to a strong
realization at the same count or with a gap-two physical flip
by simply reversing a single arc. For q=2 this uses the
sharp bound eight; for q≥3 the preceding product excludes
gap zero and gap two unless all three components are singletons.
This excludes the direct one-flip upgrade of the scalar
perfect-neighbor factorizations, not a different strong host.

## 6. Positive result: a strong E_p−15 ray

Reverse the three unmarked arcs {0,4}, {5,6}, {6,7} of the
old E_p−5 source, keeping its marks 1 and 9. The new rows are

```text
(974,1020,1016,976,929,201,528,832,608,32).
```

Its compressed signature is exactly (c_11,c_12,c_21,c_22)
=(899,551,418,256), with no other nonzero entries, so it gives

```text
H=256·2^(m+k)−94·2^m+39·2^k−15.
```

Use unequal widths m=p−5,k=p−4, for any integer p≥6. The leading
term becomes 2^(2p−1) and the linear terms combine as
(−94+78)·2^(p−5)=−2^(p−1), so the exact total is E_p−15.
The quotient has directed Hamiltonian cycle
0→1→2→3→4→7→8→6→9→5→0; run consecutively through the
transitive blocks on this cycle to obtain a cycle of the inflated
tournament. Thus it is strong, with 2p−1 vertices. At p=31 the
widths are 26 and 27 and the 61-vertex count is
**2305843008139952113**.

The independent certificate is a Python memoized compressed-word
enumeration, plus literal Hamiltonian subset DP for p=6,7,8 and
strong-connectivity checks including the actual order-61 host.
Its normal and optimized outputs coincide. This is a new uniform
strong offset, not a straddling edge at E_p.
