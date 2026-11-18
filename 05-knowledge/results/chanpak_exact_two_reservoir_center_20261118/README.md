# Exact two-reservoir centers with a signed switch

This is an all-order construction and an inverse for its own language,
not a proof of Chan–Pak. It follows the
[factorial and floor session](../chanpak_factorial_floor_20261118/README.md).
Use H for Hamiltonian paths and P_s for ordered spanning path covers
in s slots, empty slots allowed.

## Two independent profiles give the center and the offset

Take two 5-vertex templates with outgoing bit rows, starting at index 0:

    B- = (12,17,2,6,13)
    B+ = (8,17,3,6,13).

They differ only on 0--2. In each template replace vertex 0 by an
arbitrary nonempty tournament X, and vertex 4 by an arbitrary nonempty
tournament Y. Cross arcs follow the template. Both hosts are strong:
the base has a directed Hamiltonian cycle (for B-, 0,3,2,1,4,0;
for B+, 0,3,1,4,2,0), and the quotient paths/cycles route between and
within all modules.

Write h=H(X)H(Y), k=P_2(X)P_2(Y). The two literal counts are exactly

    H(B-[X,Y]) = 3k-h,     H(B+[X,Y]) = 3k+h.

Proof: contract consecutive visits to a module to one letter. The
nonempty-run multiplicities (1,1), (2,1), (1,2), (2,2) have skeleton
counts (11,6,6,3) for B-, (13,6,6,3) for B+, with no others. To check
this small finite identity, enumerate words containing 1,2,3 exactly
once, 0 and 4 at least once, following the specified arcs. Once the three singleton letters are placed, each of the
four gaps has at most one 0 and at most one 4: the arc between 0
and 4 cannot reverse. That bounds the enumeration without a width
assumption. Each skeleton is filled with independent ordered nonempty covers of X
and Y, and Q_2(Z)=P_2(Z)-2H(Z). Expansion gives the displayed formula.
If X is a singleton the hosts differ in a single arc, even when Y is
arbitrarily large. No native metagraph adjacency is asserted.

In particular, X singleton and Y transitive of width t realizes the
exact one-arc pair 3*2^(t+1) ± 1 for every t>=1. These are absolute,
unbounded +2 pairs. They bypass the fixed-module signed calibrator
rather than claiming that calibrator becomes a universal +2 repair.

A single output can point upward in H and downward in a continuation:
when X is singleton the change in the two-slot count is
`4 Q_1(Y)+2 Q_2(Y)-2 Q_3(Y)` (two labelled slots, compress visits
within each slot, with the separator preserved). For the four-vertex
source D=(2,4,9,3), this is -8, while H changes by +10.
Thus carrying only H would lose even the direction of the next update.

## Exact inverse and a target-addressable infinite family

For a target N and sign e in {-1,+1}, choose the possible signed offset
first. Since each P_2 includes two one-slot covers, k>=4h; hence h is
positive odd and at most floor(N/(12+e)). For each such h:

1. Reject unless K=(N-e h)/3 is a positive integer divisible by four.
2. Factor h=h_X h_Y and K=p_X p_Y with each p_i even and >=2h_i.
3. Check that each pair (h_i,p_i) comes from one tournament, not from
   independent witnesses for its coordinates.

Step 3 has an exact finite test. Remove singleton strong components
from the candidate profile. The non-singleton components are ordered;
their H and P_2 multiply under ordinal sum. For H=h_i, their total
order is at most (9/2)log_5 h_i by
[Busch's lower bound](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/EJC/Volume_13/PDF/v13i1n3.pdf):
each of order n>=3 contributes at least 5^((n-1)/3) and
n<=(3/2)(n-1). Appending t singleton components multiplies P_2 by
2^t without changing H. Test the quotient for an actual power of two
and append that many singletons. For h_i=1, the empty core is allowed
but at least one singleton is needed. This decides the two-reservoir
language, not the full spectrum, and is not claimed to run efficiently
with unbounded target. For prime N, gcd(h,K)=1 is necessary. Modulo 3,
P_2(Z)=(-1)^|Z| H(Z), by odd-cycle parity of its slot polynomial;
a prime target other than 3 cannot use a reservoir with H divisible by
3. Neither test replaces profile existence.

There is an explicit family with independent absolute address.
The four-vertex D above has (H,P_2)=(5,32). Use a total of j copies
of D and t singleton components in the two ordinal reservoirs. Any
k>=2 and 0<=j<=floor((k-1)/5) can be implemented with X singleton
and Y containing all j copies and k-5j-1 singletons. Thus, exactly
at center 3*2^k, the following menu is certified:

    3*2^k ± 5^j,       k>=2, 0<=j<=floor((k-1)/5).

The two hosts for each j differ by one arc. At a fixed k the offsets
are distinct; across k the enclosing intervals
[(3-1/4)2^k,(3+1/4)2^k] are disjoint, because
5^j<=2^(k-2) for these parameters (for j>=1, 5^j<=2^(5j-1)).
So for this sublanguage the target chooses k by its interval and j by
a power-of-five test on the exact residual. It still provides only
on the order of (log N)^2 counts up to N.

## Why calibrating first cannot solve arbitrary absolute placement

For a fixed finite bank of local module substitutions with
P_s(M+)-P_s(M-)=c s and nonzero c, the gap is c h where h=H(T)
for the contracted tournament T. The only parameters are local
module choice and context. After deleting singleton strong components
of T other than the marked vertex, both endpoints are unchanged,
and the compressed order is <=1+(9/2)log_5 h as above.

In fact this gives a power obstruction, not just finiteness at fixed
gap. If the bank modules have at most L vertices, remove the marked
vertex and call the remaining compressed context U. A path in an
inflated endpoint orders the <=L new vertices and puts each old
vertex into one of the <=L+1 intervals around them. Ignore constraints
at the new vertices. For each partition of U into those intervals,
the product of induced Hamilton counts is <=H(U): under the
[odd-cycle packing expansion](https://arxiv.org/abs/2307.05569),
the tuple injects weight-preservingly into packings on U. Likewise
H(U)<=h by adding a singleton to packings. For either endpoint N,

    N <= L! h^(1+(9/2)log_5(L+1)).

Consequently any nonzero gap in this bank grows at least as a fixed
positive power of its absolute endpoint N. In particular the
one-arc unit rays of Section 1 cannot be explained by a bounded
local monomial switch with arbitrarily large exterior. They are
powered by a growing center on a different slot coordinate.

## How the open-problem analogies change the search

- **Collatz inverse trees.** The ordinary accelerated odd preimage
  is (2^a y-1)/3 with integrality and oddness constraints. The earlier
  [Terras inverse-clock work](../collatz_terras_inverse_clock_bridge_20261007.md)
  distinguishes a return in a clock from an integer landing. Here a
  candidate dyadic carriage likewise demands an exact power of two,
  not congruence alone. Tao's
  [almost-all Collatz theorem](https://arxiv.org/abs/1909.03562)
  has a different scope; neither problem transfers to the other.
- **Binary Goldbach.** The analogy is at the quantifier: local
  plausibility must produce both summands at a fixed target, not
  independent lists with good residues. The
  [Green–Tao linear-forms framework](https://annals.math.princeton.edu/wp-content/uploads/annals-v171-n3-p08-p.pdf)
  explicitly separates binary Goldbach from its finite-complexity
  setting. Our corresponding coupled question is Step 3, where
  projecting to H and P_2 separately loses the connection between
  factors. Unlike Goldbach, the summands here are products of
  coupled tournament profiles, and no implication is claimed.

**Next proof target:** find a second physical selector for K that
changes the odd part of the center while retaining an offset that
is already realizable; track at least P_2 and P_3 after the switch.
There is no universal lift from the inverse sieve yet.

### A literal way to change the odd part of the center

Two six-vertex reservoirs have the same H=23 but different P_2:

    (28,53,56,50,32,1) -> (23,272) = (23,16*17)
    (14,60,56,48,33,1) -> (23,296) = (23,8*37).

Take X singleton and let Y be one of these reservoirs followed by
singletons. For every k>=5, the resulting centers are 51*2^k or
111*2^k, both with exactly the same offset 23. Choose the reservoir
from (N±23)/3 by its odd part and the padding by its 2-adic valuation;
this is a literal two-phase selector for these rays. Adding j copies
of D to Y instead changes the offset to 23*5^j, and preserves either
center after spending 5j powers of two already in its address.
This gets a physical change in the odd part without losing the
repair. It does not place a center on every possible odd part.
