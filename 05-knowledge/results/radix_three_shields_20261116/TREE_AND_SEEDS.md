# A literal Fibonacci relation, 17 and 68, and the seed reduction

## The redundant ternary digit is a repair, not a new radix

The proposed two units of notation are equivalent: choose
d_j∈{0,1,2,3} at weights 3^j, then count shift 2Σ d_j3^j.
If independently realizable, these digits fill all even shifts
through 3(3^(k+1)−1). In half-step units the lower max before
digit j is 3(3^j−1)/2, so each new translate overlaps or abuts.
The canonical ternary decoder uses {0,1,2}; the fourth choice
3 is redundant, providing a carry by replacing one low 3
with one unit of the next digit. Nothing in the arithmetic
requires the same local tournament move to survive the carry.

With arbitrary leading zeroes let a(N) count representations,
a(0)=1 and a(−1)=0. For 0≤N<3^k the higher digits cannot be
nonzero, so this is also the k-digit count. The least digit
has two possibilities only for residue 0:

    a(3N)=a(N)+a(N−1),   a(3N+1)=a(3N+2)=a(N).

Read the adjacent-count state (u,v)=(a(N),a(N−1)) from high
trit to low. The trits 0,1,2 act respectively by

    R(u,v)=(u+v,v),    L(u,v)=(u,u+v),    J(u,v)=(u,u).

The first two are Stern's inverse-parent shears. The largest
multiplicity below 3^k is exactly F_(k+1), with F_1=F_2=1:
after any step the sorted coordinates are bounded by the
next Fibonacci pair, and alternating the shears after the
first 1 attains the bound, choosing the parity to end with R.
There are exactly 2^(k+1)−1 uniquely represented targets in
[0,3^k): zero and the strings with no zero after the leading
nonzero trit. Indeed an ordinary zero below a nonzero digit
permits a borrow ending in a digit 3 there, whereas a least
digit 1 or 2 is forced and allows descent without ambiguity.
In particular a forbidden adjacent pair of physical settings
among 1 and 2 can remove some targets despite the redundancy.

## The same shears really do occur on the primitive-triple tree

For Euclid parameters p>q>0 with gcd(p,q)=1 and opposite
parity, use the chart (u,v)=(p−q,q). Then u is positive odd,
v>0 and gcd(u,v)=1. The three Berggren/Barning–Hall moves
in this chart are exactly

    L,  R²,  R S R,   where S exchanges the two coordinates.

They act as (u,u+v), (u+2v,v) and (u+2v,u+v). Their images
are the regions v>u, u>2v and v<u<2v respectively. There is
no u=2v in this parity domain, and the only coprime u=v is
(1,1). Inverting the branch dictated by the inequalities
keeps positivity, parity and gcd and decreases u+v.
Thus this is an actual invertible, exhaustive ternary coding
on its full state, not a claim that projecting to a side
length classifies triples. In contrast J above forgets v;
tournament protocols still lack the analogous state-driven
parent for a prescribed count.

The initial triple from (u,v)=(1,1), (p,q)=(2,1) is (3,4,5).
Its R² child (u,v)=(3,1), (p,q)=(4,1) is (8,15,17) up to
leg order. This gives an exact shared matrix mechanism with
the width-four arithmetic rather than the slogan “both are
ternary.” Reference: [Price, The Pythagorean Tree: A New Species](https://arxiv.org/pdf/0809.4324).

## 17 and 68, with the coordinate specified

- There is just one positive right triangle with hypotenuse 68
  up to exchanging legs: (32,60,68)=4(8,15,17). Twice applying
  the even-square parity test divides the legs by 4, and at
  hypotenuse 17 there is just the Euclid pair (4,1).
  Separately (51,68,85)=17(3,4,5) uses 68 as a leg, not the
  hypotenuse. Both scalings preserve the cone but leave the
  primitive tree.
- The archived [signed Collatz cycle](../collatz_three_mirrors_20260929.md)
  contains the exact segment −68→−34→−17 under the shortcut
  even rule. Thus the ratio 4 is a dyadic valuation there.
  It is not a multiplication rule for a count dial.
- The published [67/69 exact seams](../conditional_determinant_vitali_20261116/PAIR_67_69.md)
  put 68 at the midpoint of two realized odd counts. 68 cannot
  itself be a tournament Hamilton count. Independently the new
  [live two-dial bank](TWO_DIALS.md) has one literal port weight
  68 named in its boundary table; that is not a path-count
  midpoint.
- The [new triangle wrapper](TRIANGLE_WINDOW.md) takes a
  literal 17-count seed to 241. In general if U_m reverses the
  endpoint arc across m intervening vertices in a transitive
  order, its only nontrivial cycles use that backward arc.
  Thus H(U_m)=1+2^m. For two optional path slots, assigning the
  singleton/odd-cycle components gives
  P_2(U_m)=4·2^m+2·3^m−2 (sum over odd subsets of the m sites).
  So the triangular wrap has H=5·2^m+2·3^m−1, giving 241 at m=4.
  The two values 17 and 241 also multiply to 2^12+1 and are
  [period-24 Riesel labels](../two_three_wheels_synthesis_20261116/WHEELS.md).
  That is an exact physical link, not a proof that their
  multiplicative orders cause the tournament construction.

## The prime and 7p reduction is exact

The primes in the statement must be **odd**: all tournament
Hamilton counts are odd. Necessity of the list follows immediately
from the full conjecture. For sufficiency H(A→B)=H(A)H(B), since
a path visits the two blocks in order. Write an odd integer as
7^e m with gcd(7,m)=1. If e=0 use the odd prime factors. If e
is positive and even use 49^(e/2). If e is odd and ≥3 use
343·49^((e−3)/2). For e=1, if m has a prime p other than 3,
peel 7p and factor the rest. Otherwise m=3^a: a=0,1 gives
7,21 and a≥2 uses 63 and remaining 3's. Count one is transitive.

The three special seeds are literal in the two-apex alphabet
2[A wins]+[B wins]: 49 is 10212 with B→A; 63 is 100110
with A→B; 343 is 1012000 with B→A. Their roles are distinct:
49 is the first pure even 7-valuation, 63 handles a lone 7
with only 3's, and 343 is the first permissible pure odd
7-valuation. No factorization of 343 into two factors >1
avoids 7. This is nothing like a primitive Pythagorean triple
arithmetically; its reduction depends on factoring and the
two forbidden totals, rather than an invariant quadratic form.

## A concrete absolute seed band for a future physical carry

A reusable rule on one closed class of prepared witnesses,
with exact child counts 3H,3H+2,3H+4, would do more than
produce gap intervals at one fixed scale. Starting from all
odd H∈[L,U] in that same prepared class, its depth-k children
fill the whole odd interval

    [3^k L, 3^k U+2(3^k−1)].

This follows either by ordinary ternary digits or by the
interval join for successive parent H. The intervals at
successive depths abut or overlap iff U+2≥3L.
The same seed threshold works if a fourth child 3H+6 is
admitted. For example L=23, U=67 is the first natural seed
band chosen above 21, and its required upper boundary is
precisely 67, beside the old 68 midpoint. All 23 scalar
counts in the band are known; it is not known that they
possess the same reusable preparation. With such a prepared
band, this exact physical carry would give all odd counts
≥23, and the remaining admissible odds have small literal
witnesses. This is a conditional absolute route, not a
consequence of the two finite dials or the isolated amplifiers.
