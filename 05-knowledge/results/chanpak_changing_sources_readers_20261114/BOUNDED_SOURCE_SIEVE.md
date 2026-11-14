# A Fermat-block sieve against all bounded-count cyclic sources

**All orders, conditional only on using the indicated reader interface.**
Let B be any nonempty tournament, of unrestricted order, with
h=P_B(1), Q=P_B(2), W=h+Q. In the one-seam construction with cyclic
source C3[1,1,B] and an apex-address reader, the terminal count is of
the form

    N = W m − h ε,   m≥1, ε∈{0,1}.

See the [apex-reader construction](APEX_READER.md) for the more
restrictive *literal* shape of m and the parity that determines ε.
The theorem below even allows arbitrary m and either ε independently.
It therefore cannot confuse an address phase failure with a failure
of this broader algebraic interface.

**Theorem.** For every fixed H0, the union over all literal modules
B with h≤H0, all positive m, and both ε, misses a positive lower
density of the odd integers. This allows B to depend on N, an
unbounded source order, and an arbitrarily long address. Thus fixing
a finite list of source counts is not enough even if the source
tournaments themselves are not fixed.

## Only the singleton budget is unbounded at fixed source count

Let the strongly connected components of B occur in their ordinal
order. Both P_B(1) and P_B(2) multiply, because every color-path runs
through the components in order and may be empty in some. A singleton
multiplies Q by 2 and leaves h unchanged.

The nontrivial strong components have at least three Hamilton paths.
Their number is consequently bounded by log_3(H0). Their individual
orders are bounded by [Busch's strong-tournament minimum](
https://doi.org/10.37236/1141); only the fact that the minimum grows
without bound is used. Thus there are finitely many nontrivial
component words. For each such word let (h,q) be its two evaluation
products. Including all singleton placements gives exactly a
geometric progression of possible Q at this word:

    Q = q 2^a,   W = h + q 2^a.

Here a≥0, except the empty nontrivial word requires a≥1 because
B is nonempty. All W are odd and at least 3: P_B(1) is odd and
P_B(2) even by the odd-cycle packing formula. In addition h+1<W
since Q≥2. Multiplicities and exact isomorphism counts play no
role in this reduction. The same SCC compression is used in the
older [repunit degree work](../../../01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md).

## A progression survives all tails in positive density

The sum of (h+2)/W over this finite family of geometric rays is
finite. Pick a finite prefix of each ray so that the sum over tails
is smaller than 1/2. Let L be the lcm of all prefix W, and restrict
targets to N=1 mod 2L. No prefix W divides N or N+h. This is
automatic since the residues are 1 and h+1, strictly between 0 and W.

For a tail W, count in this same target progression. Its condition
N=0 mod W has no solutions unless gcd(2L,W) divides 1. If it does,
the relative density is 1/W. The shifted condition N=−h mod W
has no solutions unless the gcd divides h+1. If it does, the
relative density is gcd(2L,W)/W≤(h+1)/W. So the two conditions
jointly consume at most (h+2)/W of this progression.

This is also an elementary finite count, not a countable
subadditivity assumption on natural densities. For N≤X, only
W≤X+H0 can occur (since m≥1). There are O(log X) such rows across
the finitely many rays. Count each intersecting progression with
its relative density and additive O(1) rounding error, and take
a union bound. At least (1/2−o(1)) of the tested progression escapes
as X→∞. Hence the lower asymptotic density of omitted odds among
all odd integers is at least 1/(2L). The proof gives an effective
finite prefix cutoff once the strong component inventory needed
for H0 is supplied.

## A compact Fermat sieve when h=1

Here B must be transitive (any directed triangle contributes
at least 2 to H in addition to the empty packing), so W=2^a+1,
a≥1. For s≥1 put L_s=2^{2^s}−1, the product of the first s
Fermat numbers, whether prime or composite. On N≡1 mod 2L_s,
if 2^s does not divide a then 2^a+1 is divisible by the Fermat
number F_{v_2(a)}=2^{2^{v_2(a)}}+1, whereas neither N nor N+1
is. If 2^s divides a then W≡2 mod L_s and gcd(W,2L_s)=1. The
remaining relative upper density is at most

    2 sum_{j≥1} 1/(2^{j 2^s}+1) < 2/(2^{2^s}−1).

For s=2 this already leaves at least 13/15 of the progression
N≡1 mod 30 outside the transitive-source interface. This is a
testable absolute-target desert, not an asymptotic estimate on
the literal apex word lengths.

## Consequence for choosing a source

A universal theorem through this interface must repeatedly change
the h-coordinate without bound, or it must bring in readers whose
offset is not one of 0 and −h. Neither adding more singleton
vertices nor freeing the address phases can repair this. The
factorial lock in [FACTORIAL_LOCK.md](FACTORIAL_LOCK.md) explains a
different restriction: after preserving enough low coordinates,
higher moments come with fixed absolute divisors. Together these
point to a change of low-order physical profile (not just profile
complexity or order) or to a reader with a different seam.
