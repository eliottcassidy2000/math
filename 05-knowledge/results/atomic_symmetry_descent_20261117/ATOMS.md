# The exact ordinal atom reduction for Chan–Pak

Let M denote the maximal admissible arithmetic monoid, positive
odd integers except 7 and 21, and let S be the (unknown) physical
Hamilton spectrum. Since the ordinal join multiplies H and 7,21
are forbidden, S is a submonoid of M. The atom list below is for
M, not an unconditional list of every atom of S.

## Exact generators and physical consequences

The atoms of M are:
- each odd prime other than 7;
- 7p for prime p>=5 other than 7;
- 49,63,343.

Each has no proper product in M: prime is immediate, 7p and
49 need a factor 7, 63 has only the nontrivial splits 7*9 and
3*21, and 343 has only 7*49. Conversely, write a target as
7^a times a 7-free odd number. If a>=2, write a=2x+3y and
use 49s,343s and the ordinary prime factors. For a=1, peel a
7p if there is a prime p>=5; otherwise peel 63 since the
remaining target is 3^b with b>=2. The case a=0 factors freely.

Literal strong eight-vertex witnesses for 49,63,343 are in
check.py. Thus Chan–Pak is equivalent exactly to coverage of
the prime and seven-times-prime atom families above. In particular
the least additional hole, if any, belongs to one of these two
families. This does not exclude additional reducible-in-M targets
being physically primitive if their required M-atom factors are
unattained.

Every actual tournament with H equal to an M-atom has exactly
one nonsingleton SCC. All others are singletons before or after
it in the unique SCC order. For an atomic h let s_h(k) count its
strong order-k classes. Then its complete count fibre is

    a_h(n) = sum_{k<=n} (n-k+1)s_h(k).

For a converse-fixed class the strong core must itself be
self-converse and its left and right pads must have equal length;
thus a_h^sc(n)=sum_{k<=n, k≡n mod2}s_h^sc(k). This is an exact
class-level classification even if h never occurs.

## Factorization lengths, sharp elasticity and gaps

Write an integer in M as N=7^a 3^b prod_{p>=5,p!=7}p^{c_p}.
An atomic factorization is determined by x=#49, y=#343, z=#63,
and t_p=#7p, together with the implicitly remaining plain
primes. The complete equations are

    2x+3y+z+sum_p t_p=a,   2z<=b,  0<=t_p<=c_p,
    L=b+sum_p c_p+x+y-z.

All nonnegative solutions are valid, and all come this way.
This describes arithmetic decompositions into M-atoms, not all
possible physical SCC decompositions.

The ratio of maximum to minimum positive length is at most 5/2
for every N>1: assign additive weight v7/2 plus the sum of all
other prime valuations. Atom weights lie between 1 and 5/2.
Equality occurs at 63^2 = 49*3^4, where the lengths are 2 and 5.
For every m>=2 the length set of 63^m has first two entries m
and m+3: taking z=m is the minimum, z=m-1 is impossible, and
the next choice z=m-2 adds one 49 and gives m+3.

The union of adjacent-length gaps over M is exactly {1,2,3}.
To see no larger gaps can occur, any two factorizations of a
fixed N connect via the following bounded-length rewrites.
Pair any two atoms using a single factor 7 into 49 plus their
plain prime tails, changing length by at most 3. Replace pairs
of 343 by triples of 49. At this point at most one single-7
atom remains. If a>=2, absorb it with any 49 to obtain 343
plus its plain prime tail, or with 343 to obtain two 49s plus
its tail; each changes length by at most 2. This leaves the
canonical decomposition with at most one 343 and only plain
primes and 49 otherwise. If a=1, canonicalize to the least
available p>=5 using (7p)q <-> (7q)p and
63*p <-> (7p)*3^2, or use 63 when no such p is available.
Connecting via canonical form gives a sequence of jumps of
at most 3. The exact examples of gaps 1,2,3 are respectively
35^2, 63*35 and 63^2 by the displayed equations.

For a=0 or a>=2, suppose each non-7 prime divisor of N is
physically realizable. Then among **all** actual tournaments
with H=N, the maximum number of nonsingleton SCCs is
b+sum c_p+floor(a/2). For the upper bound refine every SCC
count into M-atoms (each count is admissible). For equality
use known 49,343 and all ordinary prime factors; no 7p
coverage is needed. For a=1 with some c_p>0 the upper
bound b+sum c_p is attained if and only if one of those
7p is realized (given the other ordinary primes). Indeed
an arithmetic maximal-length split at a=1 must have one
7p and no 63, and any genuine maximal split must already
be atomic. If no such 7p exists but b>=2, the next lower
bound is exactly attained using 63.

## A prime-target bridge with the compensated ×7 reader

Use the exact scaled pair from the prior cyclic-switch note,
whose lower (v,S)=(56551,3390208000) and whose upper profile
is seven times this pair. An arbitrary common ordinal
attachment C multiplies both profiles by (x,T), so the
upper output remains exactly seven times the lower one.

Choose a physical A0 of profile (u,R0) and physical C0 of
profile (x,T0). Adjoin k ordinal singleton vertices to A0
and l to C0, giving R=2^k R0 and T=2^l T0. The exact lower
target equation is

    p - 2vxu = (2^k R0-u)(2^l S T0-vx).

Consequently a prescribed p selects an eligible divisor d
of p-2vxu (when positive) iff (d+u)/R0 and
((p-2vxu)/d+vx)/(S T0) are powers of 2 with nonnegative
exponents. A hit gives both p and 7p from explicit graphs.
This is an actual Euclid-style shift and factor selection,
not just a residue condition. For fixed A0,C0 this branch
produces only O((log X)^2) counts up to X: the two factors
are at least 2^(k-1)R0 and 2^(l-1)S T0 respectively.
It does not establish that every target has an eligible
divisor for some chosen pair of cores.

Further exact primitivity sieve: write the lower readout
as p=alpha R+beta u with alpha=S-v, beta=3v-S (replace v,S
with vx,S T when appending C). Since u,R<p in any such
cyclic readout, primality demands gcd(u,alpha R)=1 and
gcd(R,beta u)=1. Thus the SCC constituents within A must
have disjoint cross prime supports between their H factors
and their P(2) factors. In particular 3 divides no H factor,
as P(2) ≡ ±H mod 3. A rich source built from powers of the
triangle cannot occur on this prime branch; primitive
profile collections with disjoint supports, such as
(11,100) and (17,448), can.
