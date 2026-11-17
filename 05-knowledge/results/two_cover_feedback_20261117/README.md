# Selecting the smaller two-cover block: frozen fibers and delayed choices

**Proved.** In the inverse 15u + 62R = N, a physical block bank can keep
the two-cover count R exactly fixed while selecting any unit residue
for u to arbitrary precision at 3, 7, and 31. A different pair of blocks
keeps *both* u and R fixed while storing as many choices as desired for
the next application of the wrapper. These are exact constructive
operations, not absolute coverage. In fact 217 is an obstruction to
the entire previous wrapper family, not merely to a bank of small inputs.

## 1. Freeze the cover count and vary the Hamiltonian count

Write P_B(q) for the optional q-slot path-cover polynomial, u = P_B(1)
and R = P_B(2). The [two-run wrapper](../two_run_ready_operator_20261117/README.md),
followed by two near-sinks, reads 15u + 62R. Ordinal juxtaposition multiplies
the whole polynomial: in each slot, concatenate the old paths in block
order; splitting by blocks inverts this operation.

| Row masks (outgoing arcs) | u | R |
|---|---:|---:|
| (2,4,9,3) | 5 | 2^5 |
| (40,1,3,6,15,30) | 19 | 2^8 |
| (2,40,11,49,71,213,143,31) | 223 | 2^12 |
| (96,21,25,51,225,198,142,15) | 599 | 2^13 |
| (0) | 1 | 2 |

Assign costs 5, 8, 12, 13 to the first four tiles. Realize a product
u = 5^i 19^j 223^k 599^l at *any prescribed* R = 2^e with
e ≥ 5i + 8j + 12k + 13l by padding with singleton tiles. This upgrades
the earlier trade P_C4(2) = P_TT5(2): the swaps form a commuting bank
and can change independently of placement.

**Frozen-fiber selection theorem.** Fix m = 7^a 3^b 31^c with a,c ≥ 1
and b ≥ 2. For all sufficiently large e, *every* unit class mod m
occurs as H(B) mod m with P_B(2) = 2^e in this bank. A finite BFS under
multiplication by 5,19,223,599 modulo 7·9·31 = 1953 reaches all 1080 units.
The replay supplies predecessor words. For the lift, the valuations of
5^30 − 1 at 7,3,31 are 1,2,1. Thus 5^30 generates the principal kernel
at each of those base depths: raising 1+p^v w to p^s gives a difference
from one of valuation v+s for the odd primes here. The three kernel
orders are pairwise coprime so one exponent suffices by CRT. The tile
cost is bounded for the fixed modulus, and singletons spend any leftover
budget.

**Applying it to a target.** Let N = 7p with gcd(p,93) = 1. For
M = 7^a 3^b 31^c, a,c ≥ 1, b ≥ 3, choose e mod 6 so N − 62·2^e
is divisible by 3 but not 9: the powers of two run through all six units
mod 9. Then the inverse specifies a unit modulo M/3:

    u_0 = [(N − 62·2^e)/3] · 5^(-1)  (mod M/3).

It is a unit at 7 since 7 divides N, and at 31 since 31 does not divide N.
Every sufficiently large e in the chosen phase can use the frozen-fiber
theorem. The replay fixes an enhanced congruence on e so that it can
construct the small part before spending the remaining padding.

## 2. Exact inverse on the bank and the all-order scale obstruction

For any physical B, R ≥ 2u (put a whole path in either slot).
With dyadic R the exact target satisfies

    62·2^e < N ≤ (62 + 15/2)·2^e.

Thus e = floor(log2(N/62)) is forced. Divide N − 62·2^e by 15 exactly,
peel the four prime factors, and check their cost. This is a complete
exact inverse **for this bank**. For example, the prescribed
N = 7·2,377,693,121 = 16,643,851,847 selects R = 2^28 and
u = 56,905 = 5·19·599, needing two singleton tiles. This gives a literal
20-vertex input and a 27-vertex finished tournament. It is an example of
exact backward selection, not a claim the value had never been hit.

No richer catalogue of *dyadic* R solves placement inside this wrapper,
even allowing arbitrary numbers t of near-sinks. When R = 2^e,

    14·2^(e+t) < (4·2^t−1)u + (16·2^t−2)R < 18·2^(e+t).

Drop the positive u term for the lower bound and use u ≤ R/2 for the
upper. The gap between adjacent scales contains [18·2^s,28·2^s].
By the prime number theorem it excludes infinitely many bona fide 7p
targets. Modular solvability in the frozen fibers does not bridge it.

**A non-dyadic bridge.** Add the literal tile with rows (24,1,3,6,14)
and profile (11,100). For fixed M take L = phi(M). The two repeatable
operations (add L singletons, or add L copies of this tile) preserve
both u and R modulo M while multiplying R by 2^L and 100^L,
respectively. The log ratio is irrational by unique factorization.
For any positive tolerance, finitely many multiples of log(100)
mod log 2 form a net of the circle. Thus one can choose the exponent
of the 100-loop from a bounded set (depending on tolerance and M),
then the exponent of the 2-loop to put 62R in any prescribed fixed
relative interval near a sufficiently large target. In the read
62R+15u the ratio u/R tends uniformly to zero because the singleton
exponent tends uniformly to infinity. Place 62R strictly inside a
smaller relative interval below or above N to accommodate the
positive 15u. Consequently for each fixed M and epsilon > 0, every
sufficiently large eligible N has physical reads in
[(1−epsilon)N,N] and [N,(1+epsilon)N], both *exactly* congruent to N
mod M. There are only finitely many base residue classes, so the
threshold is uniform in N. This removes the fixed relative gaps
while preserving the residue. It does not give absolute error < M.
A fixed finite product bank has only polylogarithmically many profiles.

## 3. Delay the choices, then let a second wrapper read them

Two order-eight blocks with outgoing rows

    A = (4,21,208,39,105,71,139,59)
    B = (128,129,27,51,99,71,15,124)

have (P(1),P(2),P(3)) equal to (345,5640,45225) and (345,5640,44505).
For any s ≥ 1, take an ordinal sum of s blocks, using i copies of A
and s−i copies of B. Its first two reads are always 345^s and 5640^s.
The first wrapper gives the same full rooted (H,g,E,O) for all s+1
choices. But their new P(2)'s are

    48·5640^s − 26·345^s + 86·45225^i 44505^(s−i).

All are distinct, and a second wrapper reads them into distinct
absolute Hamiltonian counts. Any common ordinal factors can be added
before the first wrap: they multiply the corresponding slot
polynomial values. This is an unbounded collection of different future
reads at one initialized profile, not a fixed reader menu. At s = 1
the second wrap changes H by 866,880, so the operation is not a unit
repair. The analogy with a carry register is literal about information
loss: ordinal composition stores a higher read invisible to the first
projection, and the second wrapper exposes it.

## 4. Why change the wrapper as well as the source?

N = 217 = 7·31 is impossible for the inherited wrapper with any B,t,
even for non-dyadic covers. This obstruction is from the parallel
inverse lane. The bound H_t ≥ 36·2^t−5 eliminates t ≥ 3. At t = 2,
u ≤ 1 but 31 divides u. At t = 1, u ≤ 3 and 7 divides R. Count 1 is
transitive; count 3 is a triangle with ordinal singleton padding.
Their R's are respectively powers of two and 12 times powers of two.
(The count-3 classification also follows from the [odd-cycle packing
expansion](https://arxiv.org/abs/2307.05569): there is exactly one
directed triangle; any strong block of order at least four has two
or more triangles, from vertex pancyclicity.) At t = 0, u ≤ 7 and
7 divides u, impossible. In contrast the one-apex contact word
100110010 on the transitive spine gives 217 exactly.

A sharper next selection question is to find growing, non-dyadic
iso-R trades with controlled signed u differences, or a nonordinal
loop that preserves the residue while shrinking absolute error at each
step. The delayed cohorts supply continuation multiplicity, not its
necessary spacing. Chan–Pak remains open.

Replay: `python replay.py` (also works with `-O`). Finite certificates
for the tiles and the collision use literal subset DP and optional
slot convolution independent of the wrapper transfer formula.
