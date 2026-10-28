# Perfect sources, prime divisor-count layers, and the appearances of 7 and 21

Status: PROVED elementary restricted classification. These are proper-divisor **term counts**, not values of an aliquot target or a tournament Hamiltonian-path statistic. No equivalence between those invariants is claimed.

## 1. Exactly two prime factors: a complete elementary classification

Write Delta(n)=s(n)-n. No prime power is perfect. If two distinct primes are both odd, then

\[
\sigma(p^a q^b)/(p^a q^b)<\frac p{p-1}\frac q{q-1}\le\frac32\frac54<2.
\]

For an even two-prime source n=2^a p^b, a,b>=1 and p odd prime, put S_b(p)=1+p+...+p^(b-1). Direct multiplication gives

\[
\Delta(n)+1=(2^{a+1}-p)S_b(p).
\]

Perfection forces the product on the right to equal1. Since S_b is a positive integer, S_b=1, b=1, and p=2^(a+1)-1. Such a p is prime only if r=a+1 is prime, by the usual factorization of 2^r-1 for composite r. Conversely every prime M=2^r-1 gives the perfect source n=2^(r-1)M by direct sigma multiplication.

This proves the classification for sources with at most two distinct prime factors without assuming a classification of all perfect numbers.

## 2. A complete theorem for each prime divisor-count parameter

**Theorem.** Let r be prime. A perfect source having exactly 2r-1 proper divisors exists iff 2^r-1 is prime. When it exists it is unique and equals 2^(r-1)(2^r-1).

Indeed tau(n)=2r. Every multiplicative partition of2r has at most two factors, so n has at most two distinct prime factors. The classification above applies. Its divisor count is 2r', where r' is the exponent in the Mersenne prime2^r'-1, forcing r'=r. This is a complete all-size argument for each such layer, not a scan.

The first cases:

| r | proper-divisor terms | perfect source |
|---:|---:|---:|
|2|3|6|
|3|5|28|
|5|9|496|
|7|13|8128|
|11|21|none, since2047=23*89|
|13|25|33550336|

This gives a precise role for11: it is the first prime exponent whose Mersenne candidate is composite, hence the first failure among these prime-parameter layers. The next entry25 connects the term count to8191, not to a putative prime-base90 representation:33550336=2^12*8191.

## 3. Seven terms fails for a different reason

Exactly seven proper divisors gives shapes p^7, p^3q, or pqr with distinct letters prime. The first is deficient. For p^3q,

\[
\Delta= p^3+p^2+p+1-(p^3-p^2-p-1)q.
\]

At p=2 this is15-q and cannot vanish for prime q. At p=3 it is40-14q and cannot vanish for integer q. For p>=5, even the lower bound q>=2 makes Delta negative.

For squarefree pqr with all primes odd, sigma(n)=(p+1)(q+1)(r+1) is divisible by8, while2n has exactly one factor2, ruling out perfection. Otherwise n=2pq, with odd distinct p,q, and perfection requires

\[
(p-3)(q-3)=12.
\]

The only positive even factor pair up to exchange is(2,6), giving p=5,q=9, which violates primality. If p=3 the product is zero. Thus **there is no perfect seven-term source**.

The numbers7 and21 therefore both occur as excluded term counts, but their proofs differ:7 exhausts three possible exponent shapes, whereas21 is a two-prime layer whose sole Mersenne candidate fails primality. They are not the only excluded counts. For example exactly11 terms gives tau12; this note makes no universal theorem about all other layers. Even the easy count2 is excluded, since its sources are prime squares. Numerical repetition alone supplies no map to the tournament exclusions of7 and21.

## 4. Checks and boundaries

`python3 check_perfect_layers.py` (also under `-O`) checks the prime-parameter theorem for r<=31 by exact factorization of the Mersenne candidate, verifies the named sources and their sigma/tau values, and independently tests every source through100000 with a divisor-addition sieve. Prime testing here is trial division on candidates no larger than2^31-1. The theorem's unbounded scope comes from the proof, not these finite checks.

A source equal to its aliquot target is a loop under the aliquot map. It is not an edge between distinct vertices in the simple square-sum graph. In particular6 and28 being perfect is separate from the square-sigma edges used in the1001→343→57 argument.
