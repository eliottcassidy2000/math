# Three normalized low bits from a divisor query

**Status: PROVED IN NOTE.** The proof below is elementary from the inherited
Burnside formula. **FINITE-EXACT:** independent full rational Burnside sums
at every even order 4 through 40. No external priority claim. This is a residue
query, not a full-integer counting algorithm or a Chan–Pak construction.

## Inheritance and research contract

Sources: THM-514 (`the-a000568-unit-tail-core-kernel-and-prime-tail-ladder`),
THM-516 (`burnside-core-kernel-phi-reframe`), THM-305
(`two-adic-valuation-tournament-count`), and especially THM-5098
(`precision-cycle-cutoff-and-normalized-tournament-congruences`). The latter
already proves the precision-dependent O(sqrt(b)) cycle cutoff and odd-order
mod-four/mod-eight congruences. Thus neither a low-cycle cutoff nor odd-order
valuation is claimed new here. Historical April computations
`04-computation/v2_merge_s24c.py` and `prove_hyp1712_s24c.py` already study
coprime two-cycle reciprocal sums. We extend the precision-first viewpoint
with an explicit even-order three-bit formula and a divisor-sum evaluator.

Anchor: exact A000568 queries. Niche: unit-square identities modulo 8.
Wildcard: arithmetic access to huge combinatorial indices. Concept board:
cycle multiplicities, valuation cutoff, unit-square algebra, truncated reduced
residue sums, divisor convolution. The hostile is repeated cycles: (1,1,1,1)
and (3,3) require exceptional corrections. The least-used sidecar is the sum
of representatives in the lower half of a reduced residue system.

## The theorem

Write a(n)=A000568(n). For even n>=4, U_n=a(n)/2^(n/2) is an integer.
Define

    c(n)=phi(n)/2,
    A(n)=sum_{1<=a<n/2, gcd(a,n)=1} a,
    D(n)=[3 divides n and n>=12] phi(n/3)/2.

Then

    U_n = n A(n)-c(n)+4D(n)+6[n=4]+2[n=6] (mod 8).       (1)

Furthermore, if r_d=floor((n-1)/(2d)),

    A(n)=sum_{d|n} mu(d) d r_d(r_d+1)/2.                  (2)

Only squarefree divisors need be visited. Given a prime factorization, this
uses O(2^omega(n)) arithmetic operations on O(log n)-bit integers, with
additional bit costs for the elementary arithmetic. Factoring n is NOT
included in this bound. Thus the routine accepts a factorization explicitly.
It does not enumerate partitions, tournaments, or all a<n/2.

### Proof of the surviving strata

For an odd partition lambda of n, let t be its number of cycles, m_l its
multiplicities, S=sum_l popcount(m_l), G=sum_{i<j}gcd(l_i,l_j), and z its
centralizer. The contribution 2^e/z to U_n has two-adic valuation

    delta=G-3t/2+S,

using e=(n-t)/2+G and v2(m!)=m-popcount(m). Since n is even, t is even.
For t>=6, G>=t(t-1)/2 and S>=1 imply

    delta >= t(t-4)/2+1 >= 7.

For t=4 the multiplicity patterns give:

- 1+1+1+1: delta=G-2>=4.
- 2+1+1, repeated length a: delta>=a+2>=3.
- 2+2, distinct repeated lengths a,b: delta=a+b+4gcd(a,b)-4>=4.
- 3+1, repeated length a and other length b:
  delta=3a+3gcd(a,b)-3>=3.
- 4, common length a: delta=6a-5, below 3 only at a=1.

That last exceptional partition at n=4 contributes
2^6/(4! 2^2)=2/3=6 (mod8).

For t=2 and distinct odd lengths a<b, write g=gcd(a,b). The normalized
contribution is 2^(g-1)/(ab). It survives modulo8 exactly for g=1 or g=3.
For equal lengths a,a it is 2^(a-2)/a^2, surviving at a=3 (n=6), where it
is 2/9=2 modulo8. The a=1 exception is n=2, outside the theorem's domain.
Every term is two-adically integral for n>=4. Since a(n) is an integer,
this proves the stated ordinary integrality of U_n.

The coprime two-cycle terms run over a<n/2 with gcd(a,n)=1, since n is even.
There are phi(n)/2 such a. Every odd unit u satisfies u^2=1 modulo8, so

    (a(n-a))^(-1) = a(n-a) = na-1 (mod8).

Their sum is n A(n)-c(n). The gcd-three terms are four times an odd unit;
dividing both lengths by 3 identifies them with the coprime unordered pairs
of sum n/3. There are phi(n/3)/2 of these when n>=12; at n=6 the only pair
would be equal and was already treated separately. This proves (1).
Finally inserting sum_{d|gcd(a,n)}mu(d) for the coprimality indicator and
summing a=d,2d,...,r_d d proves (2). QED.

## Immediate simplifications and queries

Modulo4 the theorem becomes

    U_n = (-1)^(n/2+1) phi(n)/2 + 2[n in {4,6}] (mod4).

Modulo2, U_n is odd exactly for n=4 or n=2 p^k with p an odd prime
congruent to 3 modulo4 and k>=1. This is a first-bit consequence, not a
claim that all higher valuations equal those of phi(n)/2.

If 8 divides n, the complete three-bit formula simplifies to

    U_n = -phi(n)/2 +4D(n) (mod8).

If n=4 modulo8 and n>4, it becomes

    U_n = 3phi(n)/2 +4D(n) (mod8),

because all lower-half representatives are odd. General n=2 modulo4 is
handled by the divisor query, retaining the lower-half residue sidecar.

The provided script tests full Burnside at orders 4,6,...,40. It also queries
n=2^1000, 2*3^1000, 4*5^1000, and 2*1000000007^2, given their factorizations.
These are exact residue calculations at huge indices, not exact full a(n).

## Connection and loss ledger

Source: cycle-index Burnside sum. Target: truncated reduced-residue arithmetic.
Map: discard terms with delta>=3, then use inverse(u)=u modulo8 and Möbius
inversion. Preserved predicate: the first three bits after removing the
forced 2^(n/2) factor. Destroyed: all higher bits, orbit identities, class
multiplicities, and Hamiltonian-path values. Sidecar: factorization plus the
lower-half reduced-residue sum; the exceptional repeated-cycle contributions
must be retained. Cheapest decisive tests: n=4 and n=6 catch omitted symmetry
denominators; n=12 catches the gcd-three correction; n=14 tests the
nontrivial lower-half representative sum.

This is a genuine bridge to elementary multiplicative arithmetic because the
map eliminates an entire partition sum for the requested observable. It does
not supply a rational global generating function: the Möbius divisor sum
retains arithmetic dependence on n. At precision four and beyond additional
cycle and gcd strata can survive. Compiling those strata into divisor queries
is a concrete next target, not an asserted completed recurrence.
