# The next two binary digits of A000568

**PROVED F1–F4 by termwise Burnside valuation and CRT. FINITE-EXACT through200.**
The modulo4 laws and bounded-cycle algorithm are inherited from
[October25 P8–P11](../port_residual_integral_20261025/precision/README.md).
This note computes modulo16 explicitly. No full-value enumeration speedup
is claimed.

Incoming commit947b6c3852, read during final integration, independently
published the [modulo8 formulas](../fibonacci_parity_precision_20261103/PRECISION.md)
and the n40 valuation hostile before this packet's first push. It has
priority here for that third-bit result and also gives a shorter factored
even-order expression modulo8. Formula(6) agrees with it. The present
extension is modulo16, including the four-cycle totient(n-2) contribution
and the weighted CRT evaluation below.

Write b_n=a_n/2^floor(n/2) for n>=3. All inverses below are taken in
Z/16, or the explicitly indicated lower-power quotient, and have odd
denominators. Let

    S_M(q)=sum_(1<=a<q/2, a odd, gcd(a,q)=1) (a(q-a))^-1 mod M

for even q>=4. Define T_n to be the distinct odd positive triples
a<b<c with a+b+c=n and all three pairwise gcds equal1.

## F1. Odd order: a triple interaction term

For every odd n>=5,

    b_n = n^-1 + 2(n-2)^-1
          +4 sum_((a,b,c) in T_n) (abc)^-1
          +8 [n>=7 and 3 does not divide n]                 mod16. (5)

The boundary is b_3=1.

**Proof.** For an odd partition lambda of n with k parts, put
G=sum_(i<j)gcd(lambda_i,lambda_j), t=sum_d v2(m_d!).
Its normalized term has valuation

    (1-k)/2+G-t.

The inherited lower bound (k-1)^2/2-v2(k!) is at least5 for odd k>=5;
all such terms vanish modulo16. For k=1 the contribution is1/n.

For k=3 with distinct parts, t=0 and the relative valuation is G-1.
All gcds are odd. Thus it is2 exactly when all are1; otherwise it is at
least4. The former terms are4/(abc), exactly the displayed sum.

For two equal parts a,a and a distinct c, t=1 and relative valuation
is a+2gcd(a,c)-2. Below4 this allows only a=1 (giving2/c),
or a=3 with gcd(3,c)=1 (giving8/(9c)=8 mod16). The latter exists
exactly when n=6+c>=7 and3 does not divide n. For three equal parts a,
the relative valuation is3a-2; only a=1 survives, at n=3.
The direct Burnside value a_3=2 gives the stated boundary. QED.

Modulo8 the formula simplifies for n>=5 to

    b_n = 3n-4+4|T_n|                              mod8.     (6)

Modulo4 it recovers the inherited b_n=n+2. Distinct coprime triple
interaction is the first extra coordinate beyond the reciprocal correction.
It is not the number of directed triangles in a tournament.

## F2. Even order: a totient shift plus a threefold divisor correction

For every even n>=4,

    b_n = S_16(n)
          +4 [n>=12 and 6 divides n] S_4(n/3)
          +8 [n>=6] phi(n-2)/2
          +6[n=4]+2[n=6]+8[n=10]                   mod16.    (7)

The guarded terms are evaluated only when their indicated domains hold.

**Proof.** For even k the relative valuation is G-k/2-t.
The inherited lower bound k(k-2)/2-v2(k!) is at least8 at k=6
and increases thereafter. Thus only k=2,4 matter.

For distinct two-part partitions(a,b), the relative valuation is
gcd(a,b)-1. Gcd1 gives S_16(n). Gcd3 contributes
4/(9uv)=4/(uv) mod16, where a=3u,b=3v and u,v are distinct coprime
odd parts of n/3; this is the second term. All larger gcds vanish.
For equal parts(d,d), relative valuation is d-2. The only surviving
d>=3 are3 and5, contributing2/9=2 at n=6 and8/25=8 at n=10.

For four parts, inspect the multiplicity pattern:

- 4: relative valuation6d-5; only four ones survive, contributing2/3=6
  at n=4.
- 3+1: relative valuation3a+3gcd(a,b)-3. Below4 only
  (1,1,1,n-3) survives, for n>=6, contributing8.
- 2+1+1: relative valuation
  a+2gcd(a,b)+2gcd(a,c)+gcd(b,c)-3. Below4 precisely
  a=1 and b<c odd>=3 coprime survive, contributing8 each.
- 2+2: the valuation a+b+4gcd(a,b)-4 is at least4.
- 1+1+1+1: the valuation is at least4.

At n>=6 the two surviving nonconstant patterns correspond exactly
to all unordered coprime odd pairs summing to n-2: the pair(1,n-3)
represents the first pattern, and the pairs with both entries>=3 the
second. There are phi(n-2)/2 such pairs. This gives the third term and
all explicit corrections. QED.

## F3. What order +2 and count doubling actually share

For n=2r+1>=3, the inherited modulo4 consequence is

    a_(n+2) = 2a_n+2^(r+2)                   mod 2^(r+3).    (8)

Indeed b_(n+2)-b_n=2 mod4 and
a_(n+2)-2a_n=2^(r+1)(b_(n+2)-b_n).
Let Q(n) be the right side of(5), with Q(3)=1. The stronger exact
low-precision recurrence is

    a_(n+2) = 2a_n+2^(r+1)(Q(n+2)-Q(n))      mod 2^(r+5).    (9)

This is a uniform recurrence for four normalized binary digits. It does
not assert a_(n+2)=2a_n, construct a tournament, or control higher digits.

## F4. Weighted CRT evaluates the pair sums without a linear scan

Let q>=4 be even, M a power of2, L=q/2-1, and R the product of the
distinct odd primes dividing q. For each odd r in1,...,M-1 and divisor
d of R, let t(r,d) be the least positive residue of r*d^-1 modulo M.
Then the exact number of admissible a in residue r is

    N_r(q)=sum_(d|R) mu(d) max(0,1+floor((floor(L/d)-t(r,d))/M)). (13)

Consequently

    S_M(q)=sum_(r odd mod M) N_r(q)*(r(q-r))^-1 mod M.        (14)

Proof: inclusion-exclusion removes multiples of every odd prime dividing
q. For each remaining term write a=d*t. The congruence a=r mod M is
equivalent to t=t(r,d) mod M because d is odd; count this progression up
to floor(L/d). The weight (a(q-a))^-1 depends only on r. Every denominator
is odd. Summing these exact counts proves both equations.

This is the same elementary prime-wheel inclusion-exclusion used in the
[preceding arithmetic wheel packet](../hourglass_triangular_dynamics_20261103/PRIME_WHEELS.md),
now with inverse-product weights and a half-interval boundary. Repeated
prime powers disappear from R but **remain in q and L**. Replacing q by
its radical would change the clock and give the wrong counts.

For M=16 the weight itself needs a smaller, provably minimal cyclic clock:
on odd a under the step a->a+2, its period is **four** if4 divides q and
**two** if q=2 mod4. Equivalently the integer-a phases are modulo8 and4.
Indeed x(a)=a(q-a) satisfies
x(a+h)-x(a)=h(q-2a)-h^2. At h=8 this vanishes modulo16 for every even q.
For4|q, h=4 gives8 modulo16 and fails. For q=2 mod4, h=4 vanishes, while
h=2 gives4 times an odd integer and fails. Inversion is a permutation of
the odd units, so it preserves these exact periods. The implementation
uses this smaller phase in the progression count but still inverts the
weight modulo16. This is the smallest clock for the **weight under the
specified step**, not for the entire coprimality predicate; the divisor
filters and interval bound are separate retained data.

Given the prime factorizations, (14) uses O(M*2^omega(R)) operations on
O(log q)-bit integers. Evaluating(7) needs the factors of n and n-2;
the n/3 term uses factors already present in n. Our implementation uses
trial division and charges its cost; no polynomial-in-log(n) factoring
claim is made. The new even-order formula replaces the inherited generic
four-cycle partition scan, which visits Theta(n^3) tuples, by an O(n)
gcd scan or this factored CRT sum. It returns only four normalized bits.

Direct weighted residue scans independently check(13),(14) for every
even q=4,...,1000, each M=2,4,8,16,32, and the additional q values
338,5510,20402,30030,1000000,1000002. These include odd prime powers
and the user-relevant modulus2*101^2. Normal and optimized receipts agree.
The standalone benchmark includes factorization and checks all three
methods agree at n=200,400,800. At800 the generic method visited454711
partitions; exact timings for all three implementations are recorded in
benchmark.json. The tiny CRT timings are illustrative, not stable
hardware-independent ratios; the operation-count reduction is the rigorous
claim. Large exact formula evaluations, not
independent full counts, are separately labelled in benchmark.json.

## Hostiles and costs

A tempting extrapolation from the first bit is
v2(a_n)=n/2+v2(phi(n))-1 for every even n. It is false.
The first failure among n=4,6,...,76 is n=40:
v2(a_40)=26 rather than23. At48 it is30 rather than27; at64 it is35
rather than36. These are finite exact rejection witnesses, not a new
global classification of the failures.

The odd formula costs O(n^2) direct triple visits; the even formula and its
CRT improvement have the costs just specified. Arithmetic is at fixed small
precision. The improvement does not reconstruct the remaining higher bits
in the full count or enumerate individual isomorphism classes.

Independent full Fraction sums check(5),(7) through40. The inherited
bounded-cycle evaluator checks them through200, including all small
exceptions. A guess matching four digits still does not certify higher
digits or the complete count.

