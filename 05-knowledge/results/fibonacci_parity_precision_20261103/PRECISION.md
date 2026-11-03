# Three normalized binary digits of A000568

**PROVED, all orders, from the inherited Davis formula and its cycle-budget
bound.** Here A(n)=A000568(n), B=floor(n/2), and u(n)=A(n)/2^B for n>=3.
At orders0,1,2 use the unnormalized value1. The divisibility by2^B is
inherited from
[fixed precision P3 through P7](../port_residual_integral_20261025/precision/README.md).
That packet already proves u modulo4 and an arbitrary fixed-precision
algorithm; the extension here is its explicit modulo8 evaluation.

## P1 The formula

Write phi for Euler's totient and rad for the product of distinct prime
divisors. For even n>=6 put

    U(n) = (n-1)phi(n)/2                         if 4 divides n;
    U(n) = m * product_(p|m)(1-p)/2, m=n/2       if n=2 mod4.

The product is mu(rad(m))*phi(rad(m)); m is odd and at least3, so U is an
integer. Then

    u(n) = U(n) + 2*[3|n]*phi(n/3)     (mod8).             (P1)

For odd n>=5, let C(n) count unordered, DISTINCT, positive odd triples
a<b<c with a+b+c=n and pairwise gcds all equal to1. Then

    u(n) = inverse(n) + 2*inverse(n-2) + 4*C(n) (mod8).   (P2)

Inverses in(P2) are modulo8. The remaining cases are u(3)=u(4)=1.
These are formulas for class counts, not Hamilton-path counts of one
tournament. C(n) retains a joint arithmetic condition, not just the number
of odd triples or the condition gcd(a,b,c)=1.

## P2 Which permutation cycles survive

For an odd partition lambda=(d_1,...,d_l) of n, put

    e=(n-l)/2+sum_(i<j)gcd(d_i,d_j),
    z=product_d d^(m_d)*m_d!.

Davis's formula is A(n)=sum_lambda 2^e/z. It is reproduced in
[OEIS A000568](https://oeis.org/A000568), with the Davis paper reference.
The inherited precision proof bounds each term separately after removing
the powers of2 in z; no unjustified cancellation or division by an even
modular denominator is used.

For even n>=6 only l=2 contributes modulo8. At l=4, multiplicity patterns
3+1 and2+1+1 have normalized valuation at least3, the distinct and2+2
patterns at least4, and four equal parts d have valuation6d-5. Only d=1
is exceptional, at n=4. At l>=6 the bound is at least8.

For odd n only l=1,3 matters modulo8. For three distinct parts the
normalized term has valuation G-1, G=sum pairwise gcds. Since every gcd is
odd, valuation2 occurs exactly when all three gcds are1; each such triple
contributes4. For (a,a,b), a!=b, the valuation is a+2gcd(a,b)-2; it is below3
exactly when a=1. This contributes2/(n-2). For three equal parts a, the
valuation is3a-2, so only n=3 is exceptional. The one-part contribution is
1/n. This proves(P2), including the separately stated small orders.

## P3 The even formula by a Mobius sum

For distinct odd parts a<b, a+b=n, the normalized contribution is

    2^(gcd(a,b)-1)/(ab).

Only gcd1 or3 survives modulo8. If gcd1, a runs over the units modulo n
below n/2. There are phi(n)/2 such representatives. Every odd unit is its
own inverse modulo8 and a^2=1 modulo8, hence their sum is

    n * sum_(0<a<n/2, gcd(a,n)=1) a - phi(n)/2  (mod8).  (P3)

If4|n, every summand a is odd; multiplying a-1 by n kills it modulo8.
This gives U=(n-1)phi(n)/2.

If n=2m with odd m>1, write R for that half-unit sum. Inclusion-exclusion
over the prime divisors of m gives an exact integer identity:

    R = sum_(d|m) mu(d)*d*((m/d-1)/2)^2
      = [m*phi(m)+product_(p|m)(1-p)]/4.                (P4)

The first equality sums odd positive multiples d*j<m: their j's are
1,3,...,m/d-2, with sum((m/d-1)/2)^2. Expanding the square uses
sum mu(d)/d=phi(m)/m, sum mu(d)=0 and sum mu(d)*d=product(1-p).
Substitute(P4) into(P3). The result is

    [(m^2-1)phi(m)+m*product_(p|m)(1-p)]/2.

The first term is divisible by8 because m is odd and phi(m) is even.
The second is U. This is the actual Mobius mechanism: it removes forbidden
prime divisors from a half-interval sum of units. No Liouville substitution
would preserve that inclusion-exclusion predicate.

For n>6 divisible by3, distinct pairs with gcd3 correspond to the
phi(n/3)/2 coprime unordered odd pairs summing to n/3. Each contributes4
modulo8, giving2phi(n/3). At n=6 the distinct bank is empty, but the equal
pair(3,3) contributes2/9=2 modulo8, the same correction2phi(2).
Equal pairs(a,a) with a>=5 have normalized valuation a-2>=3 and disappear.
This proves(P1).

## P4 An exact valuation region and its first failed extrapolation

For even n>=4 such that t=v2(phi(n))<=3,

    v2(A(n)) = n/2 + t - 1.                            (P5)

For4|n, U has valuation t-1. For n=2m, the product formula for U has the
same valuation, since m and all prime-power factors in phi(m) are odd.
When3|n and n>6, the correction2phi(n/3) has valuation at least t: removing
one factor3 changes the totient by a factor3 or2. The case n=6 is direct.
Thus the correction cannot cancel U's leading bit. Modulo8 sees that bit
exactly when t-1<3. Order4 is separately direct.

For odd n, u(n) is odd by the inherited theorem, so v2(A(n))=(n-1)/2.

**REFUTED extrapolation:** dropping t<=3 from(P5). The first failure is
n=40: v2(A(40))=26, while20+v2(phi(40))-1=23. This is an exact finite
minimum, checked at every smaller even order, not an asymptotic statement.
The mechanism is visible at normalized precision7:

    contributions from two permutation cycles:   8 mod128;
    contributions from four permutation cycles: 56 mod128;
    all other lengths:                           0 mod128;
    total:                                      64 mod128.

Thus distinct cycle-length layers cancel three candidate leading bits.
The strongest surviving statement is(P5). The missing coordinate outside
it is the higher-cycle contribution, not an extra ordinary class pairing.
The next question is whether these layers admit similarly short formulas
modulo16 and32; the inherited fixed-cycle algorithm already computes them.

## P5 Computational meaning

For even n, after obtaining the prime factorization of n, formula(P1)
uses only products over its distinct primes. It avoids the full class
integer and even the bounded four-cycle partition sum for this precision.
This is a fixed-precision arithmetic improvement, not a faster complete
A000568 enumeration or a structural classification of all tournaments.

For odd n, the present direct evaluator counts C(n) in O(n^2) gcd tests.
That part is not claimed to be an optimal algorithm. The appearance of
three jointly coprime cycle lengths explains why one should retain joint
incidence when moving beyond a leading bit.
