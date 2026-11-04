# Palindromic wings: exact clocks, integral observers and ramified fibres

**PROVED elementary statements below; CITED literature is identified
separately.** No new primality record or infinite family of primes is claimed.

## P1. A clock theorem modulo every modulus

For b>=2 and m>=0 put `P_b(m)=b^(2m+1)-b^m-1`. Its base-b digits are
`(b-1)^m (b-2) (b-1)^m` when m>=1. For b=10 these are the9-wing
palindromes, including the supplied exponent m944264.

For any M>=2, factor M=M0*M1, where all prime factors of M0 divide b and
gcd(M1,b)=1. Then the **minimal scalar** clock of P_b(m) mod M has

    exact preperiod h=max_(p|M,p|b) ceil(v_p(M)/v_p(b)),
    exact eventual period ord_M1(b),                 (P1)

with empty maximum0 and ord_1(b)=1. Two consecutive scalar observations
identify every state of this actual eventually periodic source orbit.

Proof. Write x=b^m. For p|b, `P_b(m)+1=x(bx-1)` has valuation m*v_p(b).
Each transient valuation is distinguishable until saturation, and the M0
output is then permanently -1. On M1, a candidate period a, u=b^a, gives

    P_b(m+a)-P_b(m)=(u-1)b^m [b^(m+1)(u+1)-1].

The gcd of two consecutive bracket values is1: it divides b-1 and the
first bracket is congruent u=1 modulo b-1. Thus a period requires and is
implied by M1|u-1, even after any prefix. For two different phases x=b^i,
y=b^j, equality of two successive outputs gives divisibility of
`(x-y)[b(x+y)-1]` and `(x-y)[b^2(x+y)-1]`. Their brackets have gcd1,
since x+y=2 modulo b-1. Hence M1|x-y. Combine this with the distinguishing
transient valuations and CRT. QED.

One read need not suffice: decimal phases m1,m2 both give2 modulo7, but
their next outputs differ. The result determines the minimal clock on
the source trajectory, not the memory for arbitrary linear combinations.

## P2. The full integer observer has a different size

The sequence has minimal rational recurrence polynomial
`(X-1)(X-b)(X-b^2)`. Use D=diag(1,b,b^2), C=(-1,-1,b), and arbitrary
integer initial weights. The first three observation rows C,CD,CD^2 have
Smith factors

    (1, b-1, b^2*(b-1)^2*(b+1)).                      (P2)

Indeed an entry is a unit, the gcd of all2-minors is b-1, and the absolute
determinant is b^2*(b-1)^3*(b+1). The characteristic polynomial expresses
all later rows as integer combinations of these three. Consequently the
smallest quotient of `(Z/M)^3` for this linear operation/observation has
cardinality `product_i M/gcd(M,s_i)`, with s_i from P2. Distinct roots
1,b,b^2 and nonzero coefficients prove rational minimality.

At b10 the factors are(1,9,89100). At M7 the arbitrary-weight observer
needs343 states, while the actual power orbit has only6 states and is
identified by two reads. This is a precise example of realization guards
reducing memory; it is not permission to discard correlations in other
source laws. Primes dividing b produce nilpotent loss and preperiods;
they cannot be treated as invertible Frobenius transport.

More generally, for `f_m=a*b^(2m)+c*b^m+d` and gcd(M,b)=1, t is a period iff

    M | (b^t-1)*gcd(a*(b^t+1)+c, (b-1)*c).           (P3)

This follows by taking the gcd of all successive brackets
`a*b^m*(b^t+1)+c`; the first and its difference with the second already
generate it. The constant d does not affect the period.

## P3. Quadratic norms and complete congruence charts

For x=b^m,

    P_b(m)+1 = b^m*(b^(m+1)-1),
    (2b*x-1)^2-(4b+1)=4b*P_b(m).                    (P4)

For b10 the discriminant is41, a quadratic nonresidue modulo7,11,13,101.
Thus every decimal wing with m>=1 is coprime to both1001=7*11*13 and
2222=2*11*101; parity supplies the factor2 statement. This is an infinite
family sieve exclusion, not a claim that surviving terms are prime.

Trial division proves7417 prime, and exact modular exponentiation gives
ord_7417(10)=2472. Exhausting this full clock gives

    7417 | P_10(m) iff m=1312 or2395 mod2472.          (P5)

The norm polynomial has at most two roots; the saved residue enumeration
supplies both and proves no phase was missed. The supplied m944264 is
phase2432. The decimal residue clocks at49 and343 have periods42 and294.
These are rings Z/49 and Z/343, not the fields F49 and F343 used for the
tournament observer. Their characteristic and transport operations differ.

For the Belphegor-shaped family put

    B_m=10^(2m+4)+666*10^(m+1)+1,   m>=0.

With x=10^(m+1) its exact norm identity is

    (200x+666)^2-4*7^3*17*19 = 400 B_m.               (P6)

This gives a complete valuation classification:

    v_7(B_m)=0, if m !=4 mod6;
              2, if m=4 mod6 but m !=40 mod42;
              3, if m=40 mod42.                     (P7)

Proof. The linear term vanishes modulo7 precisely when x=5 modulo7,
equivalent to m=4 modulo6. Its unique lift to zero modulo49 is x=5
modulo49, equivalent to m=40 modulo42 because10 has order42 there and
10^-1=5. If its valuation is1, the square has valuation2 while the
discriminant has valuation3. If its valuation is at least2, the square
has valuation at least4 and the difference has valuation exactly3.
Division by400 changes no7-valuation. Thus valuations1 and>=4 never
occur. In the third case `B_m/343=3 mod7`, directly from P6. QED.

Modulo7417 the discriminant443156 is a nonsquare, so no B_m is divisible
by7417. The arbitrary-weight three-mode observer for B has Smith factors
(1,18,296703000), by the same minor calculation as P2. This is separate
from the exact source clock and from primality.

## P4. What the primary sources contribute

**CITED.** Banks, Hart and Sakata, *Almost All Palindromes Are Composite*
(2004), prove asymptotic compositeness and equidistribution of the full
palindrome population under specified modulus conditions. Their prime
condition includes ord_p(b)>=3*sqrt(p); p7417,b10 satisfies it. These
results concern all palindromes, not the thin one-parameter wing family.
[Primary manuscript](https://arxiv.org/abs/math/0405056).

**CITED.** The supplied Propper--Batalov number is recorded as proven by
PrimePages, submitted in2021. Its recorded N+1 test uses a quadratic
extension with discriminant7 and a partially factored N+1. P4 displays
why this family's N+1 has accessible structure; it is not a reconstruction
of that large certificate. The polynomial discriminant41 in P4 is a
different object from the certificate's extension discriminant7.
[Primary record and verification data](https://t5k.org/primes/page.php?id=132851).

The operation transferred to tournament research is the distinction
between a small observer on a realized orbit, its integral completion
on arbitrary weights, and transport that ceases to be invertible at
exceptional primes. No map from prime density to the Hamilton spectrum
has been proved.
