# Rational anchors that certify Collatz returns

2026-10-04. **PROVED** elementary compiler and all-parameter families below;
**FINITE-EXACT** for the stated bank comparisons and replays. **OPEN** for
coverage of every positive integer or every negative basin. No novelty claim
is made for rational Collatz cycles, affine words, or inverse iteration.

The useful addition is a compiler whose repeated anchor can be a negative
**rational** cycle. The example -19/11 gives new first-descent cylinders,
separate from the explicitly named existing bank and return families. The
same chart produces compact completed certificates into each of the four
known nonzero signed cycles. These are families of certified inputs, not a
universal certificate construction.

## Inheritance and conventions

Use ordinary C(n)=3n+1 for odd n and n/2 for even n. On nonzero odd integers
write U(n)=(3n+1)/2^v2(3n+1). On negative integers, v2 means the valuation of
the absolute numerator; the division retains its sign.

The [two-copies investigation](fibonacci_two_copies_collatz_20260926.md)
already exhibits the rational growth shadow -19/11, its word (1,1,2), and
multiplier 27/16. The [shadow-error investigation](collatz_shadow_error_flp_20260927.md)
also uses this anchor. Its appearance is inherited. The new assertion here
is the source-decodable exit condition and its additional certified coverage.

Closest proved mechanism: [the recursive entry grammar](entry_20260927_recursive.md),
[the -17 return](collatz_minus17_return_20261003.md), and the incoming
[head/cycle compiler](collatz_mixed_return_compiler_20261003.md) in commit
75daccc120. The latter requires an actual negative integer cycle; the
present compiler admits rational cycles with odd denominator. The canonical
hostile is a legal long growth shadow whose chosen exit is still above its
original source. The corrected near miss is the final oddness bit in
[THM-4512, coefficient descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
The underused sidecar is the rational denominator of the anchor.

The working board is rational anchors, exact valuation guards, retained
source bounds, golden word polynomials, and completed terminal certificates.
The assigned anchor is Collatz certification; the niche is nonintegral
periodic shadows; the wildcard is symbolic powers as compressed numerals.

## A compiler for every expanding rational chart

Let w=(a_1,...,a_p), with each a_i>=1, and define

    A_i=a_1+...+a_i, A=A_p, P=3^p, Q=2^A,
    B=sum_(i=0,...,p-1) 3^(p-1-i)*2^A_i, A_0=0.

The formal composition is W(n)=(P*n+B)/Q. Require every nonempty prefix
to have expanding linear coefficient:

    3^i > 2^A_i, 1<=i<=p.                            (1)

Its fixed point is r=-h/d=-B/(P-Q), reduced with h,d positive. Both h and d
are odd, and neither is divisible by 3. In particular, d is invertible
modulo every power of 2 and 3. Indeed B is odd and is a power of 2 modulo
3, while P-Q is odd and nonzero modulo 3.

The formal rational orbit of r really has word w in the odd-denominator
extension of U. Work backwards from r_p=r_0 in the ring of rationals with
odd denominator. Each inverse equation
r_i=(2^a_i*r_(i+1)-1)/3 preserves that ring and gives an odd numerator.
The forward formula also puts r_i over d*2^A_i; canceling its powers of two
therefore puts every phase over the common odd denominator d. Each equation
3*r_i+1=2^a_i*r_(i+1) has exactly the stated valuation. All r_i
are negative: r_0 is negative and backward propagation from r_p=r_0 through
r_i=(2^a_i*r_(i+1)-1)/3 preserves negativity.

Choose any m>=1 with Q^m>h. Let t=t_m>=1 be the least integer satisfying

    2^t*(Q^m-h) > P^m-h.                              (2)

Set beta_m=h*(P^m)^(-1) mod2^t, taking its odd representative in
1,...,2^t-1. Define K=Am+t and the complete positive cylinder

    n = (beta_m*Q^m-h)*d^(-1) mod2^K.                 (3)

**Compiler theorem.** Every positive integer in (3) has its exact first
odd descent at step pm. With

    b=(d*n+h)/Q^m, z=b*P^m-h, tau=v2(z),

its word through that step is

    w^(m-1), a_1,...,a_(p-1), a_p+tau,               (4)

and its endpoint is

    U^(pm)(n)=z/(d*2^tau),  0<U^(pm)(n)<n.           (5)

**Legality.** Congruence (3) is exactly b=beta_m mod2^t. Thus b is a
positive odd integer, and tau>=t. The initial displacement from r is

    n-r=b*Q^m/d.

At a proper prefix of the m repeated blocks, the corresponding displacement
from its rational phase r_i is

    b*3^j*2^(Am-A_j)/d,

where A_j<Am. It is even in the odd-denominator ring. Consequently the
formal state remains odd. Because n is integral and d is odd, at each
intermediate division this exact rational parity also forces ordinary
integer divisibility. All earlier valuations are exact. The final formal
state is z/d, which has valuation tau; the last step absorbs those extra
divisions and yields (4), (5). Positivity also follows directly by iterating
positive affine maps from n>0.

**First descent.** Each proper prefix of w^m has multiplier >1: write it
as full blocks followed by a prefix and use (1). Its carry is positive, so
its value exceeds n. The exit budget pays the original source because

    d*(2^t*n-z)
       = b*(2^t*Q^m-P^m)-h*(2^t-1) > 0.

Here (2) says that the coefficient of b exceeds h*(2^t-1), and b>=1.
This proves strict descent without replacing the original source by a
later, larger value.

For different m the cylinders are disjoint since

    v2(d*n+h)=Am.                                    (6)

This is also a source decoder: compute that valuation, recover m, check
(2), (3), and compute the endpoint. No successful forward search is assumed.

## The family anchored at -19/11

Take w=(1,1,2). Then

    W(n)=(27n+19)/16, r=-19/11,
    r -> -23/11 -> -29/11 -> r.

For m>=2 let t_m be least with

    2^t*(16^m-19)>27^m-19.                           (7)

The first cylinders and their least positive members are:

| m | t_m | Source cylinder | Endpoint of the least source | Exact first descent |
|---:|---:|---|---:|---:|
| 2 | 2 | 999 mod 1024 | 89 | 6 |
| 3 | 3 | 21223 mod 32768 | 12749 | 9 |
| 4 | 4 | 303847 mod 1048576 | 153997 | 12 |
| 5 | 4 | 3908327 mod 16777216 | 3342643 | 15 |

Every source is 7 mod16. For m>=3 all sources lie in the single parent
class n=743 mod2048, since 11*n+19=0 modulo 2^(4m). Exact comparison against
all 171 rows of [the old table](reset_20260926_swaplift.out), whose union
normalizes to 65 cylinders, proves that the entire parent class is missed.
The m=2 cylinder does overlap that bank and is excluded from the added mass.

This is also disjoint, for all heights, from these named infinite families:

- The inherited -5 return sources are 3 mod8.
- The pure -17 sources are 15 mod16.
- The incoming mixed (1,2) head into -17 has sources 3 mod8.
- The incoming (1) head into -17 has sources 15 mod16: its decoder is
  v2(3n+35)=11m+1, so 3n+35=0 mod16.

The assertion is relative to these explicit banks and families. It is not
a claim to differ from every earlier inverse family or every possible head.

The added natural density among all positive integers is

    delta=sum_(m>=3) 2^(-4m-t_m)
         =0.000031532779875335455... .                (8)

This countable union has a natural density: its tail m>M is contained in
the single cylinder 11n+19=0 mod2^(4(M+1)), whose density tends to zero.
Also (7) implies 2^t>(27/16)^m, giving the strict numerical tail bound

    sum_(m>M) 2^(-4m-t_m) < 1/(26*27^M).              (9)

The [output](rational_anchor_returns_20261004.out) retains the exact m=3..30
partial sum and its tail bound. These are densities of starting integers,
not orbit frequencies.

**Hostile and repaired boundary.** With m=2, weakening t=2 to t=1 admits

    487 -> 731 -> 1097 -> 823 -> 1235 -> 1853 -> 695.

All six steps stay above 487. The first failed implication is that a legal
exit is a descent. The survivor is the exact word; the repair is (7), which
retains the source comparison. This does not assert that 487 never descends.

## Completed certificates into all four signed roots

The same chart gives more than partial descents. Let rho be any odd integer
not divisible by 3. For any m>=1, choose t>=1 such that

    11*rho*2^t = -19 mod27^m.                         (10)

Then define

    n(m,t,rho)
      =[16^m*(19+11*rho*2^t)-19*27^m]/[11*27^m].    (11)

For all sufficiently large t satisfying (10), n has the sign of rho and
the exact route

    U^(3m)(n)=rho,
    word=(1,1,2)^(m-1),(1,1,2+t).                    (12)

There are infinitely many such t for each fixed m and rho. Indeed 2 is a
generator of the units modulo 3^k: 2 has order 2 modulo 3, and
v3(2^(2*3^j)-1)=j+1 by the binomial expansion of 4=1+3. Thus (10) selects
one residue class modulo 2*3^(3m-1). It can be found by lifting one ternary
digit at a time, testing only three candidates at each lift.

**Integrality and exactness.** Put b=(19+11*rho*2^t)/27^m. Equation (10)
makes b an odd integer. Since 27=16 mod11, b*16^m=19 mod11, proving the
integrality of n=(b*16^m-19)/11. The proper-prefix displacement argument
above does not require b positive and proves the exact word through its
last nominal division. At that division the numerator becomes 11*rho*2^t,
so its additional valuation is exactly t. This gives (12). Taking t large
ensures the specified sign. Since positive and negative integers are each
invariant under U, all intermediate integers then have that sign as well.

For rho in {1,-1,-5,-17}, the terminal is on the corresponding known signed
cycle. In particular the certificates reach the designated root without
requiring a conjectural suffix. An earlier root visit is allowed by the
generic certificate type; a first-hit graph representation should truncate
there if needed.

Examples for rho=1:

| m | Least positive solution t | Period of the allowed t | Source size | Odd steps to 1 |
|---:|---:|---:|---|---:|
| 1 | 8 | 18 | n=151 | 3 |
| 2 | 278 | 486 | 277 binary digits | 6 |
| 3 | 11942 | 13122 | 11940 binary digits | 9 |
| 4 | 11942 | 354294 | 11939 binary digits | 12 |

The first route is 151 -> 227 -> 341 -> 1. The next two exponent classes
are t=2137706 mod9565938 and t=97797086 mod258280326; the script verifies
their guards without expanding those integers. Exponent size and number of
odd steps are different resources; the final division can be enormous.

The [JSON examples](rational_anchor_returns_20261004.json) contain twelve
fully replayed signed certificates, m=1,2,3 for each of the four roots.
Their source is stored as the power expression (11), with chart, repeated
block count, terminal exponent, and cycle-root pointer. This makes the
certificate part of the numeral's syntax. It does not make every arbitrary
input an instance of that syntax.

## Golden synchronization and the proposed family carrier

For the raw parity polynomial of one (1,1,2) block,

    H_w(z)=1+z^2+z^4, L_w=7.

If the repeated chart exits after tau additional halvings, its finite word
polynomial and length are

    H_m(z)=(1+z^2+z^4)*(1+z^7+...+z^(7(m-1))),
    L=7m+tau.

Consequently its exact parity fraction satisfies

    Theta(n)=beta*H_m(beta)+beta^L*Theta(endpoint),
    beta=phi^(-1).                                   (13)

This is the incoming compiler's position-retaining concatenation law.
The endpoint can have an open suffix, or a sealed known-cycle suffix. No
claim identifies literal base-phi digits of n with its parity itinerary.

A useful new carrier therefore stores

    (difference representative, rational chart, phase, m,
     dyadic and ternary guards, exit budget, original source,
     endpoint pointer, timed golden polynomial, proof status).

The integer projection is given by the guarded affine expression, independently
of whether a suffix has been sealed. Every integer also has a trivial open
representation, so the carrier does not discard difficult inputs by definition.
Rules for descent can link to a smaller source; rules (10)-(12) link directly
to a known terminal. Repeated blocks are explicit parameters, not fixed-size
colours. This supplies a sound, extensible certificate language.

The outstanding theorem is a total normalizer from arbitrary integers into
sealed expressions. On positives, a proof that every n>1 receives a proved
lower-source link would suffice by induction. On negatives, a certificate
must name and reach its cycle, and universal coverage by the known three
negative cycles remains OPEN. Extra coordinates or a family grammar alone
do not establish either statement.

## Verification

From the repository root:

    python3 04-computation/experiments/rational_anchor_returns_20261004.py \
      --certificates 05-knowledge/results/rational_anchor_returns_20261004.json \
      > 05-knowledge/results/rational_anchor_returns_20261004.out

The complete small-word universe is all primitive growth necklaces with
p<=6 and P>Q, each rotated to its lexicographically first all-prefix-growth
form. There are 23. The script tests twenty consecutive admissible m values
per word and four coefficient lifts (0,1,7,29), for 1840 exact first-descent
controls. A rotation with all positive prefix growth exists by starting
after a minimum of the cyclic logarithmic partial sums; the enumerator
tests the inequalities with integer powers.

For the retained word it tests m=2..30, all 171 bank rows and their 65-class
normalization, the weak-budget hostile, modular terminal guards through m=6,
four completed positive replays and twelve completed signed replays.
Direct odd iteration is checked against independent affine composition.
All checks remain active under optimized Python. The infinite claims rely
on the proofs above; finite replays provide positive and hostile controls.

The next productive test is to rank other rational anchors by added coverage
and by how often their endpoint admits an already-proved suffix. A larger
list of legal words alone does not close the coverage obligation.
