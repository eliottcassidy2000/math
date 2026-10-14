# Predecessor interfaces, Carmichael blocks, and simultaneous obligations

2026-10-14. **PROVED scoped algebra and terminating classifiers;
FINITE-EXACT complete one/two-prime catalogs for three specified roots;
CITED AGP asymptotic theorem. Universal Collatz and the other proof frontiers
remain OPEN.** The derivations below do not claim historical novelty.

The useful consequence of viewing an integer through its neighbours is a
certificate interface. For a Carmichael number N retain

    (N, lambda, e),  lambda=lcm_{p|N}(p-1),  e=(N-1)/lambda,

together with the distinct prime factors. This is a lossless description
when N is retained. The smaller pair (lambda,e) also determines N, but an
efficient verifier still needs its factor certificate. Neither a prime
support nor an isolated residue retains all these obligations.

## Inheritance and the live concept board

The closest mechanisms are the exact prime-exponent divisor balance in
[the divisor family](collatz_mod6_20260917_divisor_balance_family.md), the
primitive role/factor-choice compiler in
[THM-3334, Berggren/Gaussian torsor](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
and [actual-source Collatz gates](collatz_native_tail_recursion_20261013c.md).
The prior [Bernoulli boundary](bernoulli_boundary_20260925.md) already
separates Carmichael and Giuga conditions:561 passes one and fails the other.
That is the canonical hostile to merging adjacent divisibility predicates.
The corrected near miss is a quotient identity promoted to payment or a
terminal. The least-used sidecar here is the finite-group selection argument
in AGP, whose conductor makes its congruence condition sufficient.

The board is: predecessor periods; prime labels; collective blocks; native
source guards; finite-difference row defects; and the distinction between
counting many witnesses and covering every source. The anchor is the last
two Collatz obligations, the niche is Carmichael block classification, and
the wildcard is the Pythagorean/pyramid role change.

## 1 The ring meaning of Korselt's criterion

For n>=2, the function identity

    a^n=a mod n for every integer a

holds exactly when n is squarefree and p-1 divides n-1 for every p|n.
Add compositeness to obtain the Carmichael numbers.

Necessity of squarefreeness is visible without statistics. If p²|n, take
a=p: p^n-p has p-adic valuation1. For each prime divisor p, a generator of
F_p^* then forces p-1|n-1. Conversely these period conditions give the
identity in each F_p, including a=0; the squarefree Chinese remainder
decomposition gives it modulo n. Thus the two pieces have different jobs:
squarefreeness removes nilpotent directions, and predecessor divisibility
synchronizes the multiplicative periods of the field factors.

This is an equality of polynomial FUNCTIONS, not the assertion that the
formal polynomial X^n-X is zero. In characteristic p, X^p-X vanishes on
F_p despite formal derivative -1. Therefore this finite-field phenomenon
cannot be transported to the characteristic-zero planar Jacobian problem
by retaining only a derivative or a function table.

For the user's two supports the transition is exact:

    561=3*11*17;  560=2^4*5*7;
    lambda(561)=lcm(2,10,16)=80;  e=7.

The7 measures slack in predecessor divisibility here. It is not a prime
factor of561. More generally each prime ell has slack valuation

    v_ell(e)=v_ell(N-1)-max_{p|N} v_ell(p-1).

This keeps both supports and their exponents. The accompanying
[prime-shape/geometry packet](predecessor_prime_shapes_20261014b.md) records
the exact F=U+S definitions, the failure at p^4qr, the four-free repair, and
why561 is a valid primitive odd leg but not a primitive hypotenuse.

## 2 Complete one-prime extensions and the interface update

Let N be Carmichael and q a prime not dividing N. Then

    Nq is Carmichael
    iff q=1+d*lambda for some positive divisor d of e.

Indeed the old prime conditions are equivalent to q=1 mod lambda, since
N=1 mod lambda. The new condition q-1|Nq-1 is equivalent to q-1|N-1.
Together these say d|e. Squarefreeness is precisely the distinct-prime guard.
Every accepted edge has the exact interface update

    (N,lambda,e) -> (Nq,d*lambda,N+e/d).

Here lambda(Nq)=lcm(lambda,q-1)=d*lambda. This is a genuine directed
transition system, with a complete finite list of children at each vertex.
It is not a rooted tree covering all Carmichael numbers.

For561, e=7 gives only q=81 or561, both composite: no children.
For1105, (lambda,e)=(48,23) likewise gives only49 or1105.
For1729, (lambda,e)=(36,48) gives exactly

    q=37,73,109,433,577.

If e is prime, the only candidate other than the composite N is
lambda+1. This is a convenient stopping certificate, subject also to the
requirement that lambda+1 not already divide N.

## 3 Complete two-prime extensions: escaping a single-step obstruction

Let q<r be new primes. Put E=N-1, a=q-1, b=r-1, and lambda=lambda(N).
Korselt is equivalent to

    qr=1 mod lambda,
    a | Nb+E,       b | Na+E,

plus the distinct-prime guards. Define t=(Na+E)/b. Since b>a>=2,

    1 <= t <= N+E/2,       a | E(N+t).

The divisibility follows by multiplying Nb+E by t and using
tb=Na+E. Consequently a COMPLETE finite classifier is:

1. Enumerate the bounded integer t.
2. Enumerate positive divisors a of E(N+t), retaining a>=2.
3. Recover b=(Na+E)/t, retaining integral b>a.
4. Check the displayed divisibilities, primality, distinctness, and old
   period condition.

Every actual extension appears, and every retained output passes Korselt.
No cutoff on q or r is assumed. The proof of this finite universe makes
the following exact catalogs stronger than a search up to chosen primes:

| N | one-prime children | two-prime children |
|---:|---:|---:|
|561|0|1|
|1105|0|13|
|1729|5|21|

In particular the unique two-prime extension of561 is

    561*401*3041=684106401.

Both intermediate products561*401 and561*3041 fail Korselt. The new
primes discharge each other's obligations:

    561*401-1=74*(3041-1),
    561*3041-1=4265*(401-1).

The output interface is (684106401,15200,45007), with

    15200=2^5*5^2*19,
    684106400=15200*45007.

45007 is prime, while15201=27*563 is composite. Thus the single-prime
obstruction repeats at the child. This proves a repeated obstruction and
one successful simultaneous repair; it does not prove indefinitely many
such repairs or closure under arbitrary blocks.

## 4 A finite product-one compiler and the exact CRT guard

Choose integers L>=2,k>=1 with gcd(L,k)=1. Choose distinct primes of the
form p=dk+1, d|L, excluding p|L. If a subset S has at least two elements
and product N=1 mod L, then N is Carmichael: each p=1 mod k, hence
N=1 mod kL by CRT, while p-1|kL. These are explicit sufficient conditions,
and permutation of the selected primes preserves the certificate.

Example: L=80,k=3 has bank7,13,31,61,241. The subset7,13,31,61 gives
N=172081, with N=1 mod240. This finite example does not impose the
additional analytic hypotheses used in AGP to manufacture large banks.

The coprimality guard cannot be omitted. With L=6,k=2, primes7 and13
give91=1 mod6 and mod2, but91 is not1 mod12;13-1 does not divide90.
So even a seemingly harmless loss of conductor information destroys the
compiler. This is the exact analogue of requiring the native cell and
integer carry before accepting a Collatz block.

**CITED:** [Alford–Granville–Pomerance, *There are infinitely many
Carmichael numbers*](https://math.dartmouth.edu/~carlp/PDF/paper95.pdf),
pp.705–707 and717–719, builds many such primes using smooth shifted primes
and primes in progressions, with coprime k. It keeps the relevant group
isomorphic to (Z/LZ)^*, controlling its exponent despite introducing k.
Product-one selection and counting then yield

    C(x)>=x^(beta-epsilon) eventually for every epsilon>0,
    beta=(5/12)*(1-1/(2*sqrt(e)))=0.290306...>2/7.

Here the e in beta is Euler's number. Thus C(x)/x^(2/7) tends to infinity.
The exponent combines two analytic
inputs; it is not a consequence of the parity of7 or the slack of561.
The paper also distinguishes its arbitrarily many factors from infinitely
many examples with a fixed number of factors. Its quantitative effective
threshold result uses the smaller range alpha<25/144 (Theorem5).

The transferable design is to retain a small group whose identity is
sufficient for the target predicate. It is not enough merely to find a
small group containing many objects.

## 5 1729, cube boundaries, and pyramid differences

For Chernick's three factors the identity is

    (6m+1)(12m+1)(18m+1)-1=36m(36m²+11m+1).

If all three factors are prime, their predecessor lcm is36m, so the
product is Carmichael. At m=1 this is1729=7*13*19 and1728=12³.
Independently,1729=1³+12³=9³+10³. The shared number gives a concrete
intersection of certificates; it does not make equal-cube identities
equivalent to Korselt or establish infinitely many prime triples.

The [pyramid packet](pyramid_successor_rows_20261014b.md) treats the user's
rows through their centre and their odd finite-difference moments. It
recovers the linear and square identities, corrects the old unit-width
root bracket, and replaces it with a bounded-width exact row compiler.
This is the productive role of the successor chain: taking differences
changes polynomial degree in a controlled way, whereas primality and
finite-field periods need their own sidecars.

### Polygonal receipts at every marked prime

[Kellner–Sondow, *On Carmichael and polygonal numbers, Bernoulli
polynomials, and sums of base-p digits*](https://math.colgate.edu/~integers/v52/v52.pdf)
gives the digit criterion (Theorem2.2) and a polygonal representation at
the greatest prime factor (Corollary4.3). Its primary Carmichael numbers
satisfy s_p(N)=p at every prime divisor;1729 is the first. These are CITED
connections, not novelty claims.

Here is a direct local formulation keeping EVERY marked prime. For odd
squarefree composite N and p|N, write N=pM. Since p=1 mod(p-1),

    p-1 | N-1 iff M=1+h(p-1), h a positive integer
                 iff N=G_{2h+2}(p),
    G_r(j)=j+(r-2)*j*(j-1)/2.

Thus all prime factors supply simultaneous even-rank polygonal receipts:

    561 =G_188(3)=G_12(11)=G_6(17)=T_33,
    1729=G_84(7)=G_24(13)=G_12(19).

Conversely, satisfying this condition at every prime factor gives Korselt.
One marked prime does not suffice:91=G_6(7), but the13-condition fails.
The cheap hostile is the same91 that exposes the missing CRT conductor.

Digit sums provide an equivalent local checksum, since s_p(N)=N mod(p-1).
Squarefreeness and compositeness rule out s_p(N)=1, so Korselt makes
s_p(N)>=p. For561 the digit heights (s_p(N)-1)/(p-1), ordered by3,11,17,
are(3,1,1); for1729 they are(1,1,1). This is an exact meaning for a
three-part profile with one exceptional coordinate, not a universal rule
for prime triples.

### Carry conservation and the limitation of a digit statistic

For nonnegative a,b and base p>=2, normalize the usual digit convolution.
Let K be the sum of its nonnegative outgoing carries, counted with their
integer magnitudes. Summing the column equations gives the exact law

    s_p(ab)=s_p(a)*s_p(b)-(p-1)*K.

The proof is a telescoping sum: each internal carry enters once and leaves
with coefficient p. This supplies the retained defect between a product
of digit sums and the digit sum of the product.

For odd n, let S(n)=(3n+1)/2^v2(3n+1). Binary multiplication by3 has
s_2(3)=2, and adding1 subtracts the number c of trailing ones of3n from
the putative digit increase. Division by2 removes zero digits. Therefore

    s_2(S(n))=2*s_2(n)+1-K(n,3)-c.

This is an exact accounting identity, not a positive drift theorem.
The native path7->11->17 grows while its digit count falls from3 to2.
Any proposed digit/Fisher/entropy ranking must retain information strong
enough to distinguish that example before it can supply payment.

## 6 Transfer to Collatz and the remaining quantifier

The [Collatz successor packet](collatz_successor_holonomy_20261014b.md)
shows why importing the group argument literally fails: growing native
blocks and strictly descending native blocks can both be the identity in
the full affine group modulo any prescribed modulus coprime to6.
Keeping modular translation as well as slope still does not decide payment.

It also constructs source-guarded descending blocks with valuations
1^k followed by2^(2k), and supplies exact ROOT subfamilies. Crucially k is
fixed by a source's v2(n+1); one cannot choose a convenient k after seeing
the desired modular return. The bank has an explicitly infinite uncovered
complement. Thus the next obligation is a complementary source guard with
proved descent or a grounded terminal, not another count of neutral words.

The connection record is precise. Source: Korselt-certified prime subsets.
Target: legal Collatz blocks. Map: composition followed by a finite quotient.
Preserved predicate: quotient identity. Destroyed information: exact drift,
source legality and terminal status. Needed sidecar: native cell, integer
carry, and descent/root receipt. Cheapest decisive test: a growing Mersenne
prefix with the same quotient identity. The simultaneous-prime example
motivates simultaneous block obligations, but supplies none of those missing
Collatz receipts by itself.

## Reproduction and scope

    python 04-computation/experiments/carmichael_successor_blocks_20261014b.py
    python -O 04-computation/experiments/carmichael_successor_blocks_20261014b.py

[Script](../../04-computation/experiments/carmichael_successor_blocks_20261014b.py)
and [JSON](carmichael_successor_blocks_20261014b.json):14,272 checks, identical
normal/optimized output. LF SHA256:
`c348998d575e3afaa76c8d2766c4bc52cff2ef215db65dc78cf3c72869285e3e`.
All arithmetic and primality checks use integer trial division. The complete
catalogs use the proved quotient bound. Independent controls choose q first
and factor Nq-1 for q<=2000; those bounded controls do not supply completeness.
Other controls compare all-residue power identities through699, one-prime
certificates for the nine Carmichael numbers below20000, non-squarefree
hostiles, a bad CRT conductor,100 Chernick parameters, digit criteria
through5000, marked polygonal receipts, carry conservation in six bases,
and1000 odd Collatz sources.

Independent peer audits passed the algebra and checker interfaces of all
four linked packets. The final saved outputs contain72,658 exact checks in
total (14,272 here;6,226 Collatz;17,426 prime roles;34,734 pyramids), with
independent optimized replays matching each JSON. The AGP asymptotic is
CITED, not proved by these finite controls. No unproved child is promoted
to a terminal by the integration.
