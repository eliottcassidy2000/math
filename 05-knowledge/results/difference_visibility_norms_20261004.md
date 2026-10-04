# Difference diagonals, primitive points, golden norms, and quadratic secants

2026-10-04 (America/Denver).

**PROVED:** elementary identities and density arguments below. **FINITE-EXACT:**
the declared computational controls. **CITED:** the classical zeta value,
the physical coupling convention, and scoped inherited quadratic results.
No identity predicting the Higgs coupling is established.

The anchor is a difference family that can transport arithmetic; the niche
is primitive lattice visibility; the wildcard is the proposed lambda.
The closest mechanism is the two-register reader in
[golden carry memory](golden_primes_carries_route_compiler_20261003.md).
The hostile is two points with the same difference but different gcd and
norm. The corrected near miss is primitive coefficients versus a ring unit,
recorded in MISTAKES on October 4. The least-used coordinate is the starting
point along a difference diagonal. Live concepts: differences, translation,
content, quadratic norm, secants, and sampling measure.

## 1. The dotted diagonal and the coprime density have an exact common model

In the positive square 1<=a,b<=N, transposition exchanges the two strict
triangles and fixes the N diagonal positions. Thus N^2=N+2*C(N,2).
An off-diagonal point with b>a has unique coordinates (a,d), d=b-a>0.
For fixed d,

    gcd(a,a+d)=gcd(a,d).

Consequently exactly phi_E(d) of every d consecutive values of a are
coprime, where phi_E is Euler's totient, not the golden ratio. The
one-dimensional density within that difference family is phi_E(d)/d.
At d=1 every pair is coprime; at d=4 exactly half are. No single
fixed-d family has been assigned the two-dimensional density 6/pi^2.

Let C_N count coprime ordered pairs in the square. Counting each strict
triangle by its larger coordinate and keeping its one coprime diagonal
point (1,1) gives

    C_N=2*sum_(b=1..N) phi_E(b)-1.

Alternatively the elementary divisor identity
1_(gcd(a,b)=1)=sum_(k|a,k|b) mu(k) gives

    C_N=sum_(k=1..N) mu(k)*floor(N/k)^2.

Divide by N^2. The floor error is bounded by 2/N times the harmonic sum
plus O(1/N), and the omitted tail of sum 1/k^2 tends to zero. Hence

    lim C_N/N^2=sum_(k>=1) mu(k)/k^2
               =prod_p(1-p^-2)=1/zeta(2)=6/pi^2.

The absolutely convergent Euler product follows by expanding squarefree
divisors. The value zeta(2)=pi^2/6 is the classical identity tabulated in
[NIST DLMF 25.6.1](https://dlmf.nist.gov/25.6.E1). The limit is thus an
actual prime-by-prime statement about the user's plane. It is a density
under the specified square sampling, not a primality test for one integer.

The diagonal family d=0 is special: gcd(a,a)=a, so only a=1 is primitive.
Negative differences give the reflected triangle. Keeping the orientation
of d is necessary if later operations distinguish the two sheets.

## 2. The examples 7--11 and 8--12 expose the missing coordinate

Write a point on the difference-d diagonal in O=Z[phi] as

    alpha=a+(a+d)phi=phi^2[(a-d)+d phi].

The equality uses phi^2=phi+1, and multiplication by phi^2 is a unit
operation. Its field norm is

    N(alpha)=a^2-ad-d^2.

At the user's two gap-four points,

    7+11phi=sqrt(5)*phi^5,       N=5,  gcd(7,11)=1,
    8+12phi=4phi^4,              N=16, gcd(8,12)=4.

Here sqrt(5)=2phi-1 has norm -5 and phi has norm -1. The first principal
ideal has index 5 and the second is (4), with index 16. The same difference
therefore does not determine content, norm, ideal, or ring-unit behavior.
This is a concrete quotient-loss test, not a rejection of the family: its
faithful state is (a,d).

Translation by one place along that family adds 1+phi=phi^2. Multiplication
by phi instead sends

    (a,d)->(a+d,a).

It preserves gcd(a,d) and changes the norm's sign, while usually changing
the difference family. This is the exact Fibonacci recursion on the two
coordinates. Discarding a prevents this recursion from descending to d
alone. The [difference-operation note](difference_family_operations_20261004.md)
gives an alternative crossed multiplication that does descend to signed
differences, with a retained translation coordinate for raw pairs.

## 3. Primitive vectors and invertible golden residues are different tests

Modulo a prime p there are p^2-1 nonzero coefficient pairs. These are the
primitive pairs modulo p. The number of units in O/(p) instead is

    p^2-1            if X^2-X-1 is irreducible modp,
    (p-1)^2          if it has two distinct roots modp,
    p(p-1)           if p=5 (the repeated-root case).

These counts follow respectively from a quadratic field, a product of two
fields, or a field with a square-zero direction. The
[proved clock note](golden_prime_clocks_20261003.md) supplies the norm-unit
criterion and the split/inert distinction. For p=5 or 11, only 5/6 of the
primitive pairs are units; modulo 19 the fraction is 9/10. Modulo 2,3,7
all primitive pairs are units.

On a fixed difference d, the norm polynomial is a^2-ad-d^2. For p not
dividing d, it has zero roots at an inert prime, two at a split odd prime,
and one at p=5. For p dividing d it has exactly the one root a=0.
Completing the square gives discriminant 5d^2 for odd p; p=2 is checked
directly. These are exact local filters for the family, without claiming
that surviving values are prime or that infinitely many are prime.

## 4. The diagonal is also the derivative boundary of a pair recursion

For F_c(x)=x^2+c and a pair x,y, let s=x+y and d=y-x. Then

    F_c(y)-F_c(x)=sd,
    F_c(y)+F_c(x)=(s^2+d^2)/2+2c.

Thus the entire pair evolution is polynomial in (s,d), and transposition
is d->-d. The difference update itself does not contain c, but it needs
the sum, whose update does contain c. This identifies exactly what a
difference-only description loses.

On the dotted diagonal d=0, the transverse derivative is 2x. Away from
it, the divided difference is x+y. For a nonconstant periodic orbit
x_0,...,x_(k-1), applying the first formula to adjacent points and
telescoping their nonzero differences proves

    product_i (x_i+x_(i+1))=1.

This product is generally different from the dynamical multiplier
product_i 2x_i. For the rational cycle of x^2-29/16,

    -7/4 -> 5/4 -> -1/4 -> -7/4,

the adjacent sums are -1/2,1,-2 and their product is 1; the derivative
product is 35/8. These are inherited exact results from
[THM-4146, order-six lift and fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md),
section 2. Its projective three-cycle has an SL_2 lift of order six;
the central sign is another explicit coordinate lost by a quotient.

For x^2-7/4, the period-three cubic is

    P(x)=x^3+x^2/2-9x/4-1/8.

Exact reduction modulo P gives F_c^3(x)=x and both products equal to 1.
The remainder of P modulo the fixed-point polynomial x^2-x-7/4 is x+5/2,
and that fixed-point polynomial takes value 7 at -5/2. Thus their gcd is
one, excluding fixed roots. Discriminant(P)=49>0 proves that its three
roots are distinct and real, so they have exact period three.
The [earlier seventh-root bridge](arithmetic_braids_20260917_geometry.md)
identifies this parabolic cycle. Its coincidence with a derivative
multiplier of one does not transport to the rational -29/16 cycle.
The critical orbit starting at zero is a different object and is not
asserted to be that three-cycle. No map identifying either quadratic
system with Collatz has been established here.

There is also an exact trigonometric reading: z=x+1/2 satisfies
z^3-z^2-2z+1=0, whose roots are 2cos(k*pi/7), k=1,3,5. The map becomes
z->z^2-z-1, equal to z^3-3z on these roots, hence angle tripling. This
is a real cubic cyclotomic field with discriminant 49. It is a different
field from the quadratic golden field with discriminant 5; a three-cycle
alone does not identify those arithmetic structures.

The user's triangular split is the same polarization identity in another
chart: for x=A+1,y=B+1, the sum is N and xy=(N^2-d^2)/4. The diagonal and
two sheets can therefore carry addition, multiplication, secants, or
derivatives, provided their coordinate meanings are retained.

**Concurrent anchor bridge, with a carry-loss probe.** The incoming
[rational-anchor compiler](rational_anchor_returns_20261004.md) uses the odd
word (1,1,2), whose rational fixed point is -19/11, total halving count4,
and ordinary parity length7. The exact identity

    N(phi^7-1)/2^4=-29/16

therefore connects these two specified parameters. It follows directly from
phi^7-1=7+13phi, of norm -29. This is an algebraic parameter assignment,
not a dynamical conjugacy. The same chart's exact golden fraction is

    Theta(-19/11)=(phi^6+phi^4+phi^2)/(phi^7-1)
                 =4phi^4/(phi^7-1)=(-4+20phi)/29.

Its full denominator ideal is

    (phi^7-1)=(Phi_7(phi))=(29,phi-24).

For the first equality with the actual denominator ideal, multiplying Theta
by phi^7-1 is integral, and the proper denominator ideal therefore contains
the maximal ideal of index29; they are equal. The cyclotomic equality uses
the unit phi-1. Evaluation at phi=24 modulo29 annihilates7+13phi and has
index29, proving the final equality. This gives another exact prime-ideal
connection, while its arithmetic cycle remains rational with denominator11.

The proposed general parameter map
c(w)=N(phi^(len(w)+sum(w))-1)/2^sum(w) depends only on the two counts and
forgets the ordered carry. A cheap hostile is the pair of primitive words
(1,1,1,3) and (1,1,2,2): both give c=-121/64, while their disjoint rational
odd cycles have anchors -65/17 and -73/17. The exact script checks both
cycles. The strongest survivor is the displayed identity for the specified
word; any extension that recovers the ordered word or distinguishes these
arithmetic cycles needs more than this norm parameter. A noninjective
semiconjugacy could intentionally merge cycles, so the collision alone does
not exclude every dynamics-preserving map.
This probe supplies neither a conjugacy nor an explanation of the quadratic
cycle's stability.

## 5. The lambda expression needs its parentheses and physical convention

The expression in the current request, read literally, and the one in the
earlier note are different:

| Formula | Decimal illustration |
| --- | ---: |
| (5 phi^4)^(-1/6) | 0.5548555133444340 |
| 1/[5(2 phi^4)^(1/6)] | 0.1292805634433970 |
| 1/[5(phi^4)^(1/6)] | 0.1451125260492653 |

Each can be treated exactly. If x^6=phi^-4/C, then

    C^2 x^12-7C x^6+1=0,

since phi^4+phi^-4=7 and their product is 1. The three C values are
5,31250,15625, respectively. This is an annihilating equation; no
minimal-degree claim is required.

In the convention V(H)=m_H^2 H^2/2+lambda_3 v H^3+lambda_4 H^4/4,
the Standard Model tree-level relation is
lambda_3=lambda_4=m_H^2/(2v^2), approximately 0.13. This is stated in
the [CMS primary paper, equations 1--2](https://cds.cern.ch/record/2904902/files/2407.13554.pdf).
The older golden expression is numerically near that dimensionless value.
The cited physical model does not derive the golden expression, and none
of the lattice-density or clock arguments here supplies such a derivation.
A proposed physical relation must fix the potential convention and, beyond
tree level, its matching/renormalization prescription before comparison.

The earlier [circulant note](collatz_circulant_20260930_circulants_lucas_cubic_monotile.md)
already refutes equality of its lambda with the level-five Syracuse
Fourier maximum. That repaired boundary remains in force. The productive
prime connection in this session is the proved coprime-density and
quadratic-ring analysis, with its specified sample space.

## Reproduction

    python -X utf8 -B 04-computation/experiments/difference_visibility_norms_20261004.py
    python -X utf8 -B -O 04-computation/experiments/difference_visibility_norms_20261004.py

The [saved output](difference_visibility_norms_20261004.out) records nine
square sizes through 10000 using independent totient/Mobius counts, with
literal pair/gap enumeration through 300; 800 complete fixed-gap period
tests; 10201 signed golden-coordinate tests; unit counts at every prime
through 97 and 750 fixed-gap root controls; 3362 quadratic pair tests and
two independent exact cycle checks; and the lambda annihilating equations.
Only the explicitly illustrative lambda decimals use approximate arithmetic.
All acceptance checks remain active under optimized Python.
