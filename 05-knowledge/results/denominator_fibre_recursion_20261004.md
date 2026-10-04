# Golden denominator fibres: unique periodic lifts and recursive cycle towers

2026-10-04 (America/Denver).

**PROVED:** the all-denominator periodic-lift theorem, the binary and ternary
cycle towers with an explicit child selector, and the signed certified-root
splice formulas. **FINITE-EXACT:** the stated lattice/phase censuses and
integer-realization filters. **CITED:** the classical value zeta(2)=pi^2/6,
through the visibility note's primary reference. **OPEN:** assigning a rational golden passport
to every integer, and universal signed Collatz coverage. No historical
novelty claim or new theorem ID is assigned here.

## Inheritance and the connection being tested

The closest mechanism is the fixed-denominator golden lattice in
[the signed carrier note, sections 4–5](collatz_golden_carriers_20261004.md).
It conserves the full denominator ideal and gives complete finite gates at
q=1,2,11,76. The underused coordinate is the numerator **phase orbit**,
which retains more than the denominator ideal. The canonical hostile is
the rational Collatz cycle 1/13: its golden value has the same denominator
ideal as the -5 cycle. The corrected near miss is the boundary convention
at the golden value phi^(-1); the upper convention must be retained.

The [prime-clock proof, section 3](golden_prime_clocks_20261003.md) supplies
the exact orders at powers of 2 and 3, without a blanket prime-power law.
[THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md)
and [the route-memory splice, section 4](collatz_atom_route_memory_20261003.md)
already supply finite-head realization and explicit common-tail recursions.
The signed splice below repackages that mechanism with its golden phase;
it is not presented as a new proof of positive-prefix density.

The live board is a periodic golden point, its torus phase, denominator
ideal, a prime-power lift, an ordered arithmetic word, and an integer root.
The anchor is an exact recursive family; the niche is its finite phase
selector; the wildcard is a common arithmetic/golden fixed point for all
four certified basins. The two recursion directions have different scopes:
raising a denominator preserves golden dynamics, while certified-root
splicing preserves actual integer arrival. Their distinction is tested,
not inferred from a matching count of four objects.

## 1. The phase invariant

Write

    phi=(1+sqrt(5))/2, beta=phi^(-1), psi=1-phi=-beta,
    O=Z[phi], G=[[0,1],[1,1]].

Use the upper golden map on the closed interval:

    B(x)=phi*x-d,  d=1 if phi*x>1, and d=0 otherwise.

For x=(a+b*phi)/q, q>0, a,b integers, one step is

    (a,b) -> (b-d*q,a+b).

Thus the numerator phase v=(a,b) mod q evolves by **G alone**, independently
of the digit. Since G is invertible, the G-orbit of v is invariant along
the golden trajectory. Its full denominator ideal

    D(x)={u in O: u*x belongs to O}

is also invariant: multiplication by phi is a unit and the digit is
integral. The phase orbit is generally a strictly finer invariant.
The exact scalar denominator is q precisely when gcd(a,b,q)=1.

## 2. Every nonzero rational phase has exactly one periodic lift

**Theorem (all q).** Projection x -> x+O is a bijection from the nonintegral
periodic points of B onto the nonzero classes of Q(phi)/O. The zero class
has exactly three periodic lifts: 0 fixed, and the cycle 1 <-> beta.
For q>1, periodic points of exact scalar denominator q correspond
bijectively to primitive pairs (a,b) modulo q. The exact golden period is
the exact G-period of the pair.

Here “nonintegral” means x not in O, rather than x not in Z. Every periodic
point of B is automatically in Q(phi), by solving its finite affine return
equation, so no periodic points are omitted by the domain in the statement.

**Conjugate window.** On a periodic orbit, extend the digits periodically
into the past. Since |psi|<1,

    x'=-sum_(k>=0) psi^k d_(-1-k),
    -phi <= x' <= 1.                                 (1)

The extreme bounds follow by independently selecting the even or odd
powers; using admissibility can only reduce the range. There is one crucial
additional restriction. If x>beta, its preceding digit cannot have been 1,
because a digit-1 image lies in [0,beta]. Consequently

    x>beta  implies  x'=psi*x_previous' >= -beta.      (2)

These are two rectangles in the real/conjugate plane. The branch restriction
(2) is essential to the overlap argument; the enclosing rectangle alone
would leave false possible collisions.

**Uniqueness.** Suppose periodic x,y have the same phase class. Put
delta=x-y=a+b*phi in O. Then |delta|<=1 and |delta'|<=phi^2. Hence

    |b|=|(delta-delta')/sqrt(5)|<=phi<2.

Enumerating b=-1,0,1 gives exactly

    delta in {0, +/-1, +/-beta, +/-beta^2}.           (3)

Exchange x,y if necessary to make delta positive.

* delta=1 forces (x,y)=(1,0).
* delta=beta^2 has delta'=phi^2. Equality in (1) forces
  x'=1,y'=-phi, hence (x,y)=(1,beta).
* delta=beta implies x>=beta and delta'=-phi. If x=beta, y=0.
  If x>beta, (2) and y'<=1 force x'=-beta,y'=1. But the conjugate
  of -beta is phi, which would give x=phi>1, a contradiction.

Every nonzero collision in (3) therefore occurs only among 0,beta,1,
which are indeed periodic and lie in O. Any additional periodic point in
the zero class would collide with 0, so that exception is complete.

**Existence.** Choose a representative x in [0,1] of any rational class,
by subtracting an ordinary integer. Its scalar denominator stays fixed.
The conjugate recurrence x'_next=psi*x'-d eventually enters |x'|<=3
and remains there. Together with 0<=x<=1 and fixed denominator this is
a finite lattice set, so the orbit is eventually periodic. The phase
follows an invertible finite permutation. Rotate the terminal golden
cycle until its phase equals the original phase. This gives a periodic
representative of the original class. A nonzero initial phase can never
become zero.

Finally, if the phase has period ell, B^ell(x) is periodic and has the same
phase as x. Uniqueness gives B^ell(x)=x. Conversely any golden period is a
phase period. Thus the exact periods agree. This proves all the claims.

**Counting consequence.** For q>1 the number of periodic points with exact
scalar denominator q is the Jordan totient

    J_2(q)=q^2 * product_(prime p divides q)(1-1/p^2). (4)

The condition gcd(a,b,q)=1 excludes pairs simultaneously divisible by a
prime factor of q, proving (4) directly by inclusion-exclusion. Periodic
points with denominator dividing q number q^2+2: the q^2 torus classes
have one lift each, except zero has three. A particular phase can have
period smaller than the full matrix order; q=11 below is a required
hostile against identifying those two periods.

**Prime-density corollary.** Among all q^2 numerator phases the fraction
with exact scalar denominator q is

    J_2(q)/q^2 = product_(prime p divides q)(1-p^(-2)).

This depends on the prime support of q, rather than its exponents. Along
the primorials Q_k, the products of the first k primes, the density tends
to product_p(1-p^(-2))=1/zeta(2)=6/pi^2. The Euler-product mechanism and
the cited classical value of zeta(2) are recorded in
[the visibility note, section 1](difference_visibility_norms_20261004.md).
This is the same local two-coordinate primitivity test as visible lattice
points. It is a primorial limit, not a limit for arbitrary q tending to
infinity: prime powers retain the constant factor 1-p^(-2), while prime
denominators themselves have density tending to1.

The projection intertwines B with multiplication by phi on Q(phi)/O.
It loses the real/conjugate lift, but the theorem reconstructs that lift
uniquely on periodic nonzero phases. For arbitrary transient points a
phase alone does not specify the current point or the next parity digit.
The script exports `periodic_lift(q,phase)`, implementing the existence
proof by starting from a reduced real representative and rotating its
terminal cycle to the requested phase. Zero is its documented canonical
choice for the exceptional class.

### A commutative operation, with its necessary boundary tag

Remove the two duplicate zero-class representatives beta and1, retaining0.
The remaining periodic points identify with the additive torsion module
Q(phi)/O, isomorphic as an abelian group to (Q/Z)^2. Therefore

    x (+)_per y = the canonical periodic lift of [x+y]

is an associative, commutative group operation. Scalars in O act through
the quotient, and multiplication by phi agrees with B. This is a concrete
operation on complete periodic objects, not ordinary addition of their
real coordinates. For example (phi/2) (+)_per (phi/2)=0.

The canonical choice deliberately forgets the distinction between the
zero orbit and the upper-map cycle 1<->beta, which represents the signed
integer cycle through -1. A faithful signed carrier must keep that boundary
tag separately. The group is an O-module, **not a quotient ring**: although
1/2 and3/2 represent the same class, multiplying each by1/2 gives1/4 and3/4,
whose difference1/2 is not in O. Ordinary multiplication of arbitrary
classes is therefore not well-defined. Integral O-scalar multiplication is.

There is a stronger obstruction to repairing the product on this same
additive object. The group H=Q(phi)/O is torsion and divisible. For any
biadditive operation mu:H*H->H, choose n>=1 with n*a=0 and choose c with
n*c=b. Then

    mu(a,b)=mu(a,n*c)=mu(n*a,c)=0.

Thus every distributive product on this exact additive quotient is zero,
even if it is invented independently of ordinary multiplication. This is
a statement about H; it is not a no-go for a richer carrier.

One richer object preserves ordinary commutative ring operations. Fix a
prime p and retain a compatible sequence v_a in O/(p^a), with
v_(a+1)=v_a modulo p^a. Define sums and products at every level. Reduction
commutes with both operations, so compatible sequences form a commutative
ring with the levelwise identity. All phases, including zero and nonunits,
belong to this object. At the golden-point level, compatibility is the
parent operation “multiply by p, reduce modulo O, then take the canonical
periodic lift.” The zero-class boundary tag remains separate.
The golden shift is multiplication by the unit phi. It is an additive
automorphism, not a ring automorphism: the product of two shifted values
has a factor phi^2, whereas shifting their product has only a factor phi.

A phase point has p^2 next-level lifts. For primitive phases in the binary
and ternary towers below, passing to whole cycles groups these into p
children and discards rotation. Multiplication alone does descend to
phi-orbits, since (phi^i*u)(phi^j*v)=phi^(i+j)*uv. **Addition does not:**
modulo2, 1 and phi are in the same nonzero orbit, while 1+1=0 and
1+phi=phi^2 is nonzero. Thus full phases or anchors are needed for the
ring structure, though the full collection of phi-orbits at a fixed
modulus retains a multiplicative monoid.
Compatibility and this extra coordinate do not by themselves certify
integer Collatz realization or universal arrival.

## 3. The binary and ternary micro-to-macro towers

Both X^2-X-1 modulo 2 and modulo 3 are irreducible. A primitive pair modulo
p^a, p=2 or 3, is therefore a unit of O/(p^a): its norm is nonzero modulo
p. If phi^L*v=v for this pair, multiplication by its inverse gives
phi^L=1. Thus every primitive phase has the full clock period. The inherited
clock theorem gives

| Exact denominator | Period of every golden cycle | Number of cycles |
|---|---:|---:|
| 2^a, a>=1 | 3*2^(a-1) | 2^(a-1) |
| 3^a, a>=1 | 8*3^(a-1) | 3^(a-1) |

The counts follow from (4), divided by the common period. In particular
q=2 starts with one 3-cycle, and q=3 starts with one 8-cycle. This is an
all-height theorem for golden cycles, not a census extrapolation.

### The parent map and an explicit child selector

Fix p in {2,3}, q=p^a, and one parent phase cycle of length L. Reduction
modulo q maps the primitive phases modulo pq onto the parent system.
Every parent phase has p^2 lifts. Every child cycle has length pL, so
exactly p child cycles lie over each parent cycle, each covering it p times.

This count can be implemented without an unlabelled search. Choose an
anchor v on the parent cycle, for example its lexicographically least
representative. Write

    G^L=I+q*U,  c=U*v mod p.                         (5)

The matrix U is invertible modulo p. For p=2, the first level has
G^3-I=2G and U=G modulo 2. At the next level G^6-I=4G^3, and U=I modulo2;
successive squaring preserves that residue. For p=3, the first level is

    G^8-I=3[[4,7],[7,11]],

whose quotient determinant is 1 modulo3; successive cubing preserves
the quotient modulo3. Hence c is nonzero for a primitive v.

The p^2 lifts of the anchor are v+q*w, w in F_p^2. After L steps,

    w -> w+c.                                       (6)

They split into p parallel translation orbits, each with p points. The
explicit label

    child(w)=det(c,w) mod p                           (7)

is constant on each orbit and distinguishes all p children. Labels depend
on the chosen anchor and coefficient basis; the child cycles do not.
The script exports `lift_phase_children(prime,q,parent)` and verifies both
the affine return (6) and the complete cycle partition independently.

For the first binary split take parent v=(0,1) modulo2. Then c=(1,1),
and the child label is w_2-w_1 modulo2. The two child cycles have anchor
phases (0,1) and(0,3) modulo4, whose exact periodic golden lifts are
phi/4 and -1+3phi/4. Both have period6, and both have parent phi/2.

On periodic golden points the parent map is equally precise: take p*x
modulo O, then its unique periodic lift. It commutes with B and implements
the same p-to-1 covering of cycles. Thus cycles form a regular binary or
ternary rooted tree across these denominator levels. The zero-class
exception is avoided by starting at denominator p, rather than trying to
choose a parent among the three zero-class lifts.

**Decisive integer hostile.** The q=2 parent is the golden image of the
ordinary Collatz cycle through 1. Its two q=4 children are the golden
images of rational arithmetic cycles with odd nodes

    1/29,                 and                 5/7 <-> 11/7.

Neither child is an integer orbit. Denominator lifting preserves golden
dynamics, phase, and the parent covering relation; it does not preserve
the signed integer-realization predicate. Carrying that predicate is
necessary before calling a child a new certified integer basin.

## 4. What q=11 and q=76 retain beyond the ideal

At 11, phi has roots 4 and 8, of orders 5 and 10. A primitive phase has
two evaluation coordinates u=a+4b and v=a+8b modulo11. The possibilities
are:

| Phase support | Number of phases | Exact period | Golden cycles |
|---|---:|---:|---:|
| u nonzero, v=0 | 10 | 5 | 2 |
| u=0, v nonzero | 10 | 10 | 1 |
| both nonzero | 100 | 10 | 10 |

The first row is denominator ideal (11,phi-4). Its two cycles are exactly
the inherited -5 and 1/13 examples. For Theta(-5)=(-1+7phi)/11,
u=5; its orbit is the five nonzero quadratic residues {1,3,4,5,9}.
For Theta(1/13)=(1+4phi)/11, u=6; its orbit is the five nonresidues.
Multiplication by 4 preserves and transitively visits each set. Thus the
phase orbit separates the two examples that the full ideal merges.
This character distinguishes these two already specified period-five
components; it is not a primality test or a general Collatz classifier.

At 76=4*19, a primitive phase is a unit at 4, so its period is divisible
by 6. At 19 the roots of X^2-X-1 are 5 and 15, with orders 9 and18.
At least one component is nonzero, so its period is divisible by9.
Every primitive phase therefore has period exactly18. There are

    J_2(76)=4320 phases, hence 240 golden cycles.

Those with denominator ideal (76) number 216 cycles; the two ideals with
only one of the 19-components present give 12 cycles each. The integer
-17 cycle belongs to the first group. The finite ordered-word realization
filter, not the ideal or cardinality, identifies it among these 240.

### The user's 105: two equally numerous families of different periods

At q=105=3*5*7, primitive means that the two-coordinate phase is nonzero
in each of the three residue fields/rings. At 3 the 8 nonzero phases all
have period8. At 7 the 48 nonzero phases all have period16; this follows
from the inherited inert-prime clock and invertibility of the phases.

Modulo5, write G=3I+N, with N nonzero and N^2=0. The four nonzero phases
in ker(N), the line b=3a, have period4 because 3 has order4 modulo5.
The remaining20 phases have period20: G^4=I+3N, whose nonzero translation
on such a vector takes five iterates to return; applying N to any return
first forces its length to be divisible by4. These are exactly the
units modulo5. CRT gives

| Local phase at 5 | Primitive105 phases | Global period | Golden cycles |
|---|---:|---:|---:|
| nonzero eigenline | 8*4*48=1536 | lcm(8,4,16)=16 | 96 |
| off the eigenline | 8*20*48=7680 | lcm(8,20,16)=80 | 96 |

Thus there are exactly192 golden cycles of scalar denominator105, with
96 of each period. The primitive phase density is
9216/105^2=1024/1225. The 7680 long-period phases are the units of O/(105);
the 1536 short-period phases are primitive pairs but nonunits, their norm
being divisible by5. Keeping both coordinates and the unit test explains
why these populations differ. No claim assigns these golden cycles to
integer Collatz orbits; the independent ordered-word filter finds no
integer cycle in this exact finite denominator fibre.

## 5. Fixed-denominator recursion that really carries an integer route home

The vertical tower changes q and can lose integer realization. The
horizontal operation below retains a chosen certified root and its q.
It is the signed version of the inherited splice, with a synchronized
golden coordinate and an exact phase clock.

Let w be any finite positive valuation word, possibly empty. Put

    p=len(w), A=sum(w),
    U_w(n)=(3^p*n+S)/2^A,
    C=2^A+3S, T=2*3^p, R=2^T.

For each r in {1,-1,-5,-17}, choose a>=1 satisfying

    2^(A+a)*r=C mod3^(p+1).                         (8)

There is exactly one residue class of a modulo T. Here r and C are units
modulo3, and 2 generates the units modulo3^(p+1). A direct proof of the
latter uses 4=1+3 and induction
v3(4^(3^j)-1)=j+1; together with the sign modulo3 this gives exact order
2*3^p, the full unit count. A sufficiently large representative a0 makes
every preterminal odd node have the sign of r and absolute value >91.
This excludes padding by any of the known odd-cycle vertices.

For every k>=0 define

    a_k=a0+T*k,
    n_k=(2^(A+a_k)*r-C)/3^(p+1).                    (9)

**PROVED:** n_k realizes exactly the odd valuation word w followed by a_k
and then arrives at the chosen root r. To see legality, start from r and
invert the combined word; integrality of the full expression forces each
successive inverse division by3, and positive division exponents make the
intermediate numerators odd. For negative targets the inverses stay
negative. For the positive root the sufficiently large a0 makes every
intermediate positive. Increasing a by T preserves all divisibility and
increases these magnitudes. The last edge is
(2^a*r-1)/3 -> r, with exact valuation a.

The recurrence on actual integers is

    n_(k+1)=R*n_k+(R-1)*C/3^(p+1).                  (10)

Its common fixed point is xi=-C/3^(p+1), independent of r. The additive
constant is integral by the unit-clock identity. Each sequence diverges
in real magnitude but converges 2-adically to xi. The point xi is the
rational terminating shadow whose formal combined word ends at 0; it
is not an additional integer root.

Let eta be the finite golden value of the raw parity word for w followed
by one digit1, with zero tail. Let L_k=A+p+a_k+1. The actual certified
coordinates satisfy

    Theta(n_k)=eta+beta^L_k*Theta(r),
    Theta(n_(k+1))=eta+beta^T*(Theta(n_k)-eta).        (11)

Both multipliers in the golden operation are units in O and eta belongs
to O. Hence (11) preserves the full denominator ideal and its phase orbit.
All four root families approach the same real golden boundary eta, while
their exact denominators remain respectively 2,1,11,76. The terminating
boundary itself requires the signed carrier's boundary warning; taking
this limit does not preserve the denominator or integer-root label.

The phase of (11) advances by G^(-T). If the root phase has period ell,
the exact sibling-index period is ell/gcd(ell,T). The four root phase
periods are 3,1,5,18, respectively. For example, at head length p=1 the
sibling periods are 1,1,5,3; at p>=2 they are 1,1,5,1.

For the common head w=(1), explicit first members chosen outside all
known odd cycles are

| Root | Source | Last valuation | Exact denominator |
|---:|---:|---:|---:|
| 1 | 227 | 10 | 2 |
| -1 | -1821 | 13 | 1 |
| -5 | -285 | 8 | 11 |
| -17 | -3869 | 10 | 76 |

Each has initial valuation1 and then the displayed final edge to its root.
This generalizes to every finite head. In particular no finite valuation
head, even with a negative sign retained, can distinguish the three known
negative basin ideals: every such head has arbitrarily large-magnitude
certified instances in all three. This is a precise obstruction to a
finite-prefix-only passport rule, not an obstruction to an algorithm that
also retains the full integer source or an unbounded certificate.

## Verification, transfer contract, and unresolved target

**Concurrent integration.** Incoming commit952b93f75 independently proves the
phase correspondence by an L-shaped lattice fundamental domain in
[difference families](difference_families_20261004.md), and proves a third
tower q=4*19^k with period18*19^(k-1),240*19^(k-1) cycles, and19 children per
parent. Its identity G^18-I=76*G^9 has an invertible leading matrix modulo19,
so the proof applies to primitive nonunits as well as units. Our selector
extends verbatim: U=(G^L-I)/q is invertible modulo19, and det(Uv,w) modulo19
labels the child cycles above an anchored parent. The incoming proof and
this extension received an independent audit; no integer-realization claim
is inferred from that tower. The two notes' phase theorem is the same
mathematical result with independent elementary proofs.

Reproduce with:

    python -X utf8 -B 04-computation/experiments/denominator_fibre_recursion_20261004.py
    python -X utf8 -B -O 04-computation/experiments/denominator_fibre_recursion_20261004.py

The [script](../../04-computation/experiments/denominator_fibre_recursion_20261004.py)
and [saved output](denominator_fibre_recursion_20261004.out) compare the
complete inherited real/conjugate lattice with an independent primitive
phase enumeration, and reconstruct every periodic lattice point from its
exact rational geometric series. The universe is q=1,...,40,64,76,81,105.
The seven overlap vectors in (3) are independently enumerated in an
oversized integer box. The binary tower is checked through q=64 and the
ternary tower through q=81, including every parent/child partition at the
intervening levels and the explicit determinant selector.
The q=105 census independently checks the complete CRT split above;
the prime-product formula is checked against every phase census, and
exact primorial-density fractions are retained through last prime13.

The signed splice controls use all 40 heads of lengths0,...,3 with letters
1,...,3, all four roots, and k=0,...,4: 800 actual signed integer routes.
Independent ordinary iteration verifies every bit and endpoint. Closed
arithmetic and golden recursions, full annihilator kernels modulo q,
phase-orbit preservation, and sibling-clock periods are separately checked.
Normal and optimized Python outputs agree; no assertion is disabled by -O.
The all-q uniqueness/existence proof received an independent agent audit,
including the overlap cases, zero-phase exception, and exact-period step.
The final audit also checked the child selector and q=105 split. A boundary
API correction validates q>0 before factorization, preventing a zero-input
loop; five rejected-input controls cover the exported selector's boundary.

The target map from a periodic phase to a golden point preserves period
and the golden dynamics; it does not assert signed integer realization.
The target map from a signed splice to an integer preserves its actual
ordered route and chosen root; it does not cover an arbitrary input integer.
For signed integers whose golden coordinate is known to be rational, the
phase theory supplies a finite eventual-cycle target. Computing or proving
such a rational passport for every input remains the global obligation.
