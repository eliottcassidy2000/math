# What a small return diagram forgets

2026-10-12. **PROVED:** the return-thickness and finite boundary-budget
statements below. **FINITE-EXACT:** symbolic identities and the specified
specializations. These give tests for proposed recursive constructions;
they do not settle Collatz, planar JC, or the plane chromatic number.

## Inheritance and connection contract

The closest proved mechanism is [THM-4146, quadratic three-cycle lift and
fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md).
It distinguishes a permutation of three points, its order-six linear lift,
and an actual quadratic return. The canonical hostile is its rational cycle
with derivative multiplier35/8 despite secant product1. The corrected near
miss is equality on a divisor promoted to an identity on a neighbourhood.
The least-used sidecar here is the multiplicity of the return defect.

The anchor is the three requested proof obligations; the niche is marked
quadratic dynamics; the wildcard is the arithmetic exception63. The board
is **marked state / return defect / critical value / boundary support /
source charge / terminal**. The roots, primes, and colour classes are
different objects. The comparisons below retain their actual maps.

## 1. Exact classification of the first missing neighbourhood

Work in characteristic zero, over a field containing a chosen three-cycle
of distinct points of f_c(x)=x^2+c. Its trace sigma determines

    c=-(sigma^2+sigma+2),
    P_sigma(x)=x^3-sigma*x^2-(sigma^2+2sigma+3)*x
                 +sigma^3+2sigma^2+3sigma+1.

These formulas are inherited from THM-4146. Put eta=2sigma+1 and
q=eta^2+eta+7. The discriminant of P_sigma is q^2; distinctness implies
q!=0. Direct factorization gives

    f_c^3(x)-x=(f_c(x)-x) P_sigma(x) P_(-sigma-1)(x).       (1)

Define the return thickness as the largest integer h such that P_sigma^h
divides f_c^3(x)-x. Then

    h=2 at eta=0, c=-7/4;
    h=1 for every other distinct chosen three-cycle.       (2)

Here h=1 remembers the three points; h=2 also remembers first-order
perturbations at all three points. This is a precise use of an infinitesimal
neighbourhood, rather than a statement about an interval of real points.

**Proof.** The derivative multiplier is constant on the cycle and equals

    mu=8 product(cycle points)=1-7eta-eta^2-eta^3,
    mu-1=-eta*q.                                         (3)

At eta!=0, q!=0 makes every cycle root a simple root of f_c^3-x. At eta=0,
the two cubics in(1) coincide. The resultant of P_sigma and f_c-x is q,
which here equals7, so there is no third factor at a cycle root. Hence the
multiplicity is exactly two. QED.

For the two requested parameters this becomes

| Chosen cycle | Multiplier | Exact return thickness |
|---|---:|---:|
| c=-7/4, P=x^3+x^2/2-9x/4-1/8 |1|2|
| c=-29/16, points(-7/4,5/4,-1/4) |35/8|1|

The value -7/4 in the second row is a **point**, whereas in the first row
it is the parameter. The first cycle is cubic algebraic, the second rational.
[THM-4139, rational graph and horizontal-carrier separation](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md)
gives the complete rational graph in the second case.

**Repetition cannot buy the missing order.** At c=-7/4, write H=f_c^3 and
H-x=P^2 R, with gcd(P,R)=1. Polynomial expansion shows inductively

    H^m(x)-x = m P(x)^2 R(x) mod P(x)^3, m>=1.             (4)

The characteristic-zero coefficient m does not vanish, so every repeated
return still has thickness exactly two. At c=-29/16 the return multiplier
is (35/8)^m!=1, so every repeated return still has thickness exactly one.
This assertion is for these specified cycles. In other complex cycles a
root-of-unity multiplier can make a later return more degenerate.

Thus neither repeating the rational three-cycle nor repeating the parabolic
three-cycle creates the next missing neighbourhood. A proof requiring that
neighbourhood needs an operation changing the defect, not a longer word
representing the same return.

## 2. The affine reversal preserves points, but changes the defect

[Affine reversal of quadratic three-cycles](arithmetic_braids_20260917_geometry.md)
already proves L(x)=-x-1/2, c'=c-sigma-3/4 and eta'=-eta-1.
An exact formula now identifies the information lost:

    f_(c')(L(x))-L(f_c^2(x))=(x+sigma) P_sigma(x).          (5)

It vanishes on the cycle, but is not divisible by P_sigma^2. Consequently
this inherited pointwise reversal is not an identity on the full first
neighbourhood. Its two multipliers obey

    mu+mu'=q+2,
    mu-mu'=-(2eta+1)q,
    mu*mu'=-q^3+7q^2+q+1.                                (6)

Because q!=0, the two multipliers are equal precisely at eta=-1/2,
the chosen rational AP cycle at c=-29/16. Equal multipliers there still do
not remove the explicit nonzero defect(5).

The seventh-root cycle of x^2-2 has eta=-1 and multiplier8; its reversed
image at c=-7/4 has eta=0 and multiplier1. This gives an exact hostile to
transporting stability along a bijection of three states. The other
three-cycle at c=-2 has eta=1 and multiplier-8; its reversal has multiplier19.
The particular chosen cycle matters, not just its period.

## 3. Why63 is exceptional without supplying a return

For c=a/b reduced and b>1, the critical orbit starting at0 has nth
denominator b^(2^(n-1)), without cancellation. This follows inductively
because the new numerator is the square of a numerator coprime to b plus
a term divisible by b. Thus these critical orbits cannot be periodic.

Nevertheless, at c=-7/4,

    c(c+1)^2=-63/64,
    f_c^3(0)=c[c(c+1)^2+1]=-7/256.

The third numerator has no new prime. The stronger arithmetic explanation
is inherited from [row-braid typing, section5b](collatz_mod6_20260917_row_braid_typing.md):
the cubic norm a(a+b)^2+b^3 is1 at(a,b)=(-7,4), equivalently
rho^-14=-7+4rho^2 for rho^3=rho+1. Its cited certified Thue computation is
not rerun or promoted to an elementary proof here.

For the separate affine system u->2u+1, u_6=63=3^2*7 has only primes seen
earlier. But ord_9(2)=6, whereas ord_3(2)=2: multiplicity at an old prime
carries genuinely new order. The distinct objects are:

* finite prime support, which fails to see this increased3-adic precision;
* the expanding denominator of the quadratic critical orbit;
* the multiplier and return thickness of a different, periodic orbit.

The preserved predicate of the63 bridge is one numerator cancellation.
It does not identify the dynamics, an actual Collatz source, or an orbit
of its terminal certificate. This is why the cancellation is useful for
designing an arithmetic test but cannot itself prove convergence.

The final charged grammar in the [companion Collatz packet](collatz_charged_terminal_family_20261012a.md)
leaves127=2*63+1 as its first residual. A separate exact check supplies its
first band exit: word(1,1,1,1,1,1,2,4,3,3) ends at29 and has odd charge8.
This individual trajectory is not unknown;29 is already a supplied seed.
The unresolved family obligation is a parameter-wide paid rule. The small
Mersenne value points back to the first-reset-2 branch, while an isolated
word at127 does not discharge that variable-depth branch.

## 4. An analogous finite support budget for colour gadgets

The [boundary palette interface](boundary_palette_compression_20261011.md)
and [weighted context interface](boundary_context_weights_20261011d.md)
already prove complete equality-partition gluing. The following consequence
specifies when iterative elimination can even display partial progress.

For an unprecoloured graph with two named ports at palette6, global palette
permutations act transitively on the six equal assignments and separately
on the30 separating assignments. Thus a complete two-port relation contains
all or none of each orbit. Its only possibilities are empty, E, D, and E+D.
An equality-preserving chain cannot gradually remove some separating colour
names: at this interface its only strict successful change is E+D to E.
Extension counts may change in the meantime; positive counts are not zero.

Retaining more named ports exposes intermediate states. With six named
ports there are203 equality partitions. For a specified pair u,v, exactly
52 identify u,v (contract that pair to one of five ports), and151 separate
them. These are complete palette6 universe counts, not realizability claims.

**Conditional finite completion lemma.** Fix the same geometric boundary B,
a protected equality partition pi*, and an initial fragment. Suppose every
nonterminal stage admits a finite unit-distance fragment, compatible at B
with disjoint interiors, which preserves pi* and excludes at least one
currently surviving separating partition. Then gluing reaches a nonempty
equal-port relation after at most151 strict stages when |B|=6.

**Proof.** Gluing intersects complete supports. The finite set of separating
partitions strictly decreases, while pi* retains a colouring lift. QED.
If the distinguished ports have unit distance, adding their omitted edge
gives a non-six-colourable unit subgraph. Every fragment must use distinct
valid coordinates, and their interiors really must be disjoint; alternatively
all coincidences and cross constraints must be handled in the actual union.

The elimination premise is **OPEN**. This lemma does not manufacture the
required fragments or reduce the known existence problem to a routine
search. It tells us which finite support to retain if a candidate family
really permits incremental elimination. Three or four ports are already
more informative than two; six is a concrete test interface, not a claimed
optimal width. Any failed elimination stage must retain its surviving
partition as the next hostile, rather than discard it by compression.

## Reproduction and scope

    python 04-computation/experiments/marked_return_thickness_20261012a.py
    python -O 04-computation/experiments/marked_return_thickness_20261012a.py

The [script](../../04-computation/experiments/marked_return_thickness_20261012a.py)
checks identities symbolically over rational polynomial rings, all1..16
repetitions at four specified cycles, critical denominators through8,
the literal63 identity, the residual127 first-band word, and complete partition universes through six ports
plus the seven-port palette cutoff. The [output](marked_return_thickness_20261012a.out)
records the exact scope. The all-depth statements follow from the displayed
proofs; finite repetition checks are positive and hostile controls.
