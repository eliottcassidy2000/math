# Ternary normalization on primitive Pythagorean triples

**Date:** 2026-09-25. **Status:** PROVED elementary, scoped dynamical
isomorphisms and contraction; FINITE-EXACT independent controls; CITED
prime-search/certification status. Positive integer Collatz convergence remains
**OPEN**. No priority claim is made for the rational Collatz extension,
Pythagorean parametrization, or 3-adic repunit odometer.

## Inheritance, portfolio, and concepts

The nearest proved carriers are [THM-3333, Gaussian-square light cone](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
[THM-3756, odd-square ordinals and affine Berggren descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md),
and the [Collatz-to-Berggren edge transport](collatz_mod6_20260921_berggren_edge_transport.md).
The latter already identifies inverse fibres with sparse samples of parabolic
Berggren rays; its multiplier-changing children rule out a naive step-for-step
identification. The [inverse-fibre braid](arithmetic_braids_20260917_collatz.md)
already proves the sibling exponent's 3-adic isometry.

Canonical hostile: a primitive triple becomes nonprimitive under a raw parameter
operation unless content is retained (MISTAKE-418). Corrected near miss:
a prime-factor event does not imply a Collatz halving event. Least-used sidecar:
the **sign of the even leg**, which distinguishes reciprocal rational states.
Anchor: an actual Collatz-related dynamical system on primitive triples. Niche:
repunit clocks on residue-refinement trees. Wildcard: test the pasted decimal
prime sequence for a real carry-controlled descent, independently of primality.

| Concept | Exact retained structure | Boundary |
|---|---|---|
| Primitive triple | both signed square gaps C+B and C-B | unsigned triangles identify reciprocal states |
| Rational Collatz | numerator, odd denominator, exact halving count | integer boundary has no denominator factor to cancel |
| Ternary normalization | gcd of new numerator and denominator is 1 or 3 | prime-to-3 denominator is invariant |
| Geometric rank | hypotenuse contracts when 3 divides denominator | cannot extend the inequality to denominator 1 |
| Repunit clock | isometry of all finite ternary residue trees | ordinary integer image is sparse |
| Prime event | actual divisibility order and digit carry | no theorem of a universal geometric capacity horizon |

Research moves used: type every analogy and implication; recover the missing
second coordinate; compute a quotient's lost data before using its rank.
No new META-PATTERNS promotion is needed: the evidence instantiates those cards.

## 1. A lossless marked-triple chart

Let Q_odd be positive rational numbers s/t in lowest terms with s,t positive odd
integers. Define

    Phi(s/t)=(A,B,C)=(st,(s²-t²)/2,(s²+t²)/2).             (1)

Here A is the positive odd leg, B is an **oriented even leg**, and C is the
positive hypotenuse. A negative B records a marker on the ordinary triangle
with sides (A,|B|,C). Add the single degenerate terminal (1,0,1), representing 1.
The resulting marked primitive triples are exactly Phi(Q_odd).

The inverse is intrinsic:

    s=sqrt(C+B),       t=sqrt(C-B),
    s/t=(C+B)/A.                                        (2)

Proof: `(C+B)(C-B)=A²`; the two factors are positive, odd, and coprime in a
primitive triple, so each is a square. Conversely, coprime odd s,t make (1)
primitive and satisfy A²+B²=C². B's sign is essential: exchanging s,t changes
B to -B while replacing the rational state by its reciprocal.

This is the inherited odd-root chart, now used as a faithful state space for
a dynamical system. The new useful property comes from primitive content.

## 2. The exact content splitter

On positive odd rationals define the plus accelerated map by

    U_+(s/t)=(3s+t)/(2^k t),       k=v2(3s+t).

The numerator is divided by its entire power of 2. In reducing the remaining
fraction, the only possible common odd factor is

    g=gcd((3s+t)/2^k,t)=gcd(3,t) in {1,3}.                (3)

Indeed gcd(s,t)=1 and t is odd. Thus the next reduced pair is

    s'=(3s+t)/(2^k g),       t'=t/g.                     (4)

The exact triple operation F_+=Phi U_+ Phi^(-1) therefore exists and is
constructive. The [triple lane](ternary_triples_20260925.md) proves the signed
version, inverse guards, and all side conditions independently.

Write t=3^r q with 3 not dividing q. Equation (4) proves two statements:

* q is unchanged on every step.
* r decreases by exactly one while positive, then stays zero.

This is a genuine ternary normalization law, not an inference from how often
primes occur. For the subset

    P_r={marked primitive triples : C-B=9^r},

we obtain a graph isomorphism with odd rational Collatz states whose reduced
denominator is 3^r, and exact transitions

    P_r --> P_(r-1) for r>0;       P_0 --> P_0.            (5)

Every input in P_r reaches P_0 after exactly r accelerated steps. P_0 consists
of `(n,(n²-1)/2,(n²+1)/2)` for positive odd integers n, so its induced dynamics
is precisely integer odd Collatz. The exponent r is read from the marked
triple itself: r=v3(C-B)/2 on these strata.

## 3. Hypotenuse contraction: the subtle extra property

When 3 divides t, equation (4) simplifies to

    s'=(s+t/3)/2^k,       t'=t/3,       k>=1.

Consequently

    2C' <= (s+t/3)²/4+t²/9.

Compare this with 2C/3=(s²+t²)/3. The difference is

    [3(s-t)²+4t²]/36 > 0.

Hence, on every positive plus step with r>0,

    C' < C/3,       C'-B'=(C-B)/9.                       (6)

The hypotenuse is a strictly decreasing positive integer rank until the
integer stratum is reached. This proves finite entry; it does not assume
Collatz convergence. The contraction comes from **canceling a common factor
3 in the two spinor coordinates**, which divides their quadratic lift by 9.

Two exact trajectories distinguish the conclusions:

    (117,44,125) -> (3,-4,5) -> (1,0,1),
       13/9             1/3          1;

    (39,80,89) -> (7,24,25) -> (11,60,61),
       13/3             7            11.

In the second example, C falls from 89 to 25 during ternary cancellation,
then rises to 61 on the integer stratum. The missing contraction at t=1 is a
specific arithmetic boundary. Therefore convergence on the union of the
P_r is equivalent to integer Collatz convergence, not already a consequence
of the off-boundary rank.

## 4. An intrinsic triple update without extracting square roots

The parameter calculation can be compiled into arithmetic on A,B,C alone.
For sigma=+1, or for a positive legal sigma=-1 step, put

    R=10C+8B+6 sigma A=(3s+sigma t)²,
    k=v2(R)/2,    K=2^k,
    g=3 if 9 divides C-B, otherwise g=1.

Then

    A'=[3A+sigma(C-B)]/(K g²),
    B'=[6 sigma A+(8+K²)B+(10-K²)C]/(2K² g²),
    C'=[6 sigma A+(8-K²)B+(10+K²)C]/(2K² g²).             (7)

All divisions are exact on the stated domain. Expansion of `(3s+sigma t)²`
proves the formulas, and equation (3) proves primitivity of the output.
The coefficient 10 is 3²+1 and the coefficient 6 is the cross term; these
are algebraic coefficients, not consequences of triangular-number capacity.

For minus, positivity has an intrinsic linear guard:

    3s>t  iff  4C+5B>0.                                 (8)

On a legal minus step with 3 dividing t, the sharper inequality C'<C/4 holds.
But the domain is not invariant: the positive primitive triple

    (783,56,785) -> (45,-28,53) -> (3,-4,5)

would next give rational numerator zero. The verifier rejects that next
step. The sign and positive-domain guard must not be erased.

## 5. Why marked triples and the power-nine gap restriction are necessary

For general t, its 3-free part q cannot disappear. In particular the full
marked-triple system has noninteger fixed points

    x=1/(2^k-3),       k>=3.

At k=3, Phi(1/5)=(5,-12,13) is fixed. Yet Phi(5)=(5,12,13) goes to the terminal
(1,0,1). These are the same unmarked 5-12-13 triangle with different markings.
Even a starting triple with positive B can reach the noninteger fixed point:

    (45,28,53) -> (5,-12,13) -> (5,-12,13).

Thus neither forgetting the orientation nor extending P_r to every primitive
triple preserves the desired integer-root statement.

The [Berggren lane](ternary_berggren_20260925.md) supplies two further boundaries.
A fixed nonempty accelerated Collatz word has parameter matrix with diagonal
entries 3^L and 2^K, so its invariant

    tr(M)²/det(M)=(3^L+2^K)²/(3^L 2^K)

is noninteger. Every Berggren word has determinant +/-1 and integer trace,
hence integer invariant. No fixed projective coordinate change turns the
former into the latter. This does not exclude nonlinear encodings or adaptive
word clocks. The same lane also proves that consecutive nondegenerate plus
**edge-triangles** are never Berggren ancestor/descendant, at any distance;
those edge-triangles are a different representation from the state chart (1).

## 6. The repunit connection that really preserves ternary recursion

For an integer base b>1 with b=1 mod3, define R_b(j)=(b^j-1)/(b-1), j>=0. For j!=i,

    v3(R_b(j)-R_b(i))=v3(j-i).                            (9)

Factor out b^i and use `v3(b^h-1)=v3(b-1)+v3(h)`. Thus, modulo every 3^d,
R_b is a permutation of the exponent classes, compatible with reduction to
the preceding level. It extends to an isometric bijection of Z_3, conjugating

    j -> j+1       to       z -> bz+1.                   (10)

This is an isomorphism of rooted ternary **residue-refinement trees**. It is
not a statement that j->R_b(j) is onto the ordinary nonnegative integers.
Its inverse is constructive: choose the unique next ternary digit among three
lifts at each level.

In base 10, the pasted family is simply

    N_k=(10^k-7)/3=3 R_10(k)-2.

In the Collatz inverse fibre of odd y not divisible by 3,

    x_j=(2^(k0+2j)y-1)/3,

choose the least legal k0>=2 with x_0>y. The inherited Berggren parabolic
operation on roots is B1(s,t)=(s+2t,t), and

    (x_j,y)=B1^H(j)(x_0,y),
    H(j)=2^(k0-1) R_4(j).                                (11)

The factor 2^(k0-1) is a 3-adic unit, so H is also an isometric bijection on
each finite ternary level. Composing inverse R_10 with H gives an explicit
isomorphism between the decimal repunit-value residue tree and the Collatz
inverse-fibre height-residue tree on a primitive Pythagorean ray. The former
is itself an exponent-residue tree after using R_10 as its coordinate map. The
[Berggren lane](ternary_berggren_20260925.md) proves the height decoder and
its exact, nonuniform ordinary tree clock.

There is even a precise first nonlinear carry:

    R_4(j) = j+3 binom(j,2)+9 binom(j,3)+...,
    R_10(j)= j+9 binom(j,2)+81 binom(j,3)+... .            (12)

The triangular coefficient records when each radix departs from ordinary
exponent counting. Prime or composite values obey exactly the same identities.
For example R_4(j)=0,1,5,21,... omits height 2 over ordinary integers, despite
covering every residue class modulo every power of 3. Arithmetic height is the
sidecar that forbids turning this local isomorphism into global coverage.

## 7. What survives from the pasted prime examples

The [digit audit](ternary_digits_20260925.md) proves and exactly checks:

* N_2 through N_8 are seven primes; N_9=333333331=17*19607843 is the first
  composite in that run. Last-prime digit length 8 is different from eight
  prime steps starting at 31.
* 17 divides N_k exactly when k=9 mod16, because ord_17(10)=16. This is a
  recurring modular phase, not a one-time breakdown of the recursion.
* Every N_k with k>=2 first grows under odd Collatz. For every k>=4 it
  descends after two odd steps, independently of primality. At k>=5 the
  exact endpoint is `3*10^k/16-1`; k=4 has an extra halving.
* The cubic values 7,34,93,196,355 and F_5=641*6700417 are correct, but the
  first twelve primes sum to **197**, not 196.

The prime-related digit lengths also need their quantifiers. A multidigit
prime fixed by rotation is a repunit prime, whose length must be prime.
A nonrepunit prime with d digits has exactly d distinct rotations, but d
need not be prime: 1193 and 193939 are circular examples of lengths 4 and 6.
Circular primes and primes surviving every permutation are different notions.
The one-digit primes are fixed by rotation as well.

Finiteness of nonrepunit circular primes is conjectural; the research catalogue
records known lengths 1--6 and finite exclusions, not a theorem forbidding all
larger lengths ([De Geest's search report](https://www.worldofnumbers.com/circular.htm)).
Repunit infinitude is also open. Large certified examples go beyond the five
lengths in the inherited bounded search through 1100: Enge and Underwood
certified R_109297 and provide a checkable certificate
([Enge's report](https://enge.math.u-bordeaux.fr/blog/ecpp-109297.html)).

A repeated factor does not supply an upper bound on an increasing Collatz
prefix: the composite integers 2^(2m)-1 have arbitrarily long initial odd-step
growth. The ternary mechanisms that survived the probes are primitive-content
cancellation and compatible radix carries; neither needs a prime horizon.

## 8. Connection contract and reproduction

Source: reduced positive odd rationals with denominator a power of 3.
Target: marked primitive triples with C-B a power of 9, plus the degenerate
terminal. Map and inverse: (1)--(2). Preserved predicate: exact positive
Collatz edges and integer-boundary entry, with a proved hypotenuse rank before
entry. Lost by forgetting markers: reciprocal direction, denominator, and
root behavior. Sidecars: sign, orientation, exact halving count, primitive
content, and integer boundary. Cheapest decisive tests: 13/9's orientation
crossing, 13/3's contraction followed by integer growth, and the fixed 1/5.

A separate source/target pair uses residue classes of digit exponents and
Berggren ray heights. Equations (9)--(11) preserve all ternary prefix relations;
they destroy ordinary height density. No isomorphism of the full integer
Collatz graph with the ordinary Berggren tree has been claimed.

```
python3 04-computation/experiments/ternary_bridge_20260925.py
python3 -O 04-computation/experiments/ternary_bridge_20260925.py
```

The [output](ternary_bridge_20260925.out) records 24255 legal signed comparisons
between the intrinsic formula and independent spinor normalization, 6099 strict
contractions, 3240 exact boundary-entry cases, both finite ternary permutations
and their inverses through level 8, 189 literal ray samples (including the root target y=1), and the signed
hostiles. Companion lane proofs and independent audits are linked above.
Universal convergence on the integer stratum remains OPEN.
