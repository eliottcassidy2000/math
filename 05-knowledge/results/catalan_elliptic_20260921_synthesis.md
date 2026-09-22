# Catalan, signed Collatz, divisor statistics, and the fruit elliptic curve

**2026-09-21. Status: PROVED scoped statements; CITED classical inputs;
FINITE-EXACT computations; VERIFIED twelve concrete Lean certificates.**
Positive Collatz convergence and classification of all signed cycles remain
**OPEN**. No novelty claim is made. The ambiguous asterisks in the pasted
formulas are read as multiplication, pending clarification; results using
that reading are explicitly conditional on the notation.

This session yields three useful connections: the squarefree constant
controls the sign of the proper-divisor defect; the fruit cubic has a real
`C2 x S3` symmetry whose central involution crosses the positivity boundary;
and Catalan's unit-gap phenomenon identifies only a small part of the exact
Collatz cycle gate. The ordered carry is the coordinate that the last
analogy must retain. The pasted large integers also need a decimal repair.

## Inheritance and the live concept board

Closest proved mechanisms are the ordered affine cycle gate in
[signed parameter strata](arithmetic_braids2_20260917_signed_cycles.md),
the divisor exponent-box classification in
[divisor balance, DB1--DB3](arithmetic_braids_20260917_divisors.md), and
the order-six lift in [arithmetic seams](arithmetic_seams_20260921_dynamics.md),
which routes to [THM-4146, order-six lift and fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md).
The immediate prior [mod-six synthesis](collatz_mod6_20260917_synthesis.md)
already contains the unit-gap classification and a stronger bounded-period
record than the independent replay here. These are inherited dependencies,
not new discoveries.

Canonical hostile examples are the negative seven-cycle with raw gap -139,
`p^2` at the proper-divisor boundary, and six vectors becoming only three
projective points. The corrected near miss is assuming a quotient retains
integrality, positivity, or sampling information. The least-used relevant
coordinates here are the positive real chamber, the repetition-invariant
reduced denominator, and the Galois action on two-torsion.

Anchor: test exactly what Catalan contributes to Collatz. Niche: derive a
spatial sign law from F,S,U. Wildcard: identify the actual elliptic symmetry
behind the supplied numbers and compare the two cubic curves.

| Live concept | Exact question | Retained coordinate / decisive test |
|---|---|---|
| Signed Collatz | What does negation conjugate? | Additive parameter and sign sheet; test +/-3 |
| Catalan gap | When does an exponent word realize integers? | Ordered carry B beside `2^K-3^L` |
| F,S,U | Which measure produces each logarithm? | Uniform integers versus an orbit; endpoints1 and primes |
| Fruit symmetry | Does the involution preserve positive solutions? | Positive chamber; threshold sqrt(5) |
| Polynomial curve | Is it the same elliptic curve? | j-invariant and point count at11 |
| Formal evidence | What did the kernel actually check? | Twelve closed arithmetic witnesses, no global classification |

## 1. Catalan gives a classification, but the cycle gate is wider

For odd n and odd nonzero b, away from `3n+b=0`, let

```text
U_b(n)=(3n+b)/2^v2(3n+b).
```

Negation gives the forward conjugacy `U_b(-n)=-U_(-b)(n)`.
For b=1, write `n=s m`, with positive odd m. The signed system is exactly
`(s,m) -> (s,U_s(m))`. Identifying the sheets loses a necessary coordinate:
`abs(U_1(3))=5` but `abs(U_1(-3))=1`. Negation does not reverse arrows.

An ordered word `(k_1,...,k_L)`, with partial sums K_i and total K, has

```text
B=sum_(i=0)^(L-1) 3^(L-1-i)2^K_i,
Delta=2^K-3^L,
n=bB/Delta,
q=abs(Delta)/gcd(B,abs(Delta)).
```

**PROVED, inherited and rederived:** it is an exact signed integer cycle
word iff `q|b`. Rotating the carry proves that all intermediate nodes are
odd integers and that the specified exponents are their actual valuations.
For b=+/-1 the criterion is `Delta|B`, not `abs(Delta)=1`.

The complete unit-gap list, including exponent-one boundaries, is

```text
abs(2^K-3^L)=1, K>=L>=1
iff (K,L)=(1,1),(2,1),(3,2).
```

It produces the fixed cycles `(-b)`, `(b)` and the two-cycle `(-5b,-7b)`.
This fixed-base special case has an elementary mod-eight/factorization
proof; the general Catalan theorem supplies context rather than a needed
deep dependency. The known negative seven-cycle at b=1 instead has

```text
-17,-25,-37,-55,-41,-61,-91,
(K,L)=(11,7), Delta=-139, B=2363=17*139, q=1.
```

Moreover word repetition multiplies B and Delta by the same geometric
factor, leaving q unchanged. Even the fixed point 1 has raw gaps 1 and 7
under words `(2)` and `(2,2)`. Large raw gaps cannot be treated as a
standalone obstruction to integrality.

**FINITE-EXACT independent replay:** every signed b=+/-1 cycle has
`L<=K<=2L`. Enumerating all 250,952 words for `L<=10` therefore covers all
heights, and finds the four familiar signed cycles for each parameter.
This is an independent certificate, not an improvement of the inherited
period 23 record and not a proof excluding longer cycles or divergence.

The actual global descent obligation remains: for every positive odd n>1,
some prefix of its **actual** word must satisfy
`(2^K-3^L)n>B`. Catalan supplies neither that prefix nor its length.
Full proofs and source boundaries: [Catalan lane](catalan_elliptic_20260921_catalan.md).

## 2. A precise squarefree sign law emerges from F,S,U

Use proper nontrivial divisors throughout:

```text
F(n)=#{d|n:1<d<n},
S(n)=#{d|n:1<d<n, d squarefree},
U(n)=#{p|n:p prime, p<n},
D(n)=F(n)-S(n)-U(n).
```

For n uniform in 1,...,X, put c=6/pi^2. Exact divisor incidence and
elementary convolution, with classical prime-harmonic estimates, give

```text
P_X(n prime) ~ 1/log X,
E_X F = log X+O(1),
E_X S = c log X+O(1),
E_X U = log log X+O(1),
P_X(n squarefree) -> c.
```

The shared constant has an explicit source:
`sum mu(n)^2/n^s=zeta(s)/zeta(2s)`, whereas
`sum 2^omega(n)/n^s=zeta(s)^2/zeta(2s)`. Summing divisor incidences adds
the extra zeta factor. This is not an assumption of random orbit behavior.

**PROVED sign law:**

```text
P_X(D<0) -> c,    P_X(D=0) -> 0,    P_X(D>0) -> 1-c,
E_X D = (1-c)log X-log log X+O(1).
```

Mechanism: a squarefree composite has D=-omega(n). A nonsquarefree number
with r>=4 distinct prime factors has
`D>=2^(r-1)-r-1>0`. The exceptions lie in `omega(n)<=3`, a density-zero
set proved by a finite-prime CRT estimate. Thus the limiting majority
has negative D while the mean is eventually positive. The observable is
unbounded, so there is no contradiction.

The exact zero set is `1,p,p^3,p^2qr` with p,q,r distinct. Conditioning on
the odd residue rows gives negative-sign densities `9/pi^2,6/pi^2,9/pi^2`
for rows1,3,5 mod6. These are spatial densities; a Collatz orbit theorem
still needs an orbit-specific bridge. At X=10^6, the exact negative/zero/
positive counts are 573459/165862/260679, with independent divisor checks
through 10^4. The finite zero proportion illustrates the danger of inferring
the limiting law from a modest sample.
Proofs, constants and controls: [statistics lane](catalan_elliptic_20260921_statistics.md).

## 3. The large numbers reveal an actual central involution

The numbers identify the fruit equation

```text
a/(b+c)+b/(a+c)+c/(a+b)=4.
```

As pasted, their left side is approximately 2.318469382009154. Leaving aa
unchanged and dividing bb and cc by 10 yields a primitive positive integer
solution with the stated digit lengths 81,80,79. This exact repair is the
image of 9G, where G=(-4,28) on
`E_4:y^2=x^3+109x^2+224x`. The identification is checked by exact rational
arithmetic and the cleared identity is kernel-certified. The source is
[Bremner--Macleod's 2014 paper](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf),
with its 2025 corrigendum noted in the detailed lane; no table-wide
minimality claim is inherited.

The projective linear coordinate map is

```text
[a:b:c] -> [-28(a+b+2c):364(a-b):6(a+b)-c]=[X:Y:Z],
[X:Y:Z] -> [56Z-X+Y:56Z-X-Y:-56Z-12X].
```

Put u=b+c,v=c+a,w=a+b. The fruit equation becomes
`(u+v+w)(1/u+1/v+1/w)=14`. Reciprocation of u,v,w therefore supplies
an involution, whose a-coordinate in fruit variables is

```text
J_a=-a^2+b^2+c^2+ab+ac+bc
```

and similarly for b,c. On E_4 this is translation by the rational
two-torsion point (0,0). Coordinate permutations supply S3; the central
involution commutes with them and supplies an independent C2. This is
an actual `C2 x S3` action on the curve.

**PROVED sharp positivity boundary:** for every positive real triple
whose fruit sum is at least sqrt(5), J produces exactly one negative
coordinate. Positive rational triples with positive J-images approach
sqrt(5) from below. Thus at sum 4 the central symmetry inevitably leaves
the positive chamber. The proof retains `r=max(a,b,c)/(sum of the other
two)` and their normalized product q, rather than identifying signed and
positive solutions.

Sharpness has an explicit integer witness family: for odd n>=3, take
`p=F_(n+1),q=F_n` and `(a,b,c)=(2pq,1,2q^2-1)`. Cassini gives
`q^2+pq-p^2=1`, hence `J_a=2q^2+1>0`. The fruit sums approach sqrt(5)
from below. These varying sums do not give additional solutions at sum 4.

The six rational torsion points carry the same dihedral finite action as
the previous vector hexagon: rotation corresponds to translation by an
order-six point, reflection to negation. This is an explicit equivariant
bijection of finite sets. All six torsion points have a zero fruit
denominator, however. The finite group bridge does not transfer valid
fruit solutions, heights, or scalar quadratic dynamics.
Proofs and 72 exact translated-point controls:
[elliptic lane](catalan_elliptic_20260921_elliptic.md).

## 4. The polynomial yields a second curve, separated by an exact invariant

Under the stated multiplication reading of b,c,s,t, expansion gives

```text
c=(x+1)^2,
s=x(2x^2+4x+1),
t=-(2x+1)(x^2+x-1),
s+t=c,
(2x^3+4x^2+1)-s=1-x.
```

If one asks the additional square-value question `y^2=2x^3+4x^2+1`,
the change X=2x,Y=2y gives `E_0:Y^2=X^3+4X^2+4`. It has the nontrivial
integer point (x,y)=(10,49), obtained from 4(0,2) in the elliptic group.
Its j-invariant differs from E_4, and at the common good prime 11 the
two curves have 16 and 12 points respectively. Consequently they are not
Q-isogenous; in particular no nonconstant Q-morphism connects the two
smooth projective curves. A different kind of correspondence would
require an explicit construction.

There is another genuine S3 here: the Galois action on the three nonzero
two-torsion points of E_0. It is not the coordinate-permutation S3 of the
fruit curve. Keeping the acting group, acted-on set, and preserved
predicate together prevents a false identification.

Every integral square-value solution has even x. With x=2z,y=2k+1,
the equation becomes `k(k+1)=4z^2(z+1)`. If z>0 and z+1 is a square,
the original value is one above a nonzero square and is not a square.
The scan -2<=x<=200000 finds only (-2,1),(0,1),(10,49) with y>=0;
this is **FINITE-EXACT**, not a global integral-point classification.
Full proof: [polynomial lane](catalan_elliptic_20260921_polynomial.md).

## 5. What the synthesis preserves, and what remains to prove

| Source -> target map | Preserved predicate | Information lost | Necessary sidecar / cheapest test |
|---|---|---|---|
| Ordered Collatz word -> power gap | Affine multiplier | Translation and integrality | Carry B; compare (1,3) and(2,2) |
| Signed integer -> magnitude | Absolute size | Parameter sheet | s; compare +/-3 |
| Divisor exponents -> sign(D) | Squarefreeness outside sparse exceptions | Defect magnitude and time order | Sampling measure; compare mean and sign law |
| Fruit point -> reciprocal pair sums | Cubic equation | Positivity | Real chamber; J threshold sqrt(5) |
| Six vectors -> six torsion points | Dihedral finite action | Ambient dynamics and admissibility | Pair-sum denominators; all six vanish somewhere |
| Cubic square value -> elliptic group point | Rational equation | Affine integrality | Denominators at 5(0,2) |

The common lesson has a concrete form in each lane: a small exact model
can generate infinitely many states, but the arithmetic predicate needed
for the original problem must be carried along with that generation.
For Collatz the missing theorem is universal descent along the actual
guarded word. For the fruit example the extra predicate is positivity,
achieved here at 9G and infinitely often; global leastness is not certified.
For the polynomial the open obligation is global integral-point control.
No result here identifies these different obligations or resolves the
Collatz one.

## Reproduction and formalization scope

Run from any working directory, using the repository-relative path here:

```text
python 04-computation/experiments/catalan_elliptic_20260921_verify.py
```

The [master replay](../../04-computation/experiments/catalan_elliptic_20260921_verify.py)
runs each lane normally and under Python optimization, compares decoded
JSON, rechecks source hashes, freezes proof-note inputs, and performs a
clean Lean build through its public root. The
[manifest](catalan_elliptic_20260921_manifest.json) records actual outcomes
and finite universes. Peer review covered the Catalan gates, statistics
proofs, elliptic transformations and polynomial distinctions independently.

[CatalanEllipticAudit](../../04-computation/lean/CatalanEllipticAudit/README.md)
contains twelve closed arithmetic certificates in Lean 4.30.0, all axiom-free.
They verify the literal/repaired fruit equalities and denominator facts,
the point (10,49), and the negative seven-cycle carry and cleared steps.
They do not formalize the general density theorem, elliptic symmetry,
bounded-period exhaustiveness, Catalan's theorem, or Collatz convergence.
