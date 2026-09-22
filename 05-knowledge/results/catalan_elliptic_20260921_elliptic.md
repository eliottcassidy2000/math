# Fruit cubic: exact input audit, reciprocal symmetry, and the positivity barrier

**2026-09-21. Status: PROVED algebraic statements; CITED elliptic-curve
background; FINITE-EXACT arithmetic controls.** No novelty or global
minimality claim. The rational rank is not computed by this session.

The pasted integers do not solve the fruit equation. Dividing the second
and third integers by ten, and leaving the first unchanged, does solve it.
The corrected triple is exactly the image of `9G`, where `G=(-4,28)` on
`E: y^2=x^3+109x^2+224x`. More structurally, reciprocating the three pair
sums produces the curve's rational two-torsion translation. Together with
coordinate permutations this gives an actual `C2 x S3` action. Its central
involution cannot preserve positivity once the fruit sum reaches `sqrt(5)`.

## Inheritance, source boundaries, and the live board

The closest repository mechanism is the six-vector dihedral action in
[the arithmetic-seams dynamics note](arithmetic_seams_20260921_dynamics.md),
which inherits [THM-4146, order-six lift and fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md).
The hostile there is mistaking a finite-set lift for a scalar quadratic
six-cycle. The corrected near miss here is mistaking a signed solution
for a positive one. The least-used sidecar is the positive real locus and
the six points where the original rational expression has zero denominators.
The relevant research move is to retain action, positivity, and basepoint
when passing through a quotient or a change of coordinates.

| Lane | Object / operation | Retained predicate / omitted information |
|---|---|---|
| Anchor | Exact integer input and rational fruit sum | Equality; decimal transcription cannot be inferred |
| Niche | Symmetric cubic and pair-sum reciprocation | Cubic equation; positivity can fail |
| Wildcard | Six torsion points and the earlier hexagon | A finite group action; no identification of global dynamics |
| Height | Multiples of an infinite-order point | Rationality; first entrance into a positive arc is extra data |

Primary sources were checked directly. Bremner--Macleod provide the family
of elliptic curves, coordinate formulas, and the ninth-multiple example
in Section 2 and Remark 2.2 of
[their 2014 paper](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf).
The paper's 2025 corrigendum repairs several entries of its numerical table;
the N=4 digit counts remain 79,80,81. Consequently no unexamined table-wide
leastness claim is inherited from the original paper.
[Corrigendum record and DOI](https://publikacio.uni-eszterhazy.hu/8862/).

The standard background used below is the group law, prime-to-residue-
characteristic torsion reduction, real elliptic-curve topology, and
quadratic canonical height. These are explained in the corresponding
chapters of [Elkies's Harvard 2024 course](https://people.math.harvard.edu/~elkies/M223.24/index.html).
The proofs and exact arithmetic below are independently supplied; the
finite computation is not used to infer rank or global leastness.

## 1. Literal input and the exact decimal repair

Write

```text
A = 154476802108746166441951315019919837485664325669565431700026634898253202035277999
B = 368751317941299998271978115652254748254929799689719709962831374716372246340555790
C = 43736126779286972578612526023713901528165375581616136186214379933784234677720360
```

For `s=a+b+c`, define

```text
F_N(a,b,c) = s^3 -(N+2)s(ab+ac+bc)+(N+3)abc.
```

Multiplying the fruit equation by `(a+b)(a+c)(b+c)` gives exactly `F_N=0`.
Thus the equivalence requires all three pair sums to be nonzero.

**FINITE-EXACT.** `F_4(A,B,C)` is nonzero (negative); its full integer value
and the exact rational fruit sum are in the companion JSON. The fruit sum
is approximately `2.318469382009154`. In contrast,

```text
(a,b,c)=(A,B/10,C/10)
```

is a primitive positive integer triple with `F_4=0`; all pair sums are
positive. Its coordinate digit counts are `(81,80,79)`. This is an explicit
repair of the input, not a claim that the pasted triple already works.

## 2. The cubic-to-elliptic correspondence is projective linear

**PROVED.** Set `u=N+3`, `v=2N+5`. The linear map

```text
X = -4u(a+b+2c),   Y = 4uv(a-b),   Z = (N+2)(a+b)-c
```

satisfies

```text
Y^2 Z-X^3-(4N^2+12N-3)X^2 Z-32(N+3)X Z^2
    =64u^2 v^2 F_N(a,b,c).                         (1)
```

Its inverse, projectively, is

```text
[a:b:c]=[8uZ-X+Y : 8uZ-X-Y : -8uZ-2(N+2)X].        (2)
```

The inverse matrix times the forward matrix is `8uv I`. Hence for
`N != -3,-5/2` these are inverse projective linear transformations, not
just formulas on a dense affine chart. The additional singular parameter
`N=3/2` still has an invertible linear map, but the two cubics are singular
and must not be called elliptic curves. Their discriminant is
`2^14(N+3)^2(2N-3)(2N+5)^3`.

For `N=4`, (1)--(2) read

```text
[X:Y:Z]=[-28(a+b+2c):364(a-b):6(a+b)-c],
[a:b:c]=[56Z-X+Y:56Z-X-Y:-56Z-12X].                 (3)
```

The composite scalar is `728`, and the residual multiplier is `529984`.
The point at infinity `O=[0:1:0]` corresponds to `[1:-1:0]`.
For an affine point `(x,y)` one can use the three entries
`(56-x+y,56-x-y,-56-12x)`, clearing denominators and normalizing their
common gcd. Polynomial identity (1) proves the result; the script also
compares every coefficient in `Z[a,b,c,N]` independently of point samples.

## 3. Pair-sum reciprocation supplies the missing central involution

Put `u=b+c,v=c+a,w=a+b`. On the nonzero-pair-sum locus the fruit equation is

```text
(u+v+w)(1/u+1/v+1/w)=2(N+3).                        (4)
```

It is invariant under simultaneous reciprocation of `u,v,w`. In projective
coordinates this is the Cremona map `[u:v:w] -> [vw:uw:uv]`. In fruit
coordinates it becomes the symmetric quadratic map

```text
J_a=-a^2+b^2+c^2+ab+ac+bc,
J_b= a^2-b^2+c^2+ab+ac+bc,
J_c= a^2+b^2-c^2+ab+ac+bc.                          (5)
```

**PROVED, for every N.** With `D=(a+b)(a+c)(b+c)`, expansion gives

```text
F_N(J(a,b,c))=8D F_N(a,b,c),
J(J(a,b,c))=8D(a,b,c).                              (6)
```

Thus `J` is an involution on the valid fruit locus and commutes with every
coordinate permutation. At its missing projective chart points, its
restriction to a smooth projective cubic extends uniquely as an
automorphism. This extension does not make the fruit expression itself
defined when a pair sum vanishes.

On E at N=4 the extension is translation by `U=(0,0)`:

```text
tau_U(x,y)=(224/x,-224y/x^2),  x != 0,              (7)
tau_U(O)=U, tau_U(U)=O.
```

One can derive (7) by drawing the line through `(0,0)` and `(x,y)` and
using the Weierstrass addition law. To verify the conjugacy, let `M` denote
the map (3), and compare `M(J(a,b,c))` to
`(224XZ,-224YZ,X^2)`. Their XY cross minor is zero; their XZ and YZ cross
minors are respectively

```text
-81536(a+b+2c)F_4,   -1059968(a-b)F_4.
```

This proves equality on the curve, with extension at the chart exceptions.

## 4. Actual C2 x S3 symmetry and the finite hexagon bridge

Let `T=(56,728)`, `Q=(4,52)`, and `U=(0,0)`. Direct addition gives

```text
2T=Q,  3T=U,  6T=O,
E(Q)_tors={O,U,(4,+/-52),(56,+/-728)}.              (8)
```

The displayed points have exact orders as asserted. Completeness of (8)
uses the standard torsion reduction theorem plus the independently
counted good reductions `#E(F_11)=12`, `#E(F_17)=18`. Every primary torsion
part away from 11 and 17 divides both counts; the 11-part is excluded by
reduction at 17 and the 17-part by reduction at 11. The resulting bound is
six, attained by T. No rank computation is involved.

The cyclic coordinate permutation `(a,b,c)->(b,c,a)` is translation by Q;
the transposition `(a,b,c)->(b,a,c)` is negation. For the latter this follows
immediately from (3). For the former, it is also the order-three
translation identified in the primary source; alternatively, translate
its image of O back to O and use the elementary automorphism classification
for an elliptic curve with `j != 0,1728`. Its image of O is Q, and an
order-three automorphism cannot have the negation linear part.

It follows that the actual curve automorphisms generated by permutations
and J form

```text
<tau_Q,[-1]> x <tau_U> = S3 x C2.
```

The product is direct: U is central of order two, and its translation is
not among `tau_(iQ)` or `tau_(iQ) composed with [-1]`. Nonzero translations
are automorphisms of the curve, not endomorphisms of the pointed elliptic
curve, because they do not fix O.

**PROVED finite action comparison.** Put `V=U+Q=(56,-728)=-T`, an order-six
point. The bijection from the earlier vector hexagon is

```text
B^k v <-> kV,       k in Z/6Z.
```

Choose the earlier reflection J0 to fix v. Then rotation B corresponds to
translation by V and reflection J0 to negation. The relations are the same:
rotation has order six, reflection has order two, and reflection reverses
rotation. The central half-turn corresponds to U; coordinate cycling
corresponds to `4V=Q`. This is an equivariant bijection of finite G-sets.
It is not a conjugacy between the quadratic polynomial and an elliptic
map, nor does it transfer arithmetic height or positivity.

There is a useful additional boundary: **all six torsion points in (8)
correspond to zero denominators in the fruit expression**. Indeed setting
`a+b=0` in `F_N` gives `c(c^2-a^2)`, and the analogous three lines yield
exactly those six projective points. Thus the torsion hexagon itself lies
outside the admissible fruit locus. A translated six-point orbit of an
infinite-order point, such as `{G+kT}`, gives a different, valid locus.

## 5. The positivity threshold is sqrt(5), sharply

This strengthening was proposed by the parent lane and independently
audited here. It applies to positive real a,b,c, not only rational points.

**PROVED.** If the fruit sum is at least `sqrt(5)`, then J sends a positive
triple to exactly one negative coordinate and two positive coordinates.
The threshold is sharp: there are positive rational triples, with positive
J-images, whose fruit sums approach `sqrt(5)` from below.

Proof. Relabel so `a=max(a,b,c)`. Put

```text
r=a/(b+c) >= 1/2,  q=bc/(b+c)^2 in (0,1/4].
N=r+(r+1-2q)/(r^2+r+q),
J_a/(b+c)^2=1+r-r^2-q.                             (9)
```

The other two J-coordinates are strictly positive because `a>=b,c>0`.
Suppose `J_a>=0`. If `r>=1`, (9) forces `r<phi=(1+sqrt(5))/2`; also
`N<r+1/r<phi+1/phi=sqrt(5)`. If `1/2<=r<=1`, the largest-coordinate
condition gives `q>=r(1-r)`. N is strictly decreasing in q, since its
q-derivative is `-(2r^2+3r+1)/(r^2+r+q)^2`. Hence

```text
N <= 2r-1/2+1/(2r) <= 2 < sqrt(5).
```

The second inequality follows from convexity and the endpoint values
`3/2,2`. This proves the contrapositive.

There is an explicit primitive integer family proving sharpness. Use
`F_0=0,F_1=1` for the Fibonacci numbers. For odd `n>=3`, put
`p=F_(n+1), q=F_n` (these p,q are new variables, distinct from the normalized
q in (9)). Cassini's identity gives

```text
q^2+pq-p^2=1.
```

Take

```text
(a,b,c)=(2pq,1,2q^2-1).                             (11)
```

This triple is primitive, and a is its largest coordinate. Directly from
(5),

```text
J_a=4q^2(q^2+pq-p^2)-2q^2+1=2q^2+1>0.
```

The other two J-coordinates are positive as above. In the notation of
(9), `r=p/q -> phi` and the normalized product is
`(2q^2-1)/(4q^4) -> 0`. Thus the fruit sums of (11) tend to
`phi+1/phi=sqrt(5)`, and the proved strict bound puts every one below it.
For example `n=3` gives `(12,1,7)`. These triples have varying fruit sums;
they are not asserted to solve N=4. This supplies sharpness without a
density argument and retains the equality and failure boundary.

In particular every positive solution at integer N>=3 has this positivity
loss. For N=4 there is an even simpler bound: `4<r+1/r` implies
`r>2+sqrt(3)`, so one fruit coordinate is much larger than the other two
combined. The image in (5) necessarily has that coordinate negative.

## 6. Signed points, positive arcs, and the ninth multiple

The small point `G=(-4,28)` gives `(11,4,-1)`. Its U-translate is
`(-56,-392)`, giving `(-5,9,11)`. Both solve the rational equation but fail
positivity. Their smallness is therefore no positive-integer shortcut.

**PROVED positivity test.** For a real affine point on E, its projective
fruit triple has a representative with all coordinates positive iff

```text
x < -14/3   and   x^2+112x+784 > 0.                (10)
```

For necessity, start with positive a,b,c in (3), or use the following
sign argument. In the inverse map the sum is `14(4-x)`. If `x>4`, a positive
projective triple would require all three inverse entries negative;
then their first two sum `112-2x<0` forces x>56, while their product
`-(x-4)(x^2+112x+784)` is negative, a contradiction. At x=4 the coordinate
sum is zero. Thus positivity forces x<4. Now the third entry is positive
iff x<-14/3. The first two have positive sum, and their product equals
`(4-x)(x^2+112x+784)`; both are positive iff the quadratic is positive.
This also proves sufficiency. In particular being on the negative-x
bounded real component alone is insufficient.

**FINITE-EXACT.** The companion program computes all `nG+kT` for
`1<=n<=12`, `0<=k<=5`, checking the curve, map, fruit identity, permutations,
and central involution with exact fractions. Among `G,...,12G`, the first
positive triple is the explicitly repaired input at 9G. The complete
positivity census of these 72 points is in the JSON. This finite search
does not prove that the result is least among all positive solutions.

The small point G is of infinite order: it is not in the full six-element
torsion group (8). Standard real elliptic-curve topology therefore makes
its odd multiples dense in the bounded real component containing G. The
positive locus is open and nonempty, since it contains 9G; hence infinitely
many odd multiples give positive rational, and after scaling positive
integer, solutions. This uses topological density, not a numerical sample.

The canonical-height identity `hhat(nG)=n^2 hhat(G)` explains why a modest
multiple can acquire very large coordinates. It does not make digit count
exactly n squared, predict the first positive multiple from n alone, or
prove global leastness. A complete leastness argument additionally needs
a generator/saturation statement, all torsion cosets, and certified height
bounds or the original paper's separate argument. None is supplied by
this bounded audit.

## 7. Comparison boundary and reproduction

The other proposed square curve `y^2=2x^3+4x^2+1` becomes
`Y^2=X^3+4X^2+4` under `X=2x,Y=2y`. Its j-invariant is `-65536/91`, whereas
the fruit curve has `1408317602329/2153060`. They are not isomorphic over
an algebraic closure. At the common good prime 11 their point counts are
16 and 12, so they are not Q-isogenous either. The latter uses the standard
invariance of the Frobenius point count under isogeny, as discussed in
Elkies's Chapter VII notes. These invariants independently support the
parent [polynomial lane](catalan_elliptic_20260921_polynomial.md).

Run from the repository root:

```powershell
python 04-computation/experiments/catalan_elliptic_20260921_elliptic.py
python -O 04-computation/experiments/catalan_elliptic_20260921_elliptic.py
```

The standard-library-only script uses explicit exception checks, with no
`assert` or floating-point gates. It compares every coefficient of the
general polynomial identities, audits the literal and repaired integers,
checks the six torsion points and 72 translated multiples, tests the
Fibonacci family at all odd n from 3 through 31, and counts both curves
over F_p for p in `{3,11,17,19}` by two methods. The companion JSON
contains the source LF hash, exact universe, all check labels, and complete
finite results. Normal and optimized runs must produce identical decoded
JSON. The three-to-six comparison preserves a group action only; no Collatz,
quadratic-cycle classification, or unrelated elliptic-curve equivalence
follows.
