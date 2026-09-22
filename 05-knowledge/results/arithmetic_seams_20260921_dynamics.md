# Arithmetic seams: genuine three-to-six lifts and a twelve-element symmetry

**2026-09-21. Status: PROVED elementary statements, inherited proved results,
FINITE-EXACT controls, and a separately UNRESOLVED external reference.**
No novelty claim. No Collatz convergence, classification of all rational
quadratic cycles, or sphere construction follows.

The useful extension is an explicit `S_2 x S_3` action on a six-vector lift
of the rational three-cycle of `x^2-29/16`. The construction extends to
every rational projective three-cycle. This is an action on a specified
finite set, not an identification with an unspecified Anthropic result.

## Inheritance and the live board

The closest proved mechanism is the trace-one determinant-one lift in
[THM-4146, rational three-cycle order-six lift and fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md),
extending [THM-4139, rational three-cycle and horizontal-carrier separation](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md).
Both the rational cycle and its six-vector hexagon are inherited. The new
comparison here adds the explicit reflection, the direct-product
decomposition, a universal finite-set construction, and the complete
Chebyshev period-six conductor table.

The canonical hostile is already in THM-4139: a projective or linear
order-six construction does not produce a rational scalar six-cycle.
The corrected near miss is the attempted identification of a three-point
divisor with a single elliptic fibre; equal source sets do not identify
target dynamics. The least-used sidecar here is the action on the omitted
sheet: central sign in the vector lift, inversion in the trace quotient,
and time orientation in the three-point set.

Other inherited inputs are the
[quadratic reversal note](arithmetic_braids_20260917_geometry.md),
[Zsigmondy triad](collatz_mod6_20260917_zsigmondy_triad.md), and
[Pythagorean semicircle note](collatz_mod6_20260917_pythagorean_semicircle.md).
The last already distinguishes the conductor-seven and conductor-nine
Chebyshev three-cycles; those are not new discoveries here. The relevant
research moves are **type every analogy**, **controlled forgetting needs a
sidecar**, and **search the statement before the method**.

| Lane | Object and retained predicate | Lost coordinate and decisive test |
|---|---|---|
| Anchor | Rational quadratic cycle; exact scalar period | Field and actual map; enumerate the denominator/escape survivors |
| Niche | Six-vector lift; linear order and norm | Central sign; compare its cube with the identity |
| Wildcard | Reflection plus rotation; actual group action | Forward orientation; test `RBR=B^-1` rather than commutation |
| Bridge | Roots of unity to Chebyshev traces | Inversion sheet; distinguish `2^m=1` and `2^m=-1` modulo conductor |
| Source audit | Claimed `S_2 x S_3` construction | Spheres versus permutation groups and exact bibliographic identity |

## 1. Three different observations, with their starts specified

For `a_0=0`, the recurrence `a_(n+1)=2a_n+1` has `a_n=2^n-1`.
At `n=6`,

```text
63=3^2*7,    3 divides 2^2-1,    7 divides 2^3-1.
```

Thus no prime divisor is new at that index. In contrast, `9` is a new
prime-power divisor: `ord_9(2)=6`, although `ord_3(2)=2`. The no-new-prime
statement loses valuation depth. The inherited Bang/Zsigmondy theorem
classifies the whole sequence; the present script independently checks
indices `1..40`, finding exactly `1,6` as exceptions under the convention
that the term `1` has no prime divisor. No global theorem is inferred from
that finite check.

For `f(x)=x^2-7/4` **starting at zero**, the first three iterates are

```text
-7/4,    21/16,    -7/256.
```

The third numerator has no new prime, but the rational orbit has not
returned: denominators differ. The inherited reversal map
`L(x)=-x-1/2` transports the reversed conductor-seven three-cycle of
`x^2-2` to a parabolic algebraic three-cycle of this map. The companion
conductor-nine cycle goes instead to parameter `-11/4`, with multiplier
`19`, so a bare three-cycle is insufficient to select `-7/4`.
The primitive-divisor question concerns an orbit's numerators; the
parabolic-cycle question concerns a different starting point. See
[Krieger's primary paper](https://arxiv.org/abs/1203.2555) for the former
framework; no uniform conclusion beyond its stated hypotheses is imported.

For `g(x)=x^2-29/16` **starting at -7/4**, the rational orbit is

```text
-7/4 -> 5/4 -> -1/4 -> -7/4.                              (1)
```

Its exact period is three and multiplier is `8*(-7/4)*(5/4)*(-1/4)=35/8`.
It is the only affine rational cycle of this particular map. For
self-containment, the inherited proof is short: a preperiodic rational
point cannot have an odd denominator prime, since its negative valuation
would double at each iterate. At two, every valuation other than `-2`
eventually doubles downward. Hence a preperiodic point is `u/4` with `u`
odd. Its recurrence is `(u^2-29)/4`; every `|u|>=9` escapes, while the
eight remaining odd numerators `+-1,+-3,+-5,+-7` give exactly the graph in
THM-4139. In particular there is **no rational scalar six-cycle** at this
parameter. On the projective line infinity is an additional fixed point.

This is parameter-specific. Stoll's
[Rational 6-cycles under iteration of quadratic polynomials](https://arxiv.org/abs/0803.2836)
states a general exclusion subject to analytic and Birch--Swinnerton-Dyer
assumptions; its Theorem 7 is not an unconditional input here. We make no
all-parameter rational-period-six claim.

## 2. Every rational projective three-cycle has a dihedral twelve-element lift

**PROVED.** Let `p0,p1,p2` be three distinct points of `P^1(Q)`, and let
`M` be the unique projective linear transformation cycling them in that
order. A projective map fixing three distinct points is the identity, so
`M` has exact order three. It has a unique trace-one determinant-one
matrix representative `B`, with

```text
B^2-B+I=0,             B^3=-I,             B^6=I.         (2)
```

To see the normalization, the ratio of eigenvalues of a nontrivial
projective order-three matrix is a primitive cube root of unity. Thus
`trace(A)^2=det(A)` and the nonzero trace normalizes `A` to `B`; equivalently
these identities follow by conjugating the three marked points to
`infinity,0,1` and checking their unique cyclic interpolant.

Choose a nonzero rational vector `v` representing `p0`. The vectors
`v,Bv` are independent since their projective points differ. In this basis,

```text
B_0 = [[0,-1],[1,1]],        J_0 = [[1,1],[0,-1]].        (3)
```

Define `J` by this basis change. Then

```text
Jv=v,    J(Bv)=v-Bv=B^-1 v,
J^2=I,   JBJ=B^-1.                                      (4)
```

The six vectors `H={B^k v:0<=k<6}` are distinct: their projective images
first repeat at step three, when the vector changes sign. Formula (4)
shows that `J(B^k v)=B^(-k)v`, so `J` preserves `H`. The normal forms
`B^k J^e`, `0<=k<6`, `e in {0,1}`, are all distinct: rotations have
determinant `1`, reflections determinant `-1`, and the six rotations
are distinct. Consequently

```text
<B,J> = D_6  (the dihedral group of order 12)
      = <B^3> x <B^2,J>  isomorphic to C_2 x S_3
      = S_2 x S_3.                                      (5)
```

For the direct product, `B^3=-I` is central, the second factor is the
usual order-six group generated by a three-cycle and its reversing
involution, their intersection is trivial, and `B=B^3*(B^2)^2`.
Here the notation `D_6` counts six polygon vertices; its order is twelve.

The projection `H -> {p0,p1,p2}` forgets the central sign. The full
projective image is `S_3`, with kernel `{I,-I}`. Only its cyclic subgroup
`C_3` preserves the **directed** three-cycle. Reflections reverse arrows;
calling all of `S_3` a commuting dynamical symmetry would be false.

For a three-cycle of `x^2+c` with trace `sigma`, the inherited explicit
choice is

```text
B_sigma = [[sigma+1, -(sigma^2+sigma+1)], [1,-sigma]].    (6)
```

Its projective map agrees with the quadratic on the three points.
Equation (5) therefore applies, but it actually needs no quadratic
polynomial at all. Agreement on a finite set is not a global conjugacy:
one map has degree one, the other degree two.

## 3. The `-29/16` cycle supplies a particularly simple reflection

In coordinates `y=4x`, the inherited lift and the added reflection are

```text
B = (1/4)*[[1,-13],[1,3]],
R = [[-1,-2],[0,1]],             y -> -y-2.              (7)
```

The projective reflection is precisely `x -> -x-1/2`. Direct multiplication
gives `R^2=I`, `RBR=B^-1`. The vector orbit is

```text
(-7,1), (-5,-1), (2,-2), (7,-1), (5,1), (-2,2).          (8)
```

It is the complete integral level set of

```text
Q(X,Y)=X^2+2XY+13Y^2=(X+Y)^2+3(2Y)^2=48.              (9)
```

Indeed `|Y|<=2` and `|X+Y|<=6`; the resulting finite enumeration is exact.
Both `B` and `R` preserve (9) and permute (8), giving (5) in these concrete
coordinates. The script also independently enumerates every invertible
rational linear map determined by two hexagon vertices and finds exactly
these twelve setwise symmetries.

Put `U=X+Y`, `V=2Y`, and `w=U+i*sqrt(3)*V`. Then

```text
B: (U,V) -> ((U-3V)/2,(U+V)/2),
R: (U,V) -> (-U,V),
B: w -> ((1+i*sqrt(3))/2)*w,      R: w -> -conjugate(w). (10)
```

Thus the six-vector lift is an exact Eisenstein hexagon, with rotation
by a sixth root of unity and reflection. This is a map preserving a
specific norm and finite set, not merely a matching number of vertices.

There is also a precise exceptional statement. Among exact three-cycles
of centered monic quadratics, invariance of the underlying set under
`x -> -x-1/2` forces its trace to be `-3/4`. The inherited identity
`c=-(sigma^2+sigma+2)` then forces `c=-29/16`; conversely (1) is invariant.
This fixed affine reflection is special. The abstract dihedral lift (5)
is generic and its reflection usually requires a different rational
projective transformation.

The central-sign lift does not supply any new scalar return: (8) projects
to `-7,5,-1,-7,5,-1`. Nor is it an iteration of Gaussian squaring: (10)
multiplies by a fixed sixth root, while squaring multiplies the angle by
two. Constant rotation and angle doubling have different dynamics.

## 4. A second, distinct three-to-six mechanism: the Chebyshev trace quotient

**PROVED.** Let `pi(z)=z+z^-1` for `z!=0`. Then

```text
pi(z^2)=pi(z)^2-2.                                     (11)
```

The omitted sheet is inversion `z <-> z^-1`. Since
`pi(u)=pi(v)` iff `u=v` or `u=v^-1`, a periodic trace has a root-of-unity
lift. If `z` has odd order `N>1`, its trace period is exactly

```text
m = min{r>=1: 2^r=1 or -1 (mod N)}.                   (12)
```

When the first return has sign `+`, squaring has period `m`; a trace cycle
lifts to two separate length-`m` cycles exchanged by inversion. When the
first return has sign `-`, squaring has period `2m`; one length-`2m`
cycle folds to the trace cycle. These alternatives must not be merged.

Every nonconstant three-cycle of `x^2-2` comes from conductor `7` or `9`:
`N` must divide `2^3-1=7` or `2^3+1=9`, with the lower-period conductors
removed. The two cubics are

```text
Psi_7(x)=x^3+x^2-2x-1,
Psi_9(x)=x^3-3x+1.                                    (13)
```

At conductor seven, squaring already has period three and the two
three-cycles are paired by inversion. At conductor nine, squaring has
period six and folds to one three-cycle. The common shape of a twofold
cover over three points does not identify (11) with the vector lift (7):
the traces are cubic irrational numbers, while (1) is rational; the maps
and fields differ.

### Complete period-six conductor table for `x^2-2`

**PROVED.** Exact period six requires `N` to divide `63` or `65`, with all
divisors producing periods `1,2,3` removed. The survivors are exactly

| `N` | Squaring period | Trace period | Trace points `phi(N)/2` | Trace cycles | Multiplier |
|---|---:|---:|---:|---:|---:|
| 13 | 12 | 6 | 6 | 1 | -64 |
| 21 | 6 | 6 | 6 | 1 | 64 |
| 63 | 6 | 6 | 18 | 3 | 64 |
| 65 | 12 | 6 | 24 | 4 | -64 |

For instance the proper-divisor cases for `63` are `1,3,7,9,21`, and
only `21` survives; for `65` they are `1,5,13`, and only `13` survives.
Within any fixed conductor all primitive roots have the same period.
Distinct primitive roots give equal traces exactly for the inversion
pair, giving the displayed point counts. All traces are real. Thus
`x^2-2` has exactly nine real cycles of period six, comprising 54 distinct
points. None is rational: an irreducible cyclotomic trace has degree
`phi(N)/2>1` for all four conductors. Alternatively the elementary
rational denominator/escape proof gives its complete rational
preperiodic set `{-2,-1,0,1,2}`.

The multiplier follows directly by differentiating (11): at a periodic
trace with `z!=+-1`,

```text
(f^m)'(pi(z)) = 2^m*(z^(2^m)-z^(-2^m))/(z-z^-1),       (14)
```

which is `+2^m` or `-2^m` according to the sign in (12). In particular
the two Chebyshev three-cycle multipliers are `8,-8`, distinct from
`35/8` in (1). A local differentiable dynamical conjugacy preserving
these cycles would have to preserve their multipliers; finite-set
interpolation does not do so.

For an independent exact check, define `Psi_N` by
`Phi_N(z)=z^(phi(N)/2)*Psi_N(z+z^-1)`. Integer polynomial arithmetic gives

```text
((f^6(x)-x)*(f(x)-x))/((f^3(x)-x)*(f^2(x)-x))
  = Psi_13(x)*Psi_21(x)*Psi_63(x)*Psi_65(x),
f(x)=x^2-2.                                             (15)
```

This degree-54 factorization is an identity, not a numerical root census.
The enumeration of conductors and (12) prove exact periods independently
of possible dynatomic specialization conventions.

The `63` seam is now precise: its prime divisors have orders two and
three, but its prime-power and composite modulus can support order six.
Conductor nine folds that order down to three; conductor63 does not.
Hence absence of a new prime divisor cannot determine the observed
period until valuation depth, conductor and quotient have been supplied.

## 5. The Anthropic reference remains a separate identification question

**UNRESOLVED as of this read.** Searches for the exact `S_2 x S_3` wording
did not identify a primary recent Anthropic construction. This is a search
result, not a proof that no such source exists. The following objects are
different and are not substituted for one another:

* `S_2 x S_3` as permutation groups: the order-twelve group explicitly
  constructed in (5).
* `S^2 x S^3` as a five-dimensional manifold: the primary
  [Gauntlett--Martelli--Sparks--Waldram paper](https://arxiv.org/abs/hep-th/0403002)
  constructs Sasaki--Einstein metrics and is dated2004. It does not identify
  the requested recent Anthropic source.
* `S^6` as a six-dimensional sphere: the
  [108-page manuscript at Alpoge's site](https://alpo.ge/s6.pdf) uses a
  `(3,4,infinity)` torus construction. A newer primary exposition,
  [Philip Engel, Complex structures on S6](https://philip-engel.github.io/S6.pdf),
  is dated September13,2026 and has25 pages when checked here. Its abstract
  explicitly attributes the underlying manuscript to Claude under
  Alpoge's direction. Its opening describes an elliptic-surface and torus
  construction, followed by logarithmic modifications of multiplicities
  three and four. It is not an `S^2 x S^3` construction.

The last source updates attribution and bibliographic identity relative
to the repo's [August Hopf/S6 ledger](../reference/CORE-PAPERS-HOPF-S6-2026-08-24.md).
This lane has not independently audited its global analytic/topological
argument and imports no sphere theorem from it. Likewise the separate
[Brendle--Hung S2 x S2 ledger](../reference/CORE-PAPERS-BRENDLE-HUNG-S2XS2-2026-08-24.md)
concerns a different manifold and curvature question. A link or exact
title from the user is needed to finish identifying the intended reference.

## 6. Reproduction and scope

Run:

```text
python 04-computation/experiments/arithmetic_seams_20260921_dynamics.py
python -O 04-computation/experiments/arithmetic_seams_20260921_dynamics.py
```

The [script](../../04-computation/experiments/arithmetic_seams_20260921_dynamics.py)
uses only exact integers, fractions and polynomial arithmetic. Its
[JSON](../../04-computation/experiments/arithmetic_seams_20260921_dynamics.json)
records125 distinct rational chart parameters `t=a/d` with `|a|<=12`,
`1<=d<=8`, excluding `0,-1`; for each it checks the scalar three-cycle,
the six-vector lift and the twelve-element group. It separately checks
the complete Q=48 integer hexagon, all odd conductors `3..129`, the
degree54 factorization, the critical initial terms, and Mersenne indices
`1..40`. Checks remain active under optimized Python. The finite controls
corroborate the general proofs and do not replace their quantifiers.

The resulting connection ledger is exact: vector projectivization loses
central sign; the Chebyshev trace loses inversion; the full permutation
group loses arrow orientation; prime support loses valuation depth. Each
loss has an explicit restoring coordinate and a hostile example above.
