# Four-dimensional Kakeya: focal selectors, quadratic carriers, and atlas seams

> **Research synthesis, 2026-08-26.**  Euclidean Kakeya in `R^4` remains
> **OPEN**.  The proved outputs of this session are
> [THM-4235](../01-canon/theorems/THM-4235-finite-four-dimensional-torus-kakeya-quadratic-carrier.md),
> which is **FINITE-EXACT + VERIFIED**, and
> [THM-4236](../01-canon/theorems/THM-4236-four-dimensional-kakeya-monic-cubic-focal-selector.md),
> which is **PROVED** over the reals.  Neither gives a new general Hausdorff-
> dimension or maximal-function exponent.  This document records the
> mechanism, hostile controls, connection contracts, and the next open proof
> obligations; it is not a replacement for the canon.

## Outcome first

The session changes the preferred state variable for a selected line.  A
direction and an arbitrary basepoint are redundant; the invariant coordinate
is the line's intersection with a fixed transverse hyperplane.  On a bounded
slope chart write

```text
line(u)={(t,a(u)+tu):t in I_u},
u in U subset R^3.                                  (1)
```

The incidence map `F(u,t)=(t,a(u)+tu)` has focal determinant

```text
P_F(u,t)=det(Da(u)+tI_3),                            (2)
```

a monic cubic in the line parameter.  This yields three rigorous squeezes on
a hypothetical null Kakeya selector:

1. a Lipschitz intercept section on any positive-measure slope patch is
   impossible, because it sweeps positive four-volume;
2. an approximate polynomial atlas pays the explicit cubic complexity
   `C_delta=sum_j max(d_j,1)^3`; and
3. a two-ends shading cannot hide all its mass near the cubic's three focal
   times, changing the natural shading loss from quartic to linear.

The exact finite-field carrier supplies the hostile and positive controls for
this picture.  Over `F_61`, the same `216,000` torus directions admit affine
unions differing by a factor `7.1426...`; hence direction broadness alone is
not an affine-union invariant.  The fully quadratic placement has maximum
multiplicity eight and rank-four cubes at every maximally incident point.
Four boundary charts complete it to one line in every direction of
`P^3(F_61)`, with union size `1,868,641` and maximum multiplicity twelve.

The resulting frontier is no longer “attach some basepoints.”  It is an
atlas/seam/algebraic-parent trichotomy: prove that a small Euclidean tube union
forces either a low-complexity intercept atlas, transverse collisions between
atlas pieces, or a slowly moving narrow flag that the planebrush can consume.

## Inheritance pass

- **Closest proved mechanism.**  THM-4035 separates a full-spark twisted
  cubic from a fixed-plane hostile while preserving the same `C_60` clock.
  Analytically, endpoint four-linear Kakeya controls the broad regime and the
  corrected Katz--Zahl planebrush controls part of the narrow regime.
- **Canonical hostile at entry.**  The same finite direction clock can be
  maximally transverse or wholly planar.  Thus phase, period, direction count,
  and even many determinant statistics do not determine a Kakeya union.
- **Corrected near miss.**  THM-4035's three-clock torus supplied directions,
  not a Kakeya set: basepoints, boundary charts, multiplicity, shading, and
  scale were missing.  THM-4235 closes only the first three finite debts.
- **Least-used relevant sidecar.**  The affine-intercept section `a(u)`,
  conditioned simultaneously on direction, parent tube, scale, and shading,
  had not been made the native object.
- **External stopping boundary.**  Finite-field cardinality, multilinear
  transversality, sticky bounds, and multiplicity factoring have different
  hypotheses from the general Euclidean conjecture.  None may be imported
  without its missing sidecars.

## Live concept board

| live object | representation | invariant / observable | operation or scale | information destroyed by the cheap quotient |
|---|---|---|---|---|
| affine line section | graph `(u,a(u))` in a six-dimensional line chart | regularity and approximate polynomial degree | restrict slope patch; change transverse chart | direction-only projection loses placement and seams |
| incidence projection | `F(u,t)=(t,a(u)+tu)` | focal determinant `P_F=det(Da+tI_3)` | project line space to physical space | image union loses fibre multiplicity and incident direction flags |
| shaded segment | measurable `Y_u subset I_u` | density, all-translate two-ends, distance from focal roots | cut and rescale along a line | total shaded length loses concentration near critical times |
| polynomial atlas | pieces `(U_j,a_j,d_j)` | `C_delta=sum max(d_j,1)^3` | refine patches; cross a seam | scalar complexity loses adjacency and moving narrow planes |
| finite quadratic carrier | `s -> (s_1^2,s_2^2,s_3^2)` over `F_61` | union, multiplicity, pinned ranks, chart label | add boundary charts; shuffle placements | reduction mod 61 loses angle, thickness, two-ends, and nested parents |
| tube parent forest | child/parent indices and shadings | multiplicity product, relative Frostman cost, algebraic type | pass from `delta` to `rho` | coarse multiplicity loses exceptional parents and ruled ancestry |

The board deliberately has no tournament.  There is no intrinsic binary
orientation whose winner preserves the Kakeya target.  Fibres, critical
times, parent flags, and seams are the natural vertices here.

## Anchor / Niche / Wildcard portfolio

### Anchor -- Euclidean atlas/seam dichotomy

The main target is a scale-dependent structural theorem for the intercept
section.  Given a `delta`-tube family on a bounded slope patch, construct a
piecewise polynomial approximation `a_delta` and prove that at least one of
the following has quantitative mass:

```text
A. low atlas debt:       C_delta <= delta^(-beta);
B. transverse seams:     many cross-patch collisions have a four-volume gain;
C. coherent narrow flag: one slowly varying three-plane field captures the
                         remaining directions across parent scales.
```

In a scale-uniform neighborhood argument, branch A alone would force lower
Minkowski dimension greater than `s` only if `beta<4-s`.  Numerically, using
the corrected current Hausdorff record as a sensitivity benchmark gives

```text
4-3.059849573... = 0.940150426....                   (3)
```

This numerical comparison is not itself a Hausdorff consequence; that would
also require a multiscale induction or content argument of the correct type.
The trichotomy is intentionally not just an approximation theorem.  Without
connected patch geometry and seam adjacency, one may partition into arbitrary
measurable dust and make `C_delta` a content-free label.

### Niche -- exact finite quadratic placement

THM-4235 supplies a reproducible laboratory in which directions, placement,
multiplicity, and local rank can be changed one coordinate at a time.  The
next finite experiment is not another union count.  It is a constrained
shuffle:

```text
preserve the full-projective marginal multiplicity histogram;
preserve the pinned rank histogram;
shuffle basepoint ownership;
measure higher incidence and two-scale polynomial-parent concentration.
```

If the higher statistics also remain unchanged, the finite carrier has
reached its stopping boundary.  If they move, the first moving statistic is a
candidate sidecar for a Euclidean parent tree, not a Euclidean theorem.

### Wildcard -- caustics as the Jacobian-side analogue

The planar Jacobian work suggests looking at the critical divisor of the
incidence projection rather than at directions alone.  Here this is a typed
analogy, not a reduction to the Jacobian conjecture:

```text
source:       selected section of affine line space;
map:          F(u,t)=(t,a(u)+tu);
critical law: det(Da+tI_3)=0;
target:       physical Kakeya union.
```

The monic cubic permits only three focal times per differentiable direction.
A null image must therefore exploit roughness of `a`, concentration of
shadings near those times, or very large global multiplicity.  Constant-
Jacobian, properness, and polynomial-invertibility statements from the
Jacobian conjecture do not transfer: `a` may be merely measurable, the domain
is restricted, and the determinant is not constant.

## Pull I -- the sharp focal-selector theorem

For a selected unit segment in the slope chart, `|I_u|` is bounded below by
`ell_R=(1+R^2)^(-1/2)`.  At almost every differentiability point of a
Lipschitz `a`, the area formula gives

```text
integral N_F(y)dy
  = integral_U integral_(I_u) |P_F(u,t)|dtdu.         (4)
```

Every monic real cubic on an interval of length `L` obeys the sharp law

```text
integral_I |p(t)|dt >= L^4/64.                       (5)
```

The extremizer is the translated and scaled polynomial
`2^(-3)U_3=x^3-x/2`; its sign annihilates every quadratic.  Therefore

```text
integral N_F >= (1/64) integral_U |I_u|^4 du > 0.    (6)
```

If the image were null, the left side would vanish.  No finite multiplicity
assumption is needed for this qualitative conclusion.  The concurrent-star
hostile `a(u)=p'-cu` has `P_F=(t-c)^3`; it is one-to-one almost everywhere
but has an infinite critical fibre at the apex, so every quantitative
multiplicity statement must retain the almost-everywhere qualifier.

The scalar mechanism extends exactly to `R^n`.  For `m=n-1`, the sharp
degree-`m` law is

```text
integral_I |p(t)|dt >= L^(m+1)/2^(2m),               (7)
```

with extremizer `(L/4)^m U_m(2(t-c)/L)`.  This explains the cubic and the
constant `1/64` as the dimension-four member of a stable incidence law.

### What this actually squeezes

A null Kakeya selector cannot contain a positive-measure slope patch on which
its invariant intercept is a restriction of a Lipschitz map.  “Endpoint is
rough” is not the right invariant because shifting an endpoint along its line
changes it; `a=b'-b_0u` is invariant.

The maximal-gradient decomposition sharpens this to a **PROVED Sobolev
exclusion**.  On no nonempty open slope patch can the selected affine
intercept agree almost everywhere with a `W^(1,p)_loc` map for any `p>=1`.
After cutoff, the precise `W^(1,1)` representative is Lipschitz on the level
sets of its maximal gradient, and one such set has positive measure.  The
selector theorem then contradicts null volume.  At `p=1` this uses only that
the maximal function is finite almost everywhere, not a strong `L1` bound or
a Sobolev area formula.  The next functional-analytic frontier is `BV`,
fractional Sobolev, or quantitative Lusin-Lipschitz complexity.  No extension
to those classes is currently proved.

## Pull II -- shadings can hide only at three focal times

For a monic cubic `p=product_(j=1)^3(t-z_j)`,

```text
|{|p|<epsilon}| <= 6 epsilon^(1/3).                  (8)
```

Consequently, an arbitrary measurable shading `Y_u` of length `lambda(u)`
satisfies

```text
integral_(Y_u)|P_F|dt >= lambda(u)^4/864.            (9)
```

The quartic order is real: shade a concurrent star in an interval of length
`lambda` adjacent to its triple focal time.  Thus a density-only argument
cannot replace `lambda^4` by `lambda`.

All-translate two-ends is exactly the missing distribution predicate.  If

```text
|Y_u intersect J|
 <= A (|J|/L(u))^alpha |Y_u|                        (10)
```

for every translated interval `J`, remove three intervals around the focal
roots at radius

```text
r=(L(u)/2)(6A)^(-1/alpha).
```

At least half the shading remains and `|P_F|>=r^3` there, giving

```text
integral_(Y_u)|P_F|dt
 >= [L(u)^3/16](6A)^(-3/alpha)|Y_u|.                (11)
```

This is the clean analytic interface: two-ends converts quartic shading loss
to linear shading loss, but it does not bound how many direction/atlas fibres
project to the same physical point.

## Pull III -- polynomial multiplicity becomes an atlas debt

If `a` is polynomial of coordinate degree at most `d`, then a preimage of
`y=(y_0,y')` solves

```text
a(u)+y_0u-y'=0.                                     (12)
```

At a regular preimage these are three isolated simple equations of degree at
most `D=max(d,1)`.  Affine isolated-zero Bezout gives

```text
N_F(y)<=D^3 for almost every y.                      (13)
```

For a finite atlas, sum `D_j^3`.  If a piecewise polynomial model approximates
the invariant intercept to error `delta`, its sweep lies in `N_delta(E)` and

```text
|N_delta(E)|
 >= [integral_U |I_u|^4du]/[64 C_delta],
C_delta=sum_j max(d_j,1)^3.                          (14)
```

Under two-ends, the numerator may instead use the linear-density expression
from `(11)`.  If a hypothetical set of dimension `s` has
`|N_delta(E)| <= C delta^(4-s)` along a scale sequence, then every such atlas
pays

```text
C_delta >= c delta^(-(4-s)).                         (15)
```

This turns “the selector must be complicated” into an exponent-bearing debt.
It does not produce a useful atlas.  The missing inverse theorem must obtain
one from small union volume or route its failure to transverse seams or a
coherent narrow flag.

## Pull IV -- exact finite placement and boundary completion

Let `q=61` and use the nonzero torus directions

```text
D*={[1:u:v:w]:u,v,w in F_q^*}.
```

Changing the first `k` transverse intercept coordinates from zero to squares
gives the exact placement ladder

```text
k=0: 12,960,001
k=1:  6,696,030
k=2:  3,460,500
k=3:  1,814,460.                                    (16)
```

Every row has identical directions and global plane counts.  The fully
quadratic row has multiplicity histogram

```text
{1:480, 2:20,880, 4:302,760, 8:1,490,340}.          (17)
```

At every multiplicity-eight point, the incident directions are affinely the
cube `(1,+/-1,+/-1,+/-1)`.  Exactly `58` of its `70` quartets are transverse.
Thus abundant pinned broadness coexists with a union more than seven times
smaller than the concurrent placement using the same direction set.

Normalize all projective directions by their first nonzero coordinate.  The
four disjoint charts have sizes `q^3,q^2,q,1`; attach the recursive quadratic
lines from THM-4235.  Two independent exhaustive paths give

```text
|K_quad|=1,868,641,
max multiplicity=12,
total incidence=14,076,604=(q^3+q^2+q+1)q.          (18)
```

The density is `0.13496045...`, only a factor `1.05335...` above the
Bukh--Chao lower expression at `q=61`.  For every odd prime power `q`, the
same four-chart hierarchy has the elementary upper bound

```text
q sum_(r=0)^3 ((q+1)/2)^r=q^4/8+O(q^3).             (19)
```

The finite multiplicity-eight mechanism is the reduction of

```text
F(s,t)=(t,s_1^2+ts_1,s_2^2+ts_2,s_3^2+ts_3),
P_F=product_i(2s_i+t).                               (20)
```

The degree-two Bezout ceiling `2^3=8` is attained on the torus chart.  This is
an exact algebraic shadow of the Euclidean focal theorem.  It preserves
degree and fibre branching; it destroys the metric data needed by Kakeya.

## How each pull changes the board

| new fact | affine section | focal divisor | shading | atlas/seam | finite carrier | parent forest |
|---|---|---|---|---|---|---|
| sharp cubic L1 | forbids Lipschitz positive patches | makes three focal times exact | identifies where density may hide | supplies a volume reward for low complexity | predicts the threefold quadratic fold | asks for a differentiable representative inside parents |
| two-ends bound | leaves intercept roughness untouched | removes neighborhoods of roots | changes quartic to linear loss | improves the reward branch | has no finite analogue without scale | can feed a shading-compatible parent refinement |
| Bezout atlas | makes polynomial complexity measurable | controls regular fibre count | combines with either density law | exposes seams as the missing coordinate | explains maximum torus multiplicity eight | suggests algebraic-parent labels rather than raw multiplicity |
| exact quadratic carrier | proves placement is independent of broadness | realizes all three folds | supplies no two-ends information | offers constrained shuffle hostiles | closes directions/basepoints/multiplicity | demands a two-scale enrichment before transfer |
| ruled-quadric hostile | shows smooth algebraic structure can still be dangerous | criticality alone is not the whole obstruction | survives sparse selection | requires parent type in every patch | has no direct finite transfer | blocks convex-Wolff-only induction |

## Source-target contracts

### Finite carrier to Euclidean tubes

```text
source:     labelled lines in F_61^4
target:     delta-tubes in R^4
map:        coefficient lift of the quadratic incidence polynomial
preserves:  direction chart, polynomial degree, algebraic fibre pattern
destroys:   angular size, tube thickness, metric volume, two-ends, scale
sidecar:    nested parents, shadings, Wolff count, moving narrow flag
cheap test: compare quadratic, concurrent, and shuffled placements at two scales
```

No theorem currently supplies this transfer.

### Approximate atlas to physical union

```text
source:     pieces (U_j,a_j,d_j) and parameter shadings
target:     N_delta(E)
map:        (u,t) -> (t,a_j(u)+tu)
preserves:  direction, transverse placement up to delta, focal degree
destroys:   patch adjacency after taking the union
sidecar:    seam graph, collision angle, narrow-plane flag, essential multiplicity
cheap test: concurrent star, focal shading, and two crossing polynomial patches
```

### Buffered parent refinement to planebrush recursion

```text
source:     indexed convex parent-child families with compatible shadings
target:     one induction-on-scale recurrence
map:        refine to constant inner/global/coarse multiplicities
preserves:  multiplicity product and buffered coarse shading
destroys:   exceptional-parent mass unless separately weighted
sidecar:    relative Frostman distribution and algebraic-parent type
cheap test: one delta-child inside one rho-parent, then the ruled quadric family
```

The single-child test has relative Frostman cost comparable to
`(rho/delta)^3`, so uniform admissibility is false.

## Hostile-control battery

Every proposed lemma should be run against all five controls before it is
promoted:

1. **Concurrent star:** excellent direction spread, large union, one critical
   apex; catches confusion between global and almost-everywhere multiplicity.
2. **Focal shading:** `P_F=(t-c)^3` shaded next to `c`; proves quartic density
   loss is sharp without two-ends.
3. **Fixed-direction placement ladder:** same torus determinants, union ratio
   greater than seven; catches direction-only certificates.
4. **Ruled quadric `ad-bc=1`:** a three-dimensional line family whose thinned
   and replicated tube model is a near miss under convex Wolff axioms; catches
   omission of algebraic ancestry.
5. **One child in one parent:** catches an unjustified uniform relative
   Frostman import in multiplicity factoring.

## Ranked next proof moves

### 1. Prove a quantitative atlas/seam/flag trichotomy

This is the sharpest new route.  Use dyadic slope cubes and approximate
intercept oscillation to build a patch graph.  A useful theorem must output
one of:

- `C_delta<=delta^(-beta)` with `beta<0.940150426...`;
- a lower bound for transverse cross-patch incidence that survives projection;
- or a parent-consistent three-plane flag with controlled variation.

The seam graph must remember which physical tubes collide; boundary measure
of slope patches alone is not enough.

### 2. Replace uniform Frostman by a taxed exceptional-parent theorem

Dyadically label each parent by its relative Frostman cost.  Seek a statement
of the form: good parents admit buffered multiplicity factoring, while either
the total shaded mass of high-cost parents decays with the cost threshold or
those parents concentrate near a named low-degree ruled variety.  The second
exit must feed the algebraic-parent branch, not be discarded as an error term.

### 3. Add a focal-distance sidecar to two-ends induction

For each slope/parent, store the three roots of `P_F` when an approximate
derivative exists, or a certified roughness label when it does not.  Measure
the shaded mass outside radius-`r` focal neighborhoods.  This separates the
one-dimensional two-ends gain from the three-dimensional projection
multiplicity instead of blending both into one scalar multiplicity.

### 4. Build a ruled-parent detector

For each parent tube family, test whether its Plucker or direction-position
coordinates approximately satisfy a bounded-degree polynomial relation.  The
output must distinguish a harmless polynomial graph from a doubly ruled
carrier such as `ad-bc=1`.  A determinant histogram is insufficient; retain
the actual polynomial and its ruling dimension.

### 5. Run the constrained finite shuffle

Preserve `(17)` and the pinned cube-rank data while randomizing or
adversarially optimizing basepoint ownership.  Add a coarse/fine chart by
lifting from a smaller subfield or by nested residue boxes.  The decisive
observable is the first higher-order incidence statistic that predicts union
change after the preserved marginals are fixed.

### 6. Audit exponent sensitivity only after a typed gain exists

The 2025 restriction/maximal preprint improves the claimed maximal exponent
through a recurrence combining restriction, two-ends, and X-ray input.  A new
local inequality matters only if its gain survives the scale recursion in the
same norm and hypothesis class.  Insert the focal/two-ends gain symbolically,
differentiate the fixed-point exponent with respect to that gain, and identify
the minimum uniform saving required before attempting a long proof.  Do not
equate a better local density power with a better Hausdorff exponent.

## Literature status used by this session

- **CITED / published corrected baseline:** Katz--Zahl gives Hausdorff
  `3.059849573...` and a distinct maximal estimate `3.049570923...`; the 2025
  v3 corrects Lemma 7.1.
- **PREPRINT v1 CLAIM:** Borges--Chan--Chen--Liu--Xi--Zhan claims maximal
  `d_0=(159+sqrt(145))/56=3.054314...`; this does not exceed the Hausdorff
  record.
- **PREPRINT v1 CLAIM:** Hua--Yao--Yang gives buffered, shading-compatible
  multiplicity factoring under explicit convex and relative-Frostman
  hypotheses; it claims no new Kakeya exponent.
- **CITED hostile/survey:** Zahl's current survey supplies the ruled-quadric
  near miss and the present structural overview.
- **CITED finite-field benchmark:** Bukh--Chao gives the asymptotically sharp
  finite-field density lower bound.  It is not a Euclidean transfer theorem.

See the maintained
[primary-source audit](../05-knowledge/reference/CORE-PAPERS-KAKEYA-4D-2026-08-24.md)
for links and exact scope.

## Non-implications and stopping boundary

This session does not prove dimension four, improve any general Euclidean
dimension/maximal exponent, transfer the `F_61` construction to Euclidean
tubes, obtain a low-complexity atlas for an arbitrary selector, or verify the
relative Frostman input for every parent.  The new theorem says where a smooth
or low-degree selector must spend volume; the open problem is to force enough
of an arbitrary Kakeya selector into that regime, or make its failure usable
by broad/narrow and algebraic-parent analysis.

## Reproduction

```bash
python3 -B 04-computation/kakeya4d_quadratic_carrier_thm4235.py
python3 -B -O 04-computation/kakeya4d_quadratic_carrier_thm4235.py
sha256sum 04-computation/kakeya4d_quadratic_carrier_thm4235.py
sha256sum 05-knowledge/results/kakeya4d_quadratic_carrier_thm4235.out
```

The optimized and ordinary runs must be byte-identical to
`05-knowledge/results/kakeya4d_quadratic_carrier_thm4235.out`.
