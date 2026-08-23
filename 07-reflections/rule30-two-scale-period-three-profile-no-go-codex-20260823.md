# Exact two-scale period-three Rule 30 profile no-go

**Status:** **PROVED algebraic consequence + VERIFIED-EXACT controls**, now
subsumed by canonical THM-3778.  The full Rule 30 frontier and every prize
problem remain **OPEN**.

The inherited inputs are THM-3511
(`01-canon/theorems/THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile.md`),
THM-3512
(`01-canon/theorems/THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary.md`),
THM-3516
(`01-canon/theorems/THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge.md`),
and the 2026-08-21 adaptive-observers reflection
(`07-reflections/rule30-adaptive-observers-rank30-and-padic-selector-gaps-codex-20260821.md`).

## Statement

For an odd period `n`, let the saturated projective profile map be

```text
R(g)_j=-g_(2j)g_(2j+1)(1-g_(2j+2))/(1-g_(2j)),       (1)
```

with indices modulo `n`.  The declared saturated locus is exactly

```text
product_j g_j(1-g_j)R(g)_j(1-R(g)_j) != 0.           (1a)
```

Thus both applications of `R` and both consecutive-ratio amplitude lifts
used below are defined.  No claim is made about boundary points outside this
open locus.

Period one is the inherited scalar hostile `R(g)=-g^2`.  Period three is the
smallest odd period carrying genuine phase variation.  Over an algebraically
closed characteristic-zero field its exact two-scale locus is as follows.

> **Scratch theorem (period-three two-scale classification).**  Every
> saturated solution of `R^2(g)=g` belongs to one of two ordinary-holonomy
> components:
>
> 1. the fixed point `g=(-1,-1,-1)`;
> 2. the rational curve
>
>    ```text
>    g(t)=(2/(t+1), (1-t)/2, (t+1)/(t-1)),             (2)
>    R(g(t))=g(-t),
>    ```
>
>    with `t in P^1 \ {infinity,-1,1}`.
>
> On the curve, `t=0` is fixed and
> `t in P^1 \ {infinity,-1,0,1}` gives a genuine two-cycle.  The locus is
> rational, not elliptic.  The two nontrivial cubic-holonomy sectors contain
> no saturated point.

No algebraic genuine cycle in (2) is physically admissible over `Q_2`:
odd-unit amplitudes would force simultaneously

```text
nu_2(t+1)=nu_2(1-t)=1,                                (3)
```

which is impossible.  The remaining constant point repeats gap one at both
scales and is excluded by THM-3512's no-consecutive-gap-one law.  Therefore
there is **no physical exact two-scale period-three projective profile**.
This is a no-go for one sharply defined scheme, not a Rule 30 prize result.

This is the first genuine relaxation of the inherited period-five/seven
fixed-profile lane.  Fixedness `R(g)=g` forces `p=-1`, ordinary holonomy, and
one eigenproblem `T_1A=lambda A`.  Two-scale cycling instead permits
`h^3=1`, changes the mark by `h -> h^2`, and replaces `T_1` by the twisted
composition `T_(h^2)T_h`.  The surviving ordinary sector is a real algebraic
extension—a rational involutive family—but it still contains no elliptic
component.  Its holonomy mark also does not replace THM-3516's owner/carry
data; that information remains part of the unrealized Mealy gate.

## 1. Product and holonomy

Write `p=product_j g_j`.  Because multiplication by two permutes odd-period
indices, multiplying (1) gives

```text
product_j R(g)_j=-p^2.                                 (4)
```

Thus `R^2(g)=g` implies

```text
p^3=-1.                                                (5)
```

Choose a nonzero amplitude `A_0` and define

```text
g_j=-A_(j+1)/A_j.                                      (6)
```

At period three,

```text
A_(j+3)=h A_j,                   h=-p,                 (7)
```

and (5) says `h^3=1`.  This is the marked/twisted holonomy that a cyclic
ratio profile loses.  It is not optional: ordinary amplitudes have `h=1`,
while a nontrivial two-scale twist would have primitive cubic holonomy.

## 2. Twisted amplitude linearization

On the three-dimensional space `V_h` of sequences satisfying (7), define

```text
(T_h A)_j=A_(2j)+A_(2j+1).                            (8)
```

In the coordinates `(A_0,A_1,A_2)`,

```text
T_h = [1  1  0]
      [h  0  1]
      [0  h  h].                                      (9)
```

Put `B=T_hA`.  Direct substitution of (6) into (1) gives

```text
R(g)_j=-B_(j+1)/B_j,             B_(j+3)=h^2B_j.      (10)
```

Hence the second scale is controlled by

```text
M_h=T_(h^2)T_h
   =[h+1   1       1]
    [h^2   h^2+h   h]
    [h^3   h^3     h^3+h^2].                          (11)
```

When both profiles are saturated, `A`, `B`, and `M_hA` have no zero
coordinate.  Equality of the final and initial consecutive ratios is then
equivalent to

```text
M_h A=lambda A.                                       (12)
```

Conversely, every eigenline satisfying the same nonvanishing conditions
gives a saturated two-scale profile.  Thus no nonlinear elimination has been
discarded: (12) is an equivalence on the declared saturated locus.

The exact characteristic polynomial is

```text
det(lambda I-M_h)
 =(lambda-h)(lambda-h^2)
  (lambda-(1+h+h^2+h^3)).                             (13)
```

## 3. Complete holonomy classification

### Ordinary holonomy `h=1`

Equation (13) becomes

```text
(lambda-4)(lambda-1)^2.                               (14)
```

The `lambda=4` line is `A=(1,1,1)` and gives the constant fixed profile
`(-1,-1,-1)`.

The `lambda=1` plane has the exact parametrization

```text
A=(alpha+beta,-2beta,beta-alpha).                     (15)
```

Its first image is

```text
T_1A=(alpha-beta,2beta,-alpha-beta).                  (16)
```

All amplitude and sibling-sum coordinates are nonzero exactly when

```text
beta!=0,                    alpha!=+-beta.             (17)
```

This homogeneous boundary is authoritative.  Multiplying the rational
saturation factors after substituting the `t`-parametrization cancels zeros
against poles and simplifies spuriously; it is not a boundary test.

Putting `t=alpha/beta` yields (2), and (16) changes `t` to `-t`.  Therefore
`R(g(t))=g(-t)` and `R^2(g(t))=g(t)`.  The first image is projectively equal
to the source only at `t=0` within the saturated chart.  This proves the
ordinary-holonomy part of the theorem, including every boundary.

There is also an independent direct elimination.  Put
`g=(a,b,-1/(ab))`.  After clearing only the declared saturated denominators,
the three equations `R^2(g)_j-g_j=0` factor respectively as

```text
(a+1)(ab-a+1)=0,
(b+1)(ab-a+1)=0,
(ab-1)(ab-a+1)=0.                                    (18a)
```

The common factor gives the rational curve, while off it the residuals force
`a=b=-1`.  Thus the amplitude eigenline argument has neither lost nor added
an ordinary-holonomy component.

### Primitive cubic holonomy

If `h^2+h+1=0`, equation (13) has the three distinct eigenvalues

```text
h, h^2, 1.                                             (18)
```

Corresponding eigenvectors are

```text
(-1,1,0),              (0,-1,1),              (1,0,-h). (19)
```

Every eigenline has a zero amplitude coordinate.  None defines a saturated
ratio profile.  Thus the marked cubic twists add no hidden component.  There
is also no primitive cubic holonomy in `Q_2`: an integral root would reduce to
a root of `x^2+x+1` in `F_2`, but that polynomial has none.

## 4. Physical valuation and Mealy gates

For an actual Rule 30 scale, use the physical amplitudes

```text
A_j=U_m(t+jq_m).                                       (20)
```

Every `A_j` is an odd 2-adic unit.  THM-3512's sibling trace further says
that every coordinate of `T_1A` has the same positive valuation `d_m`, and
after division by `2^d_m` gives the next odd-unit amplitude profile.

On (15), a common scalar changes all valuations equally.  Equal unit
valuations would therefore require (3), because the middle coordinate is
`-2beta`.  But

```text
(t+1)/2 + (1-t)/2 = 1.                                (21)
```

If both valuations in (3) were one, the two summands on the left would be
odd 2-adic units.  Their sum would be even, contradicting (21).  So the whole
rational two-cycle curve fails before section realization.

The constant line does have odd amplitudes.  Its sibling sums have valuation
one at both scales, giving the forbidden gap word `(1,1)`.  It too is not
physical.

An algebraic profile that survived these tests would still owe the actual
phase owners `s_(m,t)`, the common gap, and exact Mealy-section identities of
THM-3511.  Here no candidate reaches that gate.  Saying “Mealy realization
not reached after the valuation failure” is stronger and more accurate than
reporting a bounded search miss.

## 5. Period-one boundary

At period one,

```text
R(g)=-g^2,                  R^2(g)=-g^4.              (22)
```

The saturated two-cycle equation is `g^3=-1`.  The root `g=-1` is fixed and
repeats gap one.  The two genuine algebraic roots satisfy
`g^2-g+1=0`, whose reduction `x^2+x+1` has no root in `F_2`; they are not
odd-unit `Q_2` profiles.  This recovers THM-3512's stronger scalar-collapse
hostile and explains why period three is the first phase-bearing case.

## 6. Exact controls

Reproduce with

```powershell
python3 -B 04-computation/rule30_two_scale_period3_profile_no_go_20260823.py
python3 -B -O 04-computation/rule30_two_scale_period3_profile_no_go_20260823.py
```

The script checks (4), (9)--(19), the direct factorization (18a), the complete
rational parametrization, the involution `t -> -t`, and the saturated rational
control

```text
t=2:  g=(2/3,-1/2,3),       R(g)=(-2,3/2,1/3),        (23)
```

all three primitive-holonomy eigenvectors, the period-one residue gate, and
the physical parity obstruction on every residue modulo 16.  Normal and
optimized outputs agree; their line-normalized content equals the LF frozen
output.  These finite checks audit the formulas; Sections 1--5 prove the
stated quantifiers.

Raw working-tree byte SHA-256 hashes are

```text
8cb1243fc65f81a6a2f2549d02957975ad86cea93b1ed2242383e8d57e4979fd
  rule30_two_scale_period3_profile_no_go_20260823.py
564eba6e319148a2b826dd267e485292d2c8603941720f52856b6ecb79362fd8
  rule30_two_scale_period3_profile_no_go_20260823.out
```

## 7. Strict scope and next frontier

**Canonically closed by THM-3778:** every exact finite scale-cycle at every
odd spatial period in this cyclic projective-profile ansatz.  Algebraically,
each locus is a finite union of opens in projective eigenspaces, so increasing
the period or cycle length cannot create an elliptic component here.
Physically, the determinant-valuation invoice leaves only constant amplitudes
and consecutive unit gaps, which THM-3512 excludes.

**VERIFIED-EXACT:** this companion's period-three symbolic identities,
rational control, and mod-16 hostile; THM-3778 has an independent general
odd-period companion.

**OPEN:** spatially nonperiodic profiles, nonclosing or infinite scale
evolution, genuinely nonlinear marked sidecars, adaptive overflow states,
actual finite signalizer closure, bounded innovation gaps, center
nonperiodicity, density, balance, and every Rule 30 prize question.

The next lawful algebraic target must leave this linearized finite-cycle
scheme: add a predicate-preserving nonlinear owner/carry sidecar, study a
nonclosing scale evolution, or change the carrier.  Any survivor must still
pass the physical normalization and Mealy-owner gates.
