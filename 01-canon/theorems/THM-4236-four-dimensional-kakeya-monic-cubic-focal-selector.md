---
id: THM-4236
title: "Four-dimensional Kakeya monic-cubic focal selector"
status: >
  PROVED Euclidean line-selector theorem.  A Lipschitz affine-intercept
  selector on a positive-measure three-dimensional slope patch sweeps positive
  four-volume.  The absolute line Jacobian is governed by a monic focal cubic,
  with sharp full-segment L1 constant 1/64, an arbitrary-shading bound, a
  two-ends linear-density bound, and a polynomial-atlas Bezout volume bound.
  This constrains possible counterexamples but does not prove the
  four-dimensional Kakeya conjecture or improve the general
  published/preprint exponents.  The selector, Sobolev exclusion, and sharp
  Chebyshev L1 mechanism also extend to every dimension n>=2.
source: codex-kakeya4d-broadness-multiscale-20260826 + kakeya_analytic_model + selector_lemma_referee
related:
  - THM-4035-sixty-clock-separation-and-finite-kakeya-spine
  - THM-4235-finite-four-dimensional-torus-kakeya-quadratic-carrier
  - HYP-2235-finite-field-kakeya-falconer-carrier
source_reference: 05-knowledge/reference/CORE-PAPERS-KAKEYA-4D-2026-08-24.md
---

# THM-4236 -- four-dimensional Kakeya monic-cubic focal selector

**PROVED.**  The native object is a section of affine line space, not a list
of directions.  In four dimensions a monic cubic focal determinant governs
the incidence projection's absolute Jacobian along each line.  This gives an
exact regularity obstruction for a measure-zero Kakeya construction and a
quantitative finite-complexity atlas bound.

## 1. Affine-intercept normal form

Let `U subset B_R(0) subset R^3` be measurable with
`0<|U|<infinity`.  The affine line of slope `[1:u]` has a unique
intersection `(0,a(u))` with the hyperplane `x_0=0`.  Write

```text
F(u,t)=(t,a(u)+t u).                                  (1)
```

Assume:

1. `a:U->R^3` is the restriction of a Lipschitz map;
2. `I_u subset [-T,T]` is a measurable interval and
   `D={(u,t):u in U,t in I_u}` is measurable; and
3. `F(D) subset E`.

For selected unit segments,

```text
L(u):=|I_u|=(1+|u|^2)^(-1/2)
             >= ell_R:=(1+R^2)^(-1/2).                (2)
```

The position of `I_u` along the line may be merely measurable.  Only the
transverse line intercept `a(u)` is required to be Lipschitz.  If an endpoint
is `b=(b_0,b')`, then

```text
a(u)=b'(u)-b_0(u)u,                                   (3)
```

which is invariant under shifting `b` along the line.  An arbitrary
unoriented basepoint without its permitted parameter interval is not enough
to ensure containment in `E`.

At almost every differentiability point of `a`, define the focal determinant

```text
P_F(u,t)=det(Da(u)+t I_3).                            (4)
```

This is a monic cubic in `t`.  With domain coordinates ordered as `(u,t)`,
`det DF=-P_F`; only the absolute value enters below.  For a merely measurable
`U`, take a Lipschitz extension and use approximate differentiation at density
points.  The area formula gives

```text
integral_(R^4) N_F(y) dy
  = integral_U integral_(I_u) |P_F(u,t)| dt du,        (5)
```

where `N_F` counts preimages in `D`.

## 2. Sharp cubic L1 law

For every monic real cubic `p` and every interval `I` of length `L`,

```text
integral_I |p(t)| dt >= L^4/64.                       (6)
```

The constant is sharp.  On `[-1,1]`, put

```text
P_*(x)=x^3-x/2=2^(-3) U_3(x),qquad h=sign(P_*).
```

The function `h` is odd, so it annihilates `1` and `x^2`.  Its positive
root is `c=1/sqrt(2)`, and

```text
integral_(-1)^1 x h(x) dx
  =2(-integral_0^c x dx+integral_c^1 x dx)
  =1-2c^2=0.                                         (7)
```

Thus `h` annihilates every quadratic.  For any monic cubic `P`,

```text
integral |P| >= integral hP
              = integral hP_*
              = integral |P_*|=1/4.                  (8)
```

Translation and the change of scale from an interval of length two to one of
length `L` multiply the integral by `L^4/16`, proving `(6)`.  Equality on
`[0,1]` is attained by

```text
q_*(t)=t^3-(3/2)t^2+(5/8)t-1/16,                     (9)
```

whose roots are

```text
(1-1/sqrt(2))/2,  1/2,  (1+1/sqrt(2))/2.
```

Taking `a(u)=Au` with `A` diagonal and eigenvalues equal to the negatives of
these roots realizes `q_*(t)=det(A+tI_3)`.  Sharpness therefore occurs inside
the affine-line-selector class, not merely among abstract cubics.

### Dimension-n extension

The cubic calculation is the `m=3` instance of a sharp general law.  For
every integer `m>=1`, every monic real polynomial `p` of degree `m`, and every
interval `I` of length `L`,

```text
integral_I |p(t)|dt >= L^(m+1)/2^(2m).               (9a)
```

On `[-1,1]`, the extremizer is `2^(-m)U_m`, where `U_m` is the Chebyshev
polynomial of the second kind.  It is monic and

```text
integral_(-1)^1 |2^(-m)U_m(x)|dx=2^(1-m).           (9b)
```

Indeed, under `x=cos(theta)`, the dual sign
`sigma=sign(U_m)=sign(sin((m+1)theta))` annihilates every polynomial of degree
less than `m`: after multiplication by `sin(theta)`, such a polynomial is a
linear combination of `sin(k theta)`, `1<=k<=m`, whereas the square wave
`sigma` has only odd-multiple frequencies of `m+1`.  The same dual argument
as `(8)` proves `(9a)`.  On an interval of length `L` centered at `c`, equality
is attained by

```text
(L/4)^m U_m(2(t-c)/L).                              (9c)
```

Consequently, in `R^n`, with `u in R^(n-1)` and
`F(u,t)=(t,a(u)+tu)`, the focal determinant `det(Da+tI_m)` is monic of degree
`m=n-1`, while `det DF=(-1)^m det(Da+tI_m)`.  A Lipschitz affine-intercept
selector on a positive-measure slope patch again sweeps positive `n`-volume.
A coordinate-degree-`d` polynomial patch has almost-everywhere multiplicity
at most `max(d,1)^(n-1)` by the same regular-value and affine isolated-zero
Bezout argument.  The four-dimensional constant in `(6)` is `(9a)` with
`m=3`.

Combining `(5)` and `(6)` gives

```text
integral N_F
  >= (1/64) integral_U L(u)^4 du
  >= ell_R^4 |U|/64.                                (10)
```

In particular, `F(D)` has positive four-dimensional measure.  No
multiplicity bound is needed for this qualitative conclusion: a positive
integral in `(5)` cannot be supported on a null image.

Consequently, any selected family inside a null four-dimensional Kakeya set
has no positive-measure slope patch on which its invariant affine-intercept
section is Lipschitz.

### Sobolev exclusion

There is a stronger qualitative corollary in every `R^n`, `n>=2`.  Let
`Omega subset R^(n-1)` be open, let a measurable selected family of unit
segments lie in a null set `E`, and assume the segment parameter intervals
have measurable finite endpoints.  On no nonempty open `V subset Omega` can
the affine-intercept section agree almost everywhere with a member of

```text
W^(1,p)_loc(V;R^(n-1)),qquad 1<=p<=infinity.         (10a)
```

Indeed, pass to a ball compactly contained in `V`, apply a cutoff, and extend
the Sobolev representative to `g in W^(1,1)(R^(n-1))`.  For its precise
Lebesgue representative, the standard maximal-gradient inequality gives,
off a null set,

```text
|g(x)-g(y)|
 <= C_n |x-y|[M|nabla g|(x)+M|nabla g|(y)].          (10b)
```

The maximal function is finite almost everywhere, including at the endpoint
`p=1`.  Hence the level sets of `M|nabla g|` give countably many Lipschitz
restrictions covering almost every slope.  Partition once more by bounded
segment-parameter location if needed.  Some cell has positive measure, and a
coordinatewise Lipschitz extension puts it under the dimension-`n` selector
theorem, forcing positive `n`-volume inside `E`, a contradiction.

This argument uses Lipschitz pieces, not a global Sobolev area formula; that
distinction is necessary at `p=1`.  It supplies no quantitative summability
of atlas complexity and does not automatically extend to `BV` or fractional
Sobolev sections.

## 3. Essential multiplicity and polynomial atlases

If

```text
N_F(y)<=M for almost every y,                         (11)
```

then `(10)` gives the quantitative bound

```text
|F(D)| >= [integral_U L(u)^4 du]/(64M)
       >= ell_R^4 |U|/(64M).                         (12)
```

The qualifier “almost every” is load-bearing.  The concurrent star through
`p=(c,p')` has

```text
a(u)=p'-cu,qquad P_F(u,t)=(t-c)^3.                 (13)
```

The full Jacobian has the harmless orientation sign.  The map is one-to-one
almost everywhere but has an infinite fibre at
the critical apex.

For a finite polynomial atlas, partition

```text
U=disjoint_union_(j=1)^K U_j
```

and suppose `a|_(U_j)=a_j`, where each coordinate of `a_j` is a polynomial
of total degree at most `d_j`.  Put `D_j=max(d_j,1)`.  A preimage of
`y=(y_0,y')` satisfies

```text
t=y_0,qquad a_j(u)+y_0u-y'=0.                       (14)
```

These are three equations of degree at most `D_j`.  At a regular preimage
their derivative is `Da_j(u)+y_0I_3`, so the complex zero is isolated and
simple.  Sard's theorem removes the critical values, and affine isolated-zero
Bezout gives

```text
N_F(y)<=sum_j D_j^3 for almost every y.               (15)
```

Therefore

```text
|F(D)| >=
  [integral_U L(u)^4 du]/[64 sum_j D_j^3].            (16)
```

The cubic complexity is real: coordinatewise Chebyshev maps of degree `d`
have open target boxes with `d^3` simple real inverse branches.  If polynomial
degree is assigned to an endpoint `b` rather than the invariant intercept
`a=b'-b_0u`, the safe degree is `d+1`, not `d`.

## 4. Approximate atlases and a multiscale complexity debt

Suppose instead that on `U_j`

```text
|a(u)-a_j(u)|<=delta,                                 (17)
```

with `a_j` polynomial of degree `d_j`.  Use the same interval `I_u` and
let `F_j(u,t)=(t,a_j(u)+tu)`.  For the piecewise polynomial model
`F_poly=F_j` on `D|_(U_j)`, one has

```text
F_poly(D) subset N_delta(F(D)) subset N_delta(E).
```

Applying `(16)` to the polynomial model proves

```text
|N_delta(E)| >=
  [integral_U L(u)^4 du]/
  [64 C_delta],qquad
C_delta:=sum_j max(d_j,1)^3.                          (18)
```

Thus a small-neighborhood counterexample must pay an explicit affine-atlas
complexity debt.  For example, if `|N_delta(E)|<=C delta^(4-s)`, then every
uniform `delta`-accurate atlas of this type satisfies

```text
C_delta >= c_(R,U,C) delta^(-(4-s)).                  (19)
```

This is a counterexample-space squeeze, not an existence theorem for a
low-complexity atlas.

## 5. Arbitrary shadings and the focal divisor

Let `Y_u subset I_u` be measurable, assume the joint shaded domain
`{(u,t):u in U,t in Y_u}` is measurable, and write
`lambda(u)=|Y_u|`.  If a monic cubic factors over `C` as
`p(t)=product_(j=1)^3(t-z_j)`, then

```text
{|p|<epsilon}
 subset union_j {|t-z_j|<epsilon^(1/3)}.              (20)
```

Its real length is at most `6 epsilon^(1/3)`.  The increasing rearrangement
of `|p|` on any measurable set of length `lambda` is therefore bounded below
by `(s/6)^3`.  Integration gives

```text
integral_(Y_u)|P_F(u,t)|dt >= lambda(u)^4/864.        (21)
```

Hence

```text
integral N_Y >= (1/864) integral_U lambda(u)^4 du.    (22)
```

Under essential multiplicity `M`, divide the right side by `M` to lower-bound
the shaded sweep volume.  For a polynomial atlas, divide by `sum_j D_j^3`.

Equation `(20)` is the focal law: a cubic line selector has at most three
moving focal times per direction.  An arbitrary shading can compress near
those times, explaining the quartic density loss.  The concurrent star
`(13)` with shading `[c,c+lambda]` has Jacobian mass `lambda^4/4`, so the
quartic order cannot be replaced without a distribution hypothesis.

## 6. Two-ends converts quartic density to linear density

Assume `A>=1`, `alpha>0`, and for every translated interval `J`

```text
|Y_u intersect J|
 <= A (|J|/L(u))^alpha |Y_u|.                        (23)
```

The cubic sublevel set at radius `r` is covered by at most three intervals of
length `2r`.  Choose

```text
r=[L(u)/2](6A)^(-1/alpha).                            (24)
```

By `(23)`, these intervals capture at most half the shading.  On the remaining
half, `|P_F|>=r^3`, so

```text
integral_(Y_u)|P_F|dt
 >= [L(u)^3/16](6A)^(-3/alpha)|Y_u|.                 (25)
```

Consequently

```text
integral N_Y
 >= [(6A)^(-3/alpha)/16]
    integral_U L(u)^3 |Y_u| du.                      (26)
```

The same essential-multiplicity or polynomial-atlas division converts `(26)`
to a volume bound.  Control only for midpoint-centered intervals is
insufficient; all translated intervals at the tested scale are required.

This is the precise interface with modern two-ends incidence estimates:
two-ends prevents a shading from hiding near the three-point focal divisor and
changes the density dependence from quartic to linear.  It does not itself
control the remaining affine projection multiplicity.

## 7. Relation to the finite carrier and the open frontier

THM-4235 uses the polynomial intercept

```text
a(s)=(s_1^2,s_2^2,s_3^2).
```

Then

```text
F(s,t)=(t,s_1^2+ts_1,s_2^2+ts_2,s_3^2+ts_3),
P_F=product_i(2s_i+t).                               (27)
```

The degree-two Bezout bound is `2^3=8`, exactly the maximum multiplicity of
the finite torus carrier.  Its four boundary charts have actual maximum
multiplicity twelve, below the crude four-chart bound thirty-two.  Reduction
modulo `61` preserves polynomial degree and the threefold fold, but not
Euclidean angles, tube volume, shadings, two-ends, or scale.

The source-target contract is

```text
affine-line section u |-> (u,a(u))
  -> incidence manifold (u,t)
  -> physical union F(u,t).
```

It preserves direction, transverse placement, and shading.  Projection loses
multiplicity and the direction flag of each fibre.  The necessary sidecars are
therefore:

- essential or polynomial-atlas multiplicity;
- two-ends shading distribution;
- a broad/narrow flag on the incident directions; and
- multiscale seam/algebraic-parent data.

THM-4035's one-clock twisted cubic has finite quartet full-spark but only a
one-dimensional parameter.  It lacks the three-dimensional direction form
that makes `(4)` monic.  Its three-clock torus has the right parameter count
only after discretization and still lacked `a(u)` before THM-4235.  This
explains the earlier guardrail rather than weakening it.

A complete Euclidean attack would need an atlas/seam dichotomy: either
`C_delta` is small enough for `(18)`, many patch collisions are transverse,
or the narrow direction flag varies slowly enough for a planebrush.  This
dichotomy is **OPEN**.  THM-4236 proves no new general Hausdorff or maximal
exponent, and Euclidean Kakeya in `R^4` remains **OPEN**.  **QED.**
