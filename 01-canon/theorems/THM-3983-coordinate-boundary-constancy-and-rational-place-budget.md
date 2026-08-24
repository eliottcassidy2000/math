---
id: THM-3983
title: "Coordinate boundary constancy and rational-fibre place budget"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On every
  height-tower completion X_n of THM-3973, a source polynomial coordinate
  lying in B_n restricts constantly to the boundary D. More generally, if
  f|D has separable degree d>0 and the geometric generic source fibre is a
  rational curve with r places at infinity, then d<=r-1. In particular a
  cusp restriction y^2 requires at least three places at infinity if its
  generic fibre is rational. The critical-free linear seams of THM-3978 are
  retained as hostile controls: they restrict constantly to D and have
  generic fibre G_m, but are not coordinates. No converse to boundary
  constancy and no classification of rational-fibre elements is claimed.
source: jc-zero-debt-lift + root / post-THM-3981 coordinate-fibre lane, 2026-08-24
audit: >
  PASS (root / jc-cohn3709, 2026-08-24). The audit independently checked
  geometric integrality and reducedness of the completed generic fibre,
  using purity of the Cartier curve and exclusion of a component in D. It
  rederived the d distinct transverse boundary branches, finiteness and
  affineness of normalization, and the exact identity between those added
  branches and removed points on the common P1 model. The final residual
  infinity point, coordinate endpoint, cusp endpoint, and hostile Gm seam
  all check. Normal, optimized, and frozen outputs byte-match at CHECKS=174;
  all hashes agree.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3977-simultaneous-cusp-arm-family-critical-resultant
  - THM-3978-linear-seam-submersion-rational-mate-pole-obstruction
  - THM-3981-centered-cusp-quadrature-genus-two-transcendence
  - THM-3984-boundary-generator-coupling-criticality-and-holomorphic-time-form
script: 04-computation/jc2_coordinate_boundary_place_budget_thm3983.py
output: 05-knowledge/results/jc2_coordinate_boundary_place_budget_thm3983.out
script_sha256: a1b03588747575ed81ee8f23823edcf5e378448e88dfc753f48bda0fb7603687
output_sha256: d0e0fc46ad0be6433f74d09237d80863189b4f5a3f74c5b218c47ad8a95e3483
semantic_sha256: 8f813d0377f5f996625582858894293ba8857f91c45a93471753b64fae48b7c8
hash_basis: raw LF bytes
---

# THM-3983 -- a rational fibre must pay one more infinity place than its boundary degree

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over
an algebraically closed field `k` of characteristic zero. For `n>=2`, use
the smooth affine surface of THM-3973

```text
z=1+x^n t,                  p=zt,
y=x^(n-1)zt^2,
B_n=k[x,z,p,y] subset k[x,t],          X_n=Spec(B_n),       (1)
U=X_n minus D=A2_(x,t),                D=V(x,z)=A1_y.       (2)
```

Let `f in B_n`, and write

```text
g(y)=f|D in k[y].                                           (3)
```

Assume `g` is nonconstant of degree `d`. Suppose that the geometric generic
fibre of

```text
f|U:A2_(x,t) -> A1                                         (4)
```

is geometrically integral and rational. Let `r` be its number of geometric
places at infinity: over an algebraic closure of `k(a)`, its unique smooth
projective completion is `P1`, and the source fibre is `P1` minus `r`
distinct points. Then

```text
d <= r-1.                                                   (5)
```

In characteristic zero the finite map `g:A1_y->A1` is automatically
generically separable, so no extra separability hypothesis is needed for
the displayed family. The proof below records separability explicitly
because it is the precise input that turns degree into distinct boundary
addresses.

Two consequences are immediate.

1. If `f` is a polynomial coordinate of the source ring `k[x,t]`, then
   `g` is constant.
2. If `g(y)=alpha y^2+beta y+gamma` with `alpha!=0`, then a geometrically
   rational generic source fibre has at least three places at infinity.
   In particular every element with boundary cusp profile `y^2` is not a
   source polynomial coordinate.

The theorem is a place-budget obstruction, not a converse. Constant boundary
restriction does not imply coordinatehood, and a nonconstant boundary
restriction can instead pay positive genus rather than the rational place
budget `(5)`.

## 1. The boundary is horizontal and contributes no fibre component

Let `a` be transcendental over `k`, put `K=k(a)`, and then pass to an
algebraic closure `Kbar`. Since `g` has degree `d>0`, the restriction

```text
g:D=A1_y -> A1_a                                            (6)
```

is finite of degree `d`. In particular `D` dominates the base. Thus `D`
cannot be a component of the geometric generic fibre of `f`: scheme
theoretically its intersection with that fibre is

```text
Spec(Kbar[y]/(g(y)-a)).                                    (7)
```

Because `a` is transcendental and `char(k)=0`, the polynomials `g(y)-a`
and `g'(y)` are coprime. Indeed, a common root would be a critical point of
`g`, and its image under `g` would be an element of the algebraically closed
constant field `k`, not the transcendental `a`. Hence `(7)` consists of
exactly `d` distinct geometric points.

This also rules out a hidden copy of `D` in the fibre. Equivalently, in a
local equation `s=0` for `D`, the class of `f-a` modulo `s` is the nonzero
polynomial `g(y)-a`; it is not divisible by `s`.

## 2. Every boundary address is one transverse normalization point

The boundary chart from THM-3973 is

```text
z(z-1)^2=x^(n+1)y,                D=(x,z),                (8)
F_z=(z-1)(3z-1),                  F_z|D=1.                (9)
```

Consequently `x,y` are regular local coordinates near every point of `D`,
`x=0` is a local equation for `D`, and any `f in B_n` has the local form

```text
f=g(y)+x h(x,y)                                          (10)
```

for a regular local function `h`. At a point `y=y_i` of `(7)`, one has
`g'(y_i)!=0`. Therefore the equation `f-a=0` is smooth there and meets
`D` transversely. Each of the `d` points of `(7)` supplies exactly one
branch of the geometric generic fibre closure, and distinct points of `D`
supply distinct normalization points. There is no collision or quotient
loss in passing from the degree of `g` to the number of added places.

For completeness, generic smoothness gives a smooth geometric generic fibre
on `U`. The fibre on the smooth surface `X_n` is an effective Cartier curve,
hence is pure of dimension one. No component can lie in `D` by Section 1,
so every component meets the geometrically integral dense open source fibre.
It follows that the geometric generic fibre on `X_n` is itself integral;
the transverse calculation makes it smooth at the added points. One can
instead normalize it: normalization of an affine finite-type curve over a
field is finite, hence affine, and `(10)` shows that its inverse image over
each boundary address is still one point. Either route gives the same curve.

## 3. The intermediate-open argument on `P1`

Write the geometric generic source fibre as

```text
C_U=P1_Kbar minus S,                  |S|=r.              (11)
```

Let `C_X` be the geometric generic fibre on `X_n`, or equivalently its
normalization as in Section 2. The open immersion `U subset X_n` gives

```text
C_U subset C_X                                       (12)
```

with the same rational function field. By Section 2, the difference
`C_X minus C_U` consists of exactly the `d` transverse boundary points.
The unique smooth projective model of their common function field is
`P1_Kbar`. Hence there is a subset `S_X subset S` such that

```text
C_X=P1_Kbar minus S_X,              |S_X|=r-d.            (13)
```

But `C_X` is affine and one-dimensional. It cannot be the whole proper
curve `P1`, so `S_X` is nonempty. Equation `(13)` therefore gives

```text
r-d>=1,                                                     (14)
```

which is exactly `(5)`. This last `1` is the unavoidable residual infinity
place: an affine completion can absorb all `d` boundary addresses but cannot
absorb the final projective place.

## 4. Coordinates and cusp profiles

If `f` is a coordinate of `k[x,t]`, then the generic fibre of `(4)` is
`A1`, whose projective completion has exactly one missing point. Thus `r=1`.
If `g` were nonconstant, `(5)` would give `d<=0`, a contradiction. Therefore

```text
f a source coordinate and f in B_n   =>   f|D in k.       (15)
```

For a quadratic boundary profile, `d=2`, so `(5)` gives `r>=3`. This is the
precise rational alternative to the genus-two payment in THM-3981: the
present theorem does not assert that a cusp-profile fibre is rational, only
that rationality cannot reduce its infinity support below three.

## 5. Exact hostile controls and sharp scope

The boundary-constancy conclusion is compatible with genuine coordinates.
The element `x in B_n` is a source coordinate, has generic fibre `A1`, and
restricts to zero on `D`.

It is not a coordinate criterion. For `c in k^*`, the THM-3978 linear seam

```text
A_c=x+c(z-1)=x+c x^n t                                  (16)
```

has no affine critical point, because

```text
(A_c)_t=cx^n,
(A_c)_x=1+cnx^(n-1)t,                                   (17)
```

and `(A_c)_t=0` forces `x=0`, where `(A_c)_x=1`. It restricts to the
constant `-c` on `D`. For a generic value `a!=0`, its fibre ring is

```text
K[x,t]/(x+c x^n t-a)=K[x,x^-1],                          (18)
t=(a-x)/(c x^n),
```

so its generic fibre is `G_m`, with two infinity places rather than one.
THM-3978 proves independently that it has no polynomial Jacobian mate.
Thus even a critical-free rational-fibre element with constant boundary
profile need not be a coordinate.

The source coordinate `t` is a second hostile boundary check: it does not
contradict `(15)` because `t` is not regular on `X_n`; its boundary order is
`-n`. Finally, the theorem says nothing about elements with nonrational
generic fibre, about a converse to `(5)`, or about arbitrary affine
completions with singular/nonreduced boundary. Its exact geometric inputs
are the smooth one-line boundary `(8)--(9)`, generically separable
restriction `(6)`, and rationality of the geometric generic source fibre.

## 6. Exact companion and audit state

The companion checks the determinantal family and smooth boundary chart for
heights `2<=n<=10`; the coordinate and critical-free seam boundary rows;
the exact `G_m` generic-fibre parametrization; quadratic cusp restriction;
and separable degree/address and place-budget controls. These computations
are hostile controls for the geometric proof, not a finite census standing
in for it.

Reproduce with

```bash
python3 04-computation/jc2_coordinate_boundary_place_budget_thm3983.py
python3 -O 04-computation/jc2_coordinate_boundary_place_budget_thm3983.py
sha256sum 04-computation/jc2_coordinate_boundary_place_budget_thm3983.py \
  05-knowledge/results/jc2_coordinate_boundary_place_budget_thm3983.out
python3 agents/check_docs.py
```

The geometric generic-fibre passage, transverse normalization count, and
intermediate-open argument have all passed independent hostile audit.
**QED.**
