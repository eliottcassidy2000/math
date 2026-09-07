# A finite quartic map and a global relative-primitive obstruction

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
September 6, 2026. This concerns the single explicit DG surface and maps
with a prescribed first-coordinate class. It is not a general Keller
exclusion or an identification of a possible finite envelope.

## 1. Object and inherited distinction

Use the smooth affine surface and charts in
[the DG surface note](planar_jc48_sep06_dg_surface.md):

```text
W=(P1_s x P1_z)\{z=s^2},
U0=A2_(x,t),       Uinfty=A2_(r,b),
x=1/r,            t=-r^2-r^4b,
b=-x^2-x^4t,      D={r=0},      W\D=U0.
```

It has `Pic(W)=Z[D]`, `K_W=2[D]`, and `chi_c(W)=2`. Its regular form
`omega=dx wedge dt=r^2 dr wedge db` is globally exact. The inherited
carrier `k[u,p,y]`, with `u=x^2t,p=t(1+u),y=xtp`, collapses `D`, while
the global function `b` separates it. This addendum turns that separator
into an actual finite map, then identifies exactly what fails.

The board is: boundary separation; proper finite morphisms; ramification
inside the source; exact global forms; and a prescribed first coordinate.
The corrected near miss is to treat separation of the boundary as enough
for an etale source map. The restoring sidecar is the divisor of the
first coordinate, together with regularity of a relative primitive.
The surface itself is classical, as cited in the supplier; no external
priority claim is made for the elementary calculations below.

## 2. Explicit finite flat degree-four morphism

Both `t` and `b` are globally regular. Define

```text
f:W -> A2_(T,B),        (T,B)=(t,b).                      (1)
```

In fact `W` is exactly the closed subscheme

```text
T X^4+X^2 R^2+B R^4=0
       in P1_[X:R] x A2_(T,B).                           (2)
```

On `R!=0`, put `x=X/R`; the equation solves
`B=-x^2-Tx^4`, giving `U0`. On `X!=0`, put `r=R/X`; it solves
`T=-r^2-Br^4`, giving `Uinfty`. The overlaps are precisely the inherited
transition. This verifies an isomorphism, not just a birational comparison.

The projection is projective. Each fibre is a nonzero binary quartic,
because its `X^2R^2` coefficient is one, so every fibre has finite length
four. Hence the projection is proper and quasi-finite, therefore finite.
The following free basis also proves flatness and degree four directly.

Put `A0=C[T,B]`. The global functions

```text
u=x^2t,      j=xt,      v=x(1+x^2t)
```

have boundary-chart expressions
`-1-r^2b`, `-r(1+r^2b)`, and `-rb`, respectively. Then

```text
O(W)=A0*1 direct_sum A0*u direct_sum A0*j direct_sum A0*v, (3)
u^2=-u-TB,       j^2=Tu,       v^2=-B(1+u),
uj=Tv-j,        uv=-Bj,       jv=-TB.                    (4)
```

Here is a check that (3) is the full integral module, rather than only
a generic-field basis. Push forward the exact sequence of the quartic
divisor in the projective line bundle:
`0->O(-4)->O->O_W->0`. The usual two-chart Laurent calculation gives
`H^1(P1_A0,O(-4))=A0*x^-1+A0*x^-2+A0*x^-3` and `H^1(O)=0`.
Writing `F=Tx^4+x^2+B`, the differences of the two regular lifts are

```text
u: Tx^2-(-1-B/x^2)=F/x^2,
j: Tx-(-1/x-B/x^3)=F/x^3,
v: x+Tx^3-(-B/x)=F/x.
```

Their connecting classes form that complete Laurent basis. Together
with the constant function they give (3). The multiplication identities
(4) are direct in either original chart.

## 3. Full discriminant and exact ramification

The trace pairing in the basis `(1,u,j,v)` is

```text
[ 4    -2          0       0    ]
[-2   2-4TB        0       0    ]
[ 0     0         -2T     -4TB  ]
[ 0     0         -4TB    -2B   ].
```

Its determinant, independently equal to the binary quartic discriminant,
is

```text
disc(f)=16TB(1-4TB)^2.                                   (5)
```

Thus the branch locus is exactly the two axes and the hyperbola `4TB=1`.
This includes `T=0`, which would be missed by treating the affine
degree drop in `Tx^4+x^2+B` as loss of a source point.

The two literal Jacobians are

```text
dt wedge db=2x(1+2x^2t) dx wedge dt          on U0,
dt wedge db=-2r(1+2r^2b) dr wedge db         on Uinfty.   (6)
```

Let `Z={x=0}` and `C={1+2x^2t=0}` inside `U0`. They are closed in
`W`, are disjoint from `D`, and all three curves are smooth. Equation
(6) proves that the ramification divisor is exactly

```text
Ram(f)=D+Z+C,
```

with generic ramification index two along each. There are no additional
components: the two charts cover `W` and (6) gives their entire zero loci.
The complete generic prime profiles, with `(e,f_res)` meaning ramification
index and residue degree, are

| Target branch | Ramified prime | Other prime, if present |
|---|---|---|
| `T=0` | `D: (2,1)` | `L: (1,2)`, where `L={t=0}` in `U0` |
| `B=0` | `Z: (2,1)` | closure of `{1+x^2t=0}: (1,2)` |
| `4TB=1` | `C: (2,2)` | none |

For the first row `b|L=-x^2`; for the second, `t=-r^2` on the other
component `b=0` in the second chart. For the third,
`t=-1/(2x^2)` and `b=-x^2/2`, so the residue map has degree two.
The identity `1-4tb=(1+2u)^2` verifies the square pullback in (5).
The origin has two distinct fibre points, each of length two, as the
binary quartic becomes `X^2R^2`. In particular (1) is a real finite
surface map, but it is ramified along `Z` and `C` in the retained affine
plane and is not a Keller map.

These divisors are compatible with the canonical invoice. As rational
divisors `div(x)=Z-D` and `div(1+2u)=C`; using `div(omega)=2D` gives
`div(dt wedge db)=Z+D+C` and class `2[D]`. Matching that class does
not eliminate the forbidden interior ramification.

## 4. An unbounded obstruction with first coordinate in `C[t]`

The geometric supplier of the obstruction is

```text
div_W(t)=L+2D.                                           (7)
```

The two components are disjoint affine lines. At a generic point of
`L`, `t` is the original simple coordinate. At `D`, its expression
`-r^2(1+r^2b)` has order exactly two. There are no other zero components.

**Relative-primitive theorem.** There is no global regular function
`g in O(W)` and nonzero constant `c` such that

```text
dt wedge dg=c omega.                                    (8)
```

Indeed restriction to `U0` forces `-g_x=c`. Polynomial integration in
characteristic zero gives `g=-cx+h(t)` for a polynomial `h`. But on
the second chart this is

```text
-c/r+h(-r^2-r^4b),
```

whose simple pole cannot cancel. It is not a global regular function.
Thus `omega` has a regular global primitive `x dt`, while no regular
primitive of the prescribed relative form (8) exists. The globally
regular product `x dt` does not make its coefficient `x` regular.
This directly connects exactness with the earlier descent/differential
questions while retaining the missing regularity coordinate.

More generally, for every nonconstant polynomial `F` and global `g`,
the source Jacobian of `(F(t),g)` is `-F'(t)g_x`. If it were a nonzero
constant, both factors would be units of `C[x,t]`; `F` would be linear
and (8) would follow after rescaling. Hence that Jacobian is never a
nonzero constant. In particular **every finite map `(F(t),g):W->A2`
ramifies somewhere in `U0`**: its Jacobian is a nonzero polynomial by
generic separability, and a nonunit polynomial over `C` has a zero.

There is also an exact boundary index consequence. A finite map with
these two coordinates cannot make `g|D` constant, since that would
collapse `D`. It maps `D` onto the vertical target line with first
coordinate `F(0)`. If `m=ord_(t=0)(F(t)-F(0))`, then (7) gives

```text
e_D=ord_D(F(t)-F(0))=2m.                                 (9)
```

Thus the boundary index is even throughout this fixed-first-coordinate
class; it cannot supply the index-three prime of the alternating
envelope reduction. This is a scoped exclusion of these coordinate
choices, not a classification of all finite maps from `W`.

## 5. Scope and exact reproduction

The source is the explicit two-chart surface, the target is an actual
finite quartic morphism and a prescribed-coordinate obstruction. The
map keeps the boundary separator and all projective fibre points. The
ramification calculation retains the interior divisor that the bare
Picard/canonical/Euler data forget. The missing coordinate in the relative
primitive problem is the pole of `x`, despite regularity of `x dt`.

The positive control (1) proves that finite maps exist outside the
boundary-collapsing carrier. The strongest negative conclusion is (8)–(9)
for first coordinates polynomial in `t`. No arbitrary coordinate pair,
finite envelope identification, or general Jacobian conclusion follows.

The [standalone source](../../04-computation/planar_jc48_sep06_dg_finite_map.py)
and [output](planar_jc48_sep06_dg_finite_map.out) check the actual two
quartic charts, their nonzero fibre equation, all module relations,
independent multiplication-matrix traces and binary discriminant, the
three Laurent connecting classes, both Jacobians, branch profiles, and
the exact uncancelled relative-primitive pole. The all-degree statements
use the proofs above; no parameter census replaces them.

```sh
python3 -B 04-computation/planar_jc48_sep06_dg_finite_map.py
python3 -B -O 04-computation/planar_jc48_sep06_dg_finite_map.py
```

Both ordinary and optimized runs pass **41 exact gates** and agree
byte for byte with the frozen output. Replay pins:

```text
source SHA256: e1f4474f580b9cd076f1c00ebdb34aa9451226ed8716a74e044bf5775e31315f
output SHA256: f12e3a3f7d0b558fef894aab4d3addf3284b25ac2b0559e0014e08c4409c6b62
semantic SHA256: a2d4e154df72764f0e53027b12eb6587800c416ccc5357ce2acc47475f1c89cc
```

The [independent analytic and source audit](planar_jc48_sep06_dg_finite_map_audit.md)
passes. Source and output are frozen; the main DG bundle is unchanged
by this addendum.
