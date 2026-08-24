---
id: THM-3981
title: "Centered cusp quadrature has a genus-two transcendence obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the generic
  centered boundary-cusp slice used by THM-3979, the canonical quadrature
  curve is birational to v^2=(u^3-Y)^2+2u^2, a smooth genus-two curve over
  k(Y). Its quadrature differential is exactly du/v, a nonzero holomorphic
  form with one simple zero at each infinity. It is not rationally exact,
  and trace descent forbids a primitive after every finite algebraic
  extension. Hence the canonical transverse coordinate X, and already the
  canonical first coordinate A_D=Y^2+2LX, are transcendental over the
  generic D-color function field. On each exceptional scalar slice the
  normalization has genus one and the form has two opposite nonzero
  residues, so canonical nonalgebraization survives every scalar
  specialization. This does not cover every alternative formal Darboux
  gauge; unrestricted Darboux pairs and JC(2) remain open.
source: root / jc-extra-debt-local / post-THM-3979 generic-slice audit, 2026-08-24
audit: >
  PASS (root probe-3981-genus2 and audit-3979-formal, 2026-08-24). Two
  independent derivations checked the centered projection, birational
  inverse, sextic discriminant, generic genus, differential sign, infinity
  divisor, and trace-descent obstruction. The optimization-safe companion
  freezes the two exceptional discriminant slices and the first formal
  D-chart coefficients. Normal, optimized, and frozen outputs byte-match at
  CHECKS=25; both raw hashes and the semantic hash agree.
depends_on:
  - THM-3979-two-color-formal-cusp-darboux-lifting
related:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3980-all-height-canonical-split-formal-cusp-nonalgebraization
script: 04-computation/jc2_centered_cusp_quadrature_genus_two_thm3981.py
output: 05-knowledge/results/jc2_centered_cusp_quadrature_genus_two_thm3981.out
script_sha256: 9ba76e156c092410e755d2680e5c7413721d729957802bd03da1dc7ef20ff3b8
output_sha256: 3f29a9da90581085fae01be3fe8887bdf294f5a1100f247b0278c179ac765b32
semantic_sha256: f88ec164728e05717827b54c573afed10668a4c4c448aae3f27f0bd8bcbb15b7
hash_basis: raw LF bytes
---

# THM-3981 -- the canonical cusp quadrature is transcendental

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. In the height-two
completion of THM-3973, the boundary-color chart used by THM-3979 is

```text
z(z-1)^2=x^3y,
w=1/((z-1)(3z-1)),
dx wedge dt=xw dx wedge dy.                              (1)
```

Center the cusp coordinate by `Y=y-x/2`, and regard `Y` as transcendental.
Over `K=k(Y)`, put

```text
F=K(x,z),             z(z-1)^2=x^3(Y+x/2).              (2)
```

THM-3979's canonical transverse coordinate is the unique formal series
`X=x+O(x^2)` defined by

```text
X^2=2 integral_0^x s w(s,Y+s/2) ds.                     (3)
```

The following three statements hold.

1. The curve with function field `F` is birational over `K` to

   ```text
   v^2=(u^3-Y)^2+2u^2,                                  (4)
   ```

   a smooth genus-two curve.

2. Under this birational map, the quadrature form is

   ```text
   xw dx=du/v.                                           (5)
   ```

   The right side is a nonzero holomorphic differential whose divisor is
   the sum of the two points at infinity.

3. No element in any finite algebraic extension of `F` has differential
   `du/v`. Consequently `X` is transcendental over `F`. For any
   `L in k^*`, the canonical cusp first coordinate

   ```text
   A_D=Y^2+2LX                                           (6)
   ```

   is transcendental as well.

4. For every scalar specialization `Y=Y_0`, the specialized canonical
   quadrature remains nonalgebraic. Away from the discriminant locus the
   same genus-two argument applies; on its two exceptional loci, the
   normalized form has opposite nonzero residues.

Thus the exact formal construction of THM-3979 cannot terminate or
algebraize in its canonical quadrature gauge. The conclusion is deliberately
gauge-specific: it does not exclude a different formal Darboux coordinate,
a different first coordinate in `B_2`, or a polynomial Keller pair.

## 1. The centered curve is generically genus two

In the function field `F`, put

```text
u=(z-1)/x.                                               (7)
```

Substitution of `z=1+xu` into `(2)` and cancellation in the function
field give

```text
x^2+2(Y-u^3)x-2u^2=0.                                   (8)
```

Complete the square with

```text
v=x+Y-u^3.                                               (9)
```

Then `(8)` becomes `(4)`. Conversely,

```text
x=v-Y+u^3,                  z=1+u(v-Y+u^3),             (10)
```

so `(7)--(10)` identify the two function fields rather than merely
constructing a dominant map.

The sextic on the right of `(4)` has exact discriminant

```text
Disc_u((u^3-Y)^2+2u^2)
  =-512Y^2(729Y^4+128).                                 (11)
```

This is nonzero in `K`. Hence the sextic is squarefree of degree six, and
the smooth projective model of `(4)` is hyperelliptic of genus two. It has
six finite branch points and two unramified points at infinity.

## 2. The quadrature is the first holomorphic differential

Differentiate `v^2=(u^3-Y)^2+2u^2` at fixed `Y`. Using
`x=v-Y+u^3` gives

```text
v dx/du=u(3ux+2).                                       (12)
```

On the other hand, `z-1=xu` and `3z-1=2+3xu`. Therefore

```text
xw dx
 =dx/[u(2+3xu)]
 =du/v,                                                  (13)
```

which proves `(5)` including its sign.

At a finite branch point, `v` is a local parameter and `du/v` is regular.
At either infinity, take `s=1/u`; since `v/u^3` tends to `+1` or `-1`,

```text
du/v = -/+ s ds + higher terms.                         (14)
```

Thus the form has one simple zero at each infinity and no other zero or
pole. This is the canonical degree-two divisor of a nonzero holomorphic
differential on the genus-two curve.

## 3. Trace descent forbids algebraic quadrature

A nonzero holomorphic differential on a smooth projective curve cannot be
the differential of a rational function in characteristic zero. Indeed, a
pole of order `m` in a rational function produces a pole of order `m+1`
in its differential; a function with holomorphic differential therefore has
no poles and is constant.

The obstruction persists after algebraic extension. Let `E/F` be finite.
Characteristic zero makes it separable, and for every `H in E`,

```text
d Tr_(E/F)(H)=Tr_(E/F)(dH).                              (15)
```

If `dH=du/v`, the right side is `[E:F]du/v`, so division by the nonzero
integer `[E:F]` would make `du/v` exact already in `F`, contradicting
the preceding paragraph.

Finally, differentiating `(3)` at fixed `Y` and using `(5)` yields

```text
d(X^2)=2xw dx=2du/v.                                    (16)
```

If `X` were algebraic over `F`, then `X^2` would lie in a finite
extension and violate `(15)--(16)`. This proves the transcendence of
`X`; equation `(6)` then proves the same for `A_D`.

## 4. Failure boundary and scope

The discriminant `(11)` vanishes exactly at

```text
Y=0                  or                  729Y^4+128=0.    (17)
```

At `Y=0`, write `v=uV`; the normalization is

```text
V^2=u^4+2,
du/v=du/(uV).                                            (18)
```

It has genus one, and the form has simple poles at the two points over
`u=0`, with residues whose squares are `1/2`. On the other exceptional
locus the unique repeated root is `a=4/(9Y)`, with

```text
R''(a)/2=-4.                                             (19)
```

Writing `v=(u-a)V` again gives a genus-one normalization, now with two
simple residues whose squares are `-1/4`. In both cases the nonzero residues
forbid rational exactness; trace descent forbids exactness after every finite
algebraic extension. Thus the canonical obstruction survives every scalar
slice, while `(11)` marks the change from the generic holomorphic mechanism
to a logarithmic one.

There is no contradiction with THM-3979. Formal power-series completion can
contain elements transcendental over the algebraic function field, and
THM-3979 asserts formal solvability, not algebraicity. The present theorem
turns that distinction into an exact obstruction for its chosen quadrature.
THM-3980 remains reserved for all-height or gauge-invariant
nonalgebraization. Alternative Darboux pairs and `JC(2)` remain open.

## 5. Exact companion

The assertion-free companion checks the birational equations, discriminant
and resultant, both degeneration controls, differential identity, and the
first formal D-chart terms with 25 optimization-safe gates. Reproduce with

```bash
python3 04-computation/jc2_centered_cusp_quadrature_genus_two_thm3981.py
python3 -O 04-computation/jc2_centered_cusp_quadrature_genus_two_thm3981.py
```

Both runs byte-match the frozen output after LF normalization. **QED.**
