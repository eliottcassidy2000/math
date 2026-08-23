---
id: THM-3832
title: "The nonlinear cubic surface has a triangular root-ratio chart"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  On the dense k!=0 chart of the THM-3811 etale surface, z=h/k=A/omega
  and C are birational coordinates.  In them the target map is triangular,
  A=zC(1+z^2C)/(3z^3+7z^2+1), and its Jacobian density is exactly 1/(kr).
  All intrinsic generators and the precise polynomialization/Keller passport
  are explicit.  This constructs a rational chart, not a polynomial plane
  atlas or a Jacobian counterexample.
source: root / nonlinear-cubic plane-atlas constructive reframe, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The deterministic companion reconstructs the
  three Delone--Faddeev multiplication laws, characteristic cubic, different,
  SL2 determinant and both lift laws from the chart; verifies both birational
  directions, the factored target derivative, the source-chain Jacobian law,
  and the two cancellation branches over every cubic spectral root.  Normal
  and optimized runs byte-match the frozen transcript; independent hostile
  audit remains pending.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
related:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
script: 04-computation/jc2_nonlinear_cubic_root_ratio_chart_thm3832.py
output: 05-knowledge/results/jc2_nonlinear_cubic_root_ratio_chart_thm3832.out
script_sha256: 3e21279f5d1deac6960d6d3a7b17d9b8d7accfb9428120b518627192737c5460
output_sha256: c3aafc9a27ae4daea7c82ec9811a8dcfdeee6bec56093857d225092be7846c52
semantic_sha256: a5445ccd279daadc45016d001be9b28233c0c8c306e11b2e10a752e9a1179853
hash_basis: raw LF bytes
---

# THM-3832 -- the marked cubic root makes the target map triangular

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**  Work over an algebraically closed field `K` of characteristic zero.
Let `U=Spec B` be the THM-3811 surface, with intrinsic functions

```text
A,C,omega,theta,D,       h=A/D,       k=omega/D,
m=3theta+14A.                                               (1)
```

On the dense chart `k!=0`, put

```text
z=h/k=A/omega.                                               (2)
```

Define the four one-variable polynomials

```text
r=3z^3+7z^2+1,                 q=7z^2+3,
b=6z^3+7z^2-1,                 s=z^2 q C-b.                  (3)
```

Then `K(U)=K(z,C)`, and every intrinsic function is the following rational
function of `(z,C)`:

```text
A     = z C(1+z^2 C)/r,
omega =   C(1+z^2 C)/r,
theta = -z C(7Cz^2+C-3z)/r,

D = C^2(1+z^2 C)s/r^2,
k = r/(Cs),                         h = zr/(Cs),
m = -zC(7Cz^2+3C-9z-14)/r.                              (4)
```

The target morphism is therefore triangular in these birational coordinates:

```text
(z,C) |--> (A(z,C),C).                                     (5)
```

Its exact Jacobian density is

```text
partial A/partial z = C s/r^2 = 1/(k r).                    (6)
```

Consequently, for any rational source functions `z=z(x,y)` and `C=C(x,y)`
landing in this chart,

```text
Jac_(x,y)(A,C)=Jac_(x,y)(z,C)/(k r).                         (7)
```

Thus the Keller equation with value `lambda in K*` becomes the single
weighted area equation

```text
Jac_(x,y)(z,C)=lambda k r.                                  (8)
```

Equation `(8)` is an equality of rational functions.  It does not erase the
separate requirement that the expressions `(4)` extend polynomially on the
whole source plane.

## 1. Derivation from the marked-root cubic

The characteristic polynomial of `omega` in THM-3811 is

```text
omega^3-C omega^2+7A^2 omega+3A^3-A^2C^2=0.                 (9)
```

Substitute `omega=A/z`, divide in the function field, and multiply by `z^3`.
The result is

```text
A(3z^3+7z^2+1)=zC(1+z^2C),                                 (10)
```

which is the first formula in `(4)`.  It also gives `omega=A/z`.

The first cubic multiplication law gives

```text
theta=(C omega-7A^2-omega^2)/A.                             (11)
```

Moreover the intrinsic different is the derivative of `(9)` at its marked
root:

```text
D=C omega-3A theta-14A^2
 =3omega^2-2Comega+7A^2.                                   (12)
```

Substitution of `(10)--(11)` in `(12)` gives the factored `D` in `(4)`.
Now `k=omega/D`, `h=A/D`, and `m=3theta+14A` give the remaining formulas.

Conversely, direct substitution of `(4)` verifies all three original cubic
multiplication laws, `(9)`, `(12)`, and

```text
Ck-mh=1,
D(7h^2+3k^2)=1+2Ck,
hD(9h+14k)=km+3hC^2.                                      (13)
```

It also gives `h/k=z`.  Hence every generator of `K(U)` belongs to `K(z,C)`
and the reverse inclusion is immediate from `(2)`.  This proves the
birational assertion without mistaking it for an affine isomorphism.

## 2. The cubic spectral split is a resolved chart pole

The five slopes in THM-3827 divide into roots of `q` and roots of `r`.  The
elementary identity

```text
b+q=2r                                                       (14)
```

shows why the three cubic slopes are different.  If `r(alpha)=0`, then
`alpha q(alpha) b(alpha)!=0`, while the numerator of `A` on `z=alpha` is

```text
alpha C(1+alpha^2 C).                                      (15)
```

Thus the pole `r=0` has exactly the two cancellation addresses

```text
C=0,                         C=-1/alpha^2.                  (16)
```

These are the two target-line values that appear as the two intrinsic
components of the cubic spectral fibre.  The two quadratic roots `q=0` do
not create this two-address cancellation.  This chart therefore explains,
rather than merely records, the three-versus-two split behind THM-3827's
disconnected-fibre passport.

Statement `(16)` is a boundary reading of the birational chart; regularity
at either address must be checked in the intrinsic surface chart.  No claim
that `(4)` itself is defined at `r=0` is made.

## 3. Exact counterexample-construction passport

A polynomial plane atlas through this chart must simultaneously arrange:

```text
(i)   A=zC(1+z^2C)/r is polynomial;
(ii)  k=r/(Cs) and h=zk are polynomial;
(iii) theta,D,m in (4) are polynomial;
(iv)  Jac(z,C)=lambda k r with lambda!=0;
(v)   one cubic spectral value reaches both cancellation branches (16);
(vi)  the required codimension-one coverage survives the chart boundary. (17)
```

This replaces a large coefficient search by a denominator-and-area-form
problem on a rational surface.  It also identifies what the rational chart
destroys: regularity at `C=0`, `s=0`, and `r=0` is not visible until the
corresponding numerator cancellation and intrinsic sidecar are restored.

No polynomial atlas and no planar Jacobian counterexample is constructed.
**QED, pending independent hostile audit.**
