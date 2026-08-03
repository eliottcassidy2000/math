---
id: THM-3327
title: "Positive-moment-tensor covering-arc and factorial phase dispersion"
status: >
  PROVED + INDEPENDENTLY AUDITED.  For a fixed order m>=2, if every entry of
  the symmetric moment tensor on a chosen basis is a strictly positive real,
  then a null m-th moment forces the nonzero coefficient phases to have
  shortest closed covering-arc width strictly greater than pi/m.  The bound
  is sharp in the positive-tensor class; strict positivity and m>=2 are
  load-bearing.  In the standard monomial basis of the factorial functional,
  the hypothesis holds at every order.  More strongly, L(P)=L(P^2)=0 forces
  covering width greater than pi and at least three noncollinear coefficient
  rays.  On a primitive cyclic eigenspace of order q>=3 that standard-basis
  gate is automatic: every occupied monomial orbit already contains q equally
  spaced coefficient phases.  The result is basis-dependent and proves
  neither FC(n) nor HFC(n).
audit: >
  An independent proof audit checked the closed-arc endpoint case, the pure
  and mixed multinomial terms, the m=1 opposite-ray boundary, a nonnegative-
  tensor equality witness, strict-positive sharpness by simple-root
  perturbation, the factorial real-square argument, and the basis-change
  warning.  No computation is used.
source: root/creative-passports-ii/2026-08-03
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
related:
  - THM-3022-two-slot-moment-dichotomy-for-weight-sequences
  - THM-3301-symmetry-vanishing-is-mathieu-compatible
  - THM-3310-degree-four-cyclic-eigenspace-on-the-triangle
  - THM-3321-hesse-moment-kernel-and-cyclic-quartic-support-four-exclusion
---

# THM-3327 -- positive-moment-tensor covering-arc and factorial phase dispersion

**PROVED + INDEPENDENTLY AUDITED.**

This theorem extracts THM-3310's cubic phase argument from its Fourier chart
and gives the stronger first-two-moment boundary in the standard factorial
monomial basis.  It is a coamoeba prefilter, not a factorial-conjecture proof.

## 1. Covering width

For a finite nonempty set of phases `Theta` on the circle, define

```text
w(Theta)=min{ell in [0,2*pi]:
             Theta lies in a closed circular arc of length ell}.      (1)
```

This is the **shortest closed covering-arc width**.  It is not the usual
pairwise circular distance, which is always at most `pi`.  Zero coefficients
are omitted from the phase set.

Let `A` be a commutative complex algebra, let `Lambda:A->C` be complex-linear,
fix `m>=2`, and choose `u_1,...,u_r in A`.  Assume that for every weak
composition of `m`,

```text
Lambda(product_i u_i^k_i) is a strictly positive real,
k_i>=0,                     sum_i k_i=m.                  (2)
```

For `P=sum_i c_i u_i`, put `Theta(P)={arg(c_i):c_i!=0}`.

**Positive-tensor sector theorem.**  If `P!=0` and

```text
Lambda(P^m)=0,                                               (3)
```

then

```text
w(Theta(P))>pi/m.                                           (4)
```

## 2. Half-plane proof

Suppose instead that the width is at most `pi/m`.  A common rotation puts all
occupied coefficient phases in `[0,pi/m]`; it only multiplies `(3)` by a
nonzero scalar.  The multinomial expansion is

```text
Lambda(P^m)=
 sum_(k_1+...+k_r=m) multinomial(m;k)
   c_1^k_1 ... c_r^k_r Lambda(u_1^k_1 ... u_r^k_r).        (5)
```

Every term in `(5)` lies in the closed upper half-plane.  A zero sum of such
terms forces every term onto its boundary.  The pure terms `c_i^m
Lambda(u_i^m)` force each occupied coefficient phase to be `0` or `pi/m`.

If both endpoints occur, a mixed term `c_i^(m-1)c_j` has phase `pi/m` or
`pi-pi/m`, strictly between `0` and `pi`; its tensor entry is positive by
`(2)`.  If only one endpoint occurs, every term in `(5)` lies on one nonzero
real ray.  Both cases contradict `(3)`, proving `(4)`.

## 3. Sharp boundaries

The order condition is exact.  At `m=1` there is no mixed term.  Positivity of
all `Lambda(u_i)` only proves

```text
Lambda(P)=0  =>  w(Theta(P))>=pi,                          (6)
```

and two opposite rays attain equality.

Strict positivity of the mixed entries is also essential.  On `C[x,y]`,
define a linear functional at degree `m` by

```text
Lambda(x^m)=Lambda(y^m)=1,
Lambda(x^(m-k)y^k)=0,              1<=k<=m-1.              (7)
```

Extend it arbitrarily on monomials of other degrees.
For `zeta=exp(i*pi/m)`, the polynomial `P=x+zeta*y` has
`Lambda(P^m)=1+zeta^m=0` and width exactly `pi/m`.  Thus nonnegative tensors
permit equality.

The constant remains sharp even with strict positivity.  Replace the mixed
zeros in `(7)` by `epsilon>0`.  The corresponding two-variable moment
polynomial is

```text
F_epsilon(z)=1+z^m
 +epsilon*sum_(k=1)^(m-1) binomial(m,k)z^k.                (8)
```

The root `zeta` of `F_0` is simple, so for small positive `epsilon` it
continues to a root `z_epsilon -> zeta`.  Every degree-`m` tensor entry is now
strictly positive, while the phase width of `x+z_epsilon y` tends to `pi/m`.
No uniformly larger constant follows from `(2)` alone.

## 4. Factorial specialization

For the factorial functional in `n` variables,

```text
L(x^alpha)=alpha!,                                         (9)
```

take any finite support of distinct standard monomials `u_i=x^alpha_i`.
Every product in `(2)` is another monomial and

```text
L(product_i u_i^k_i)=(sum_i k_i alpha_i)!>0.              (10)
```

Consequently, for every `m>=2`,

```text
L(P^m)=0  =>  w(Theta(P))>pi/m                             (11)
```

in standard monomial coordinates.  THM-3018 supplies the equivalent product-
exponential integral and, on homogeneous forms, the simplex-moment model.

The first two factorial moments give a stronger joint boundary.  Suppose the
phases lie in a closed semicircle and rotate them into `[0,pi]`.  Since

```text
L(P)=sum_alpha alpha!*|c_alpha|*exp(i arg(c_alpha)),       (12)
```

has positive weights, `L(P)=0` forces every occupied phase to an endpoint.
Thus `P=exp(i*phi)R` for a nonzero real polynomial `R`.  But

```text
L(R^2)=int_([0,infinity)^n) R(x)^2 exp(-sum_i x_i) dx>0.   (13)
```

Therefore

```text
L(P)=L(P^2)=0  =>  w(Theta(P))>pi.                        (14)
```

The same proof for a homogeneous form uses the strictly positive simplex
integral of a nonzero real polynomial square.  Hence `(14)` applies to both
full `FC(n)` and `HFC(n)` in their standard monomial coordinates.

Equation `(12)` writes the origin as a strict convex combination of the
occupied unit phase vectors.  The strict inequality `(14)` says those rays
are not contained in any closed half-plane through the origin.  The origin is
therefore in the two-dimensional interior of their convex hull, and at least
three noncollinear coefficient rays are necessary.

## 5. Basis boundary and cyclic-quartic instance

Positivity in `(2)` belongs to the labelled basis.  It does not survive an
arbitrary change of coordinates.  In particular, THM-3310 uses the Fourier
basis

```text
zbar, z^2, z zbar^2, z^3 zbar, zbar^4.                   (15)
```

There is a sharp symmetry-vacuity boundary before changing basis.  Suppose an
order-`q` automorphism `tau` permutes a real monomial basis, preserves the
functional, and

```text
tau(P)=chi*P,                  chi primitive of order q.   (16)
```

Functional invariance gives the usual selection rule

```text
L(P^j)=L(tau(P^j))=chi^j*L(P^j).                          (16a)
```

For `q>=3`, the moments `j=1,2` therefore vanish automatically.
Every occupied monomial orbit has length `q`: an orbit of shorter length
`ell` would force its coefficient to equal `chi^(-ell)` times itself.  Along a
full orbit the coefficients differ successively by `chi^(-1)`, so its phase
set is a rotated copy of all `q`-th roots of unity.  Its covering width is

```text
2*pi*(q-1)/q.                                             (17)
```

For `q>=3`, `(17)` is already greater than `pi`.  Thus the standard-basis
first-two-moment gate `(14)` is automatic on a primitive cyclic eigenspace; it
does not by itself prune that symmetry cell.  At `q=2`, every occupied orbit
is antipodal.  If the total phase set has width `pi`, all antipodal orbit axes
coincide, so `P` is global-phase real and `(13)` gives `L(P^2)!=0`; if distinct
axes occur, its width is already greater than `pi`.

The cyclic quartic has `q=3`, so each occupied barycentric monomial orbit
already has width `4*pi/3`.  The standard-basis theorem cannot be transported
through the orbit Fourier transform as an additional exclusion.  What
THM-3310 verifies anew in the Fourier basis is `(2)` at `m=3`: all `35` cubic
tensor entries are positive.  Its strict `pi/3` barrier is exactly `(4)`;
its lopsided-axis boxes use additional quantitative tensor data.

The theorem is only a necessary phase gate.  It neither produces a zero
outside the forbidden arc nor makes different bases physically identical.
`FC(n)`, `HFC(n)`, and the full-support-five cyclic quartic remain open.

QED.
