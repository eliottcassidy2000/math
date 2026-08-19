---
id: THM-3549
title: "Low-complexity correction no-go for the torus collision quotient"
status: >
  PROVED + VERIFIED-EXACT / HOSTILE-AUDITED.  In collision coordinates the
  THM-3543 quotient has Jacobian -2w^2, generic degree three, and one
  contracted line.  No one-sided polynomial correction can make its
  Jacobian a nonzero constant.  Arbitrary corrections of total degree at
  most two give the unit coefficient ideal; through degree three the
  constant is forced to zero, with or without the inherited collision.
  Corrections affine in w fail in every u-degree.  Any collision-preserving
  Keller correction must instead give sorted final w-degrees at least (4,5).
source: boxeph-2026-08-18-jacobian-dephasing
depends_on:
  - THM-3543-torus-quotient-ramification-square-no-go
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-2071-quadratic-fiber-square-parity-gate
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
related:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-3546-invariant-graph-keller-descent-criterion
script: 04-computation/torus_quotient_correction_no_go_thm3549.py
output: 05-knowledge/results/torus_quotient_correction_no_go_thm3549.out
script_sha256: 4b60d1630431ea85ece7f66719ebb9e138674e0eb8df6d0e5524d70de6d10b9b
output_sha256: bde9c5edd7934ed528e02614ff26175e2e70661ada1967158a5f3229c48d3921
hash_basis: LF-normalized bytes
---

# THM-3549 -- low-complexity correction no-go for the torus collision quotient

**PROVED + VERIFIED-EXACT / HOSTILE-AUDITED.**  This theorem tests the most
literal attempt to repair the ramified but noninjective polynomial plane map
of THM-3543.  It closes several large correction classes and leaves a sharp
mixed high-transverse-degree cell.  It does not exclude arbitrary polynomial
corrections or construct a planar Keller counterexample.

All structural statements hold over a characteristic-zero field.  The exact
coefficient calculations were performed over `QQ`, hence remain valid over
`C`.

## 1. The quotient seed in collision coordinates

Put `w=2-3v-t` in the source quotient of THM-3543 and reverse the target
order.  The map is

```text
R=w(4v+6-3(v+1)^2 w),
S=w^2((v+1)(v+2)-(v+1)^3 w),                           (1)
```

with

```text
Jac_(v,w)(R,S)=-2w^2.                                  (2)
```

The inherited points

```text
p=(0,2),                 q=(-3/2,0)                    (3)
```

both map to `(0,0)`.  More precisely,

```text
(R,S)^(-1)(0,0)=V(w) union {p}                         (4)
```

set-theoretically, so `q` lies on a whole contracted line.

For target coordinates `(rho,sigma)`, elimination gives the essential
fibre cubic

```text
-2w^3+(4-3rho)w^2+rho^3-rho^2-18rho sigma
       +27sigma^2+16sigma=0.                           (5)
```

At `(rho,sigma)=(1,1)`, an exact lexicographic basis is

```text
650v-34w^2-133w+675,
(2w-5)(w^2+2w+5),                                     (6)
```

whose cubic discriminant is `-67600`.  Thus the generic degree is three;
the fibre `(4)` is exceptional.  Equations `(1)`--`(6)` are exact controls
for the correction calculation below.

## 2. One-sided corrections fail in every degree

Let

```text
P=R+A(v,w),                 Q=S+B(v,w).                (7)
```

The seed has

```text
grad R(q)=0,                grad S|_(w=0)=0.            (8)
```

Therefore keeping `R` (`A=0`) makes `Jac(P,Q)` vanish at `q`, for every
polynomial `B`.  Keeping `S` (`B=0`) makes the Jacobian vanish all along
`w=0`, for every polynomial `A`.  A nonzero constant Jacobian requires
simultaneous changes to both outputs.  No collision hypothesis is used.

## 3. Exact low-total-degree coefficient ideals

For `D=1,2,3`, let `A` and `B` range over all polynomials of total degree at
most `D` with arbitrary coefficients.  Constant terms are omitted in the
calculation because they affect neither the Jacobian nor the equal-value
collision equations.  Introduce an independent scalar `kappa` and take all
coefficients in

```text
Jac(R+A,S+B)-kappa.                                    (9)
```

Two exact universes were tested:

1. no value constraints; and
2. the collision-preserving equations

   ```text
   A(p)=A(q),                    B(p)=B(q).             (10)
   ```

For `D=1` and `D=2`, the Groebner ideal is the unit ideal in both universes.
For `D=3`, it contains `kappa` in both universes.  Consequently

```text
D<=2:  no constant-Jacobian correction exists at all;
D=3:   every constant Jacobian is zero.                (11)
```

This is `FINITE-EXACT`, not an extrapolation beyond degree three.  The
companion prints monomial counts, equation counts, basis sizes, and the
remainder verdict for all six cases under normal and optimized execution.

## 4. Corrections affine in the transverse variable fail in every degree

Set `u=v+1` and allow arbitrary polynomials `f,a,g,b in k[u]`:

```text
A=f(u)+a(u)w,                  B=g(u)+b(u)w.            (12)
```

For the corrected pair `(7)`, the coefficient of `w^3` in the Jacobian is

```text
-3u^2(ua'-a).                                          (13)
```

Constancy forces `ua'-a=0`, hence `a=cu` for a scalar `c`.  After this
substitution, the coefficient of `w^2` is

```text
cu-3u^3f'+6u^2b'-6ub-2.                               (14)
```

Its value at `u=0` is `-2`, a contradiction.  Therefore no correction
affine in `w`, regardless of its degree in `u`, has constant Jacobian.
Neither the collision condition nor a degree bound is used.

## 5. The surviving transverse-degree cell begins at `(4,5)`

Suppose `(7)` has nonzero constant Jacobian and preserves the collision
`(3)`.  It is noninjective, hence not a polynomial automorphism.  Apply
THM-2063, THM-2071, and THM-2118 to the source fibre coordinate `w` and to
the two target-pencil members `P,Q`.  They give

```text
deg_w P>=4,                       deg_w Q>=4.            (15)
```

Both outputs therefore require new terms of `w`-degree at least four, since
the seed degrees are two and three.

They cannot both have `w`-degree exactly four.  If their leading terms are
`A_4(u)w^4` and `B_4(u)w^4`, the coefficient of `w^7` in their Jacobian is

```text
4(A_4'B_4-A_4B_4').                                  (16)
```

Its vanishing makes the two nonzero polynomial coefficients constant-
proportional.  A nonzero target-pencil combination then cancels the
`w^4` term and has fibre degree at most three, contradicting the same three
canonical fibre gates.  After ordering,

```text
(deg_w P,deg_w Q)>=(4,5).                              (17)
```

Even the boundary case has a rigid top.  Write its transverse leading terms
as

```text
P=a(u)w^4+lower w-rows,       Q=b(u)w^5+lower w-rows.  (18)
```

The coefficient of `w^8` in the Jacobian is

```text
5a'b-4ab'.                                                (19)
```

Constancy forces `(a^5/b^4)'=0`.  Unique factorization and
`gcd(4,5)=1` then give

```text
a=c h(u)^4,                    b=d h(u)^5               (20)
```

for one polynomial `h` and nonzero constants `c,d`.  Thus the first open
box is already a common-power transverse skeleton, not a generic pair of
quartic/quintic coefficients.  Its lower rows must cancel the inherited
`-2w^2` and preserve the collision.  Equations `(18)`--`(20)` are a further
necessary gate, not an existence claim.

## 6. Search interpretation and failure boundary

The quotient already supplies polynomiality, a genuine collision, generic
degree three, and a simple square conductor.  The no-gos show why subtracting
that conductor is not a low-order perturbation problem.  A viable repair must
simultaneously:

1. alter both outputs;
2. introduce genuinely nonlinear transverse layers in both;
3. reach at least sorted `w`-degrees `(4,5)`;
4. retain the two-point collision without retaining the contracted divisor;
5. cancel every positive `w` coefficient of the Jacobian globally.

The theorem does not exclude corrections beyond total degree three with
nonlinear `w`-dependence, rational or algebraic changes of variables, or the
coordinate-hypersurface descent lane of THM-3546.  In particular it neither
proves JC(2) nor turns the THM-1300 map into a planar counterexample.
