---
id: THM-2723
title: "Split exact-square-prefix rational-primitive pole capacity"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  For every
  polynomial split exact-square prefix and every constant-coefficient reduced
  Faber mate of arbitrary degree and parity, the three-observable Hamiltonian
  identity gives Phi_Q'=Psi_Q'=0 and U R_Q'=kappa.  The rational-primitive
  lemma forces U to be constant or a power of one linear factor.  Therefore
  q=A_src/U has at most one finite pole and at most two pole points on the
  source P1.  This is a source-side capacity theorem, not an emptiness result;
  q=0 and all odd Faber seeds remain allowed.
source: root/split-rational-primitive-pole-capacity-2026-07-28
audit: coordinate-first-audit-2026-07-28 (independent derivation; final-text line audit; hostile q=0 and parity-boundary checks)
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
related:
  - THM-2071-quadratic-fiber-square-parity-gate
  - THM-2202-uniform-all-degree-quartic-pole-closure
  - THM-2713-split-prime23-component-divisor-budget-and-perfect-power-normal-form
  - THM-2726-a21-transverse-integral-split-response-three-pole-closure
---

# THM-2723 -- one rational primitive permits at most two `q`-pole points

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**

The third Faber observable is usually retained as a continuation sidecar.
On a polynomial split exact-square prefix it already imposes a global
restriction on the source line: its rational primitive has only one possible
finite denominator place.  The argument is independent of the mate's degree,
parity, and first two flux values.

## 1. Statement

Work over `C`.  Put `R=C[x]` and `K=C(x)`.  Let

```text
P=H^2+L in R[z],
H=V z^2+Bz+C_0,                L=A z+E,              (1)
```

where `V!=0` and `V` is a square in `K`.  Let `Q in R[z]` satisfy

```text
J_(x,z)(P,Q)=kappa in C*.                             (2)
```

Choose the polynomial square root `U in R\{0}` with `U^2=V` and set

```text
w=Uz+B/(2U).                                         (3)
```

Then

```text
P=w^4+p w^2+q w+r,                 p,q,r in K,
q=A/U.                                                (4)
```

Assume that in `K[w]` the target has a finite constant-coefficient reduced
Faber normal form

```text
Q=J(P)+sum_(m in M) c_m E_m,
J in C[T],                 c_m in C,
M subset {m>=1: 4 does not divide m},                 (5)
```

where `E_m=Pol_w(P^(m/4))`.  No parity restriction is imposed on `M`.
For the sums of the THM-2129 observables attached to `(5)`, one has

```text
Phi_Q'=0,                 Psi_Q'=0,
U R_Q'=kappa.                                         (6)
```

Moreover exactly one of the following holds:

```text
U in C*,                                                (7a)

U=u_0(x-a)^m,        u_0 in C*, a in C, m>=2.          (7b)
```

Consequently the rational function `q=A/U` has at most one finite pole and
at most two distinct pole points on `P1_x`.  Cancellation between `A` and
`U` can only lower this number.  If `q!=0`, `T=q^2` has the same pole support;
if `A=0`, both are zero and have no poles.

## 2. The split root is polynomial

Initially the split hypothesis gives a square root in `K`.  Write it as
`a/b` in lowest terms.  From

```text
a^2=b^2 V,                     V in C[x],             (8)
```

every irreducible factor of `b` would divide `a`, a contradiction.  Thus
`b` is a unit and, after absorbing a constant square, the root may be chosen
as `U in C[x]`.

Completing the square in `(1)` gives

```text
H=w^2+C_0-B^2/(4U^2),
L=(A/U)w+E-AB/(2U^2),                                (9)
```

which proves `(4)` and in particular the exact physical identity `q=A/U`.
All depressed coefficients lie in `K`; no polynomiality of `w` is asserted
or needed.

## 3. Coefficient comparison spends the third observable

The change `(x,z)->(x,w)` has determinant `U`, so `(2)` becomes

```text
J_(x,w)(P,Q)=kappa/U.                                (10)
```

For each reduced Faber seed THM-2129 proves

```text
J_(x,w)(P,E_m)
 =(w^2+p/4)Phi_m'+w Psi_m'+R_m'.                     (11)
```

The target shear `J(P)` contributes zero.  Since every `c_m` in `(5)` is a
constant, summing `(11)` gives the polynomial identity

```text
(w^2+p/4)Phi_Q'+w Psi_Q'+R_Q'=kappa/U                (12)
```

in `K[w]`.  Comparison of the `w^2`, `w`, and constant coefficients proves
`(6)`.  This uses neither nonsplit deck parity nor a prescribed value of
`Phi_Q` or `Psi_Q`; it uses only that their derivatives occur in distinct
coefficients.  Odd seeds merely change the three finite sums.

The Laurent recurrence defining each Faber seed has coefficients in
`Q[p,q,r]`.  Therefore every `Phi_m,Psi_m,R_m`, and hence `R_Q`, belongs to
`K`.  Thus the last equation in `(6)` is an honest rational-primitive
equation on the source line.

## 4. Rational primitives have one finite denominator place

Apply the elementary rational-primitive lemma of THM-2214, Section 7.3, to

```text
U in C[x]\{0},              S=R_Q in C(x),
U S'=kappa in C*.                                      (13)
```

It says that either `U` is constant and `S` is affine-linear, or, after a
translation `X=x-a`,

```text
U=u_0 X^m,              S=s_0+s_1 X^(1-m),
m>=2,                   u_0 s_1!=0.                  (14)
```

For completeness, a rational derivative has no simple pole.  If a
nonconstant `U` of degree `D` had `h` distinct roots, a rational primitive
of `1/U` would have map degree `D-h`, while its fibre over its value at
infinity has multiplicity `D-1`.  Hence `D-1<=D-h`, so `h=1`; direct
integration then gives `(14)` and excludes a simple root.

This proves `(7)`.  Since the numerator `A` in `(4)` is polynomial, a finite
pole of `q` can occur only at the sole root in `(7b)`.  The only additional
possible pole point on the projective source is infinity.  This proves the
capacity bound.

## 5. Sharp boundary and scope

The hypothesis `A!=0` must not be inserted.  The pair

```text
P=z^4+x=(z^2)^2+x,                 Q=z=E_1           (15)
```

has `J(P,Q)=1`, `U=1`, `A=0`, `q=0`, and `R_1=x`.  It is the minimal
hostile control showing that `(6)` does not itself make the split chart
empty and that an odd seed may carry the third flux at `q=0`.

The theorem supplies only a source-side pole budget.  To exclude a physical
trajectory one still needs a geometrically integral target response curve,
a nonconstant map from the source, and at least three distinct target pole
points of the same function `q`.  It does not derive the normal form `(5)`
for an arbitrary Keller pair, treat a nonsplit or nonpolynomial exact prefix,
or prove `JC(2)` or `DC(2)`.

The load-bearing assumptions are characteristic zero, `kappa!=0`, a
polynomial split exact-square prefix, and constant coefficients in the
reduced Faber expansion.  The target shear `J(P)` is harmless; allowing odd
seeds is harmless; cancellation or `A=0` only decreases pole support.

An independent hostile audit rederived the polynomial split root, the
Jacobian factor `U`, all three coefficient comparisons, `K`-rationality of
the observables, and the rational-primitive classification.  It checked the
constant-`U`, cancellation, odd-seed, and `q=0` boundaries against `(15)` and
then certified this final text line by line.
