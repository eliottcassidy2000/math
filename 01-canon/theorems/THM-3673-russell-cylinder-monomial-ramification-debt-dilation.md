---
id: THM-3673
title: "Russell-cylinder monomial ramification debt dilation"
status: >
  PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.  For each of
  the two THM-3642 zero-second-debt collision polynomials Q_6 and Q_*, every
  monomial vertical fold q=Q(x)+alpha t^k, w=t with alpha nonzero and k>=2
  inherits the quadratic universal debt identity after the exact index
  dilation (J_0,J_2,J_4)->(J_0,J_k,J_2k).  The forced debts are respectively
  alpha^2*365888/6561 and alpha^2*5440/81.  Hence not even an arbitrary
  target two-form, and therefore no target pair, has nonzero constant source
  Jacobian on either fold.  The affine k=1, nonmonomial H, other Q, and
  nonordinary tangent-collision boundaries remain open.
source: kps-s194 / ramification row-space continuation, 2026-08-21
depends_on:
  - THM-3642-russell-cylinder-zero-debt-actual-lift-and-fourth-stable-closure
related:
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
script: 04-computation/jc2_russell_cylinder_monomial_ramification_debt_dilation_thm3673.py
output: 05-knowledge/results/jc2_russell_cylinder_monomial_ramification_debt_dilation_thm3673.out
script_sha256: 81d8237c9fefae07176d82bffb7d9b84763be3c13bb72755da2e0530d90db573
output_sha256: e07c9dd7986e3c4e45b6bff4bbf6947043b3fd110d804bf93931d777408170ae
hash_basis: raw LF bytes
---

# THM-3673 -- monomial ramification dilates the zero-debt obstruction

**PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**  The
quadratic fourth-stable obstruction of THM-3642 is not peculiar to the
quadratic exponent.  The stable character decomposition transports it to
every monomial ramification order.  This closes a genuine part of the
nonquadratic-fold frontier while retaining the affine boundary as a hostile
control.

All rings, germs, and two-forms are over `C`.

## 0. Setup and statement

Use the exponent-two Russell-cylinder compiler

```text
D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),
y=c/3,                z=e+3.                              (1)
```

The two fixed collision polynomials are

```text
Q_6=Q_1-(259/36)x^2(x^2-1)^2,

Q_1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,

Q_*=-x^7-(27/4)x^6+3x^5+18x^4-3x^3
       -(27/2)x^2+x-3/4.                                 (2)
```

Fix

```text
k>=2,                     alpha in C*,
q=Q(x)+alpha t^k,          w=t.                           (3)
```

Write an arbitrary target two-form in the completed regular coordinates at
the retained target point as

```text
Omega=A(y,z,w) dy^dz+B(y,z,w) dy^dw+C(y,z,w) dz^dw.      (4)
```

This is a strict enlargement of forms `dF^# wedge dG^#` from actual target
pairs.  Let its source pullback be

```text
E_(Q,k,alpha)^*Omega
     =j_(k,alpha)(x,t) dx^dt,

j_(k,alpha)=sum_(n>=0) t^n J_n(x).                       (5)
```

Put

```text
Lambda(P)=(5P(-1)-18P(0)+13P(1))/18.                    (6)
```

Then the identities in Sections 2 and 3 hold for every `(4)`.  In
particular, if `(5)` were the nonzero constant `kappa`, they would give

```text
Q_6:  0=alpha^2 kappa (365888/6561),
Q_*:  0=alpha^2 kappa (5440/81),                         (7)
```

which is impossible.  Thus neither fold in `(2)` admits even an arbitrary
target two-form with nonzero constant pullback, and a fortiori admits no
actual target pair with nonzero constant source Jacobian.

## 1. The exact character-decimation map

Expand the coefficient functions in the stable coordinate:

```text
A=sum_(r>=0) A_r(y,z)w^r,
B=sum_(r>=0) B_r(y,z)w^r,
C=sum_(r>=0) C_r(y,z)w^r.                               (8)
```

Let

```text
M=Jac_(x,q)(y,z),
Y_x=(partial_x+Q'(x)partial_q)y,
Z_x=(partial_x+Q'(x)partial_q)z,                        (9)
```

evaluated at `q=Q(x)+alpha t^k`.  Direct pullback gives

```text
j_(k,alpha)
 =k alpha t^(k-1) A(y,z,t)M
   +B(y,z,t)Y_x+C(y,z,t)Z_x.                            (10)
```

All occurrences of `y,z,M,Y_x,Z_x` outside the explicit powers of `t` are
series in `t^k`.  Therefore only the characters

```text
A_(k ell+1),             B_(k ell), C_(k ell)           (11)
```

can contribute to a coefficient `J_(kn)`.  Introduce a quadratic stable
variable `u` and define

```text
T_k(Omega)=A~ dy^dz+B~ dy^du+C~ dz^du,

A~=(k/2) sum_(ell>=0) A_(k ell+1)(y,z) u^(2ell+1),
B~=      sum_(ell>=0) B_(k ell)(y,z) u^(2ell),
C~=      sum_(ell>=0) C_(k ell)(y,z) u^(2ell).          (12)
```

Pull `(12)` back through the scaled quadratic fold

```text
q=Q(x)+alpha u^2,                w=u,                   (13)
```

and write its density as `sum_n u^n K_n(x)`.  The `A` prefactors agree
because

```text
2 alpha u * (k/2)u^(2ell+1)=k alpha u^(2ell+2),         (14)
```

while the `B,C` terms agree without a prefactor.  Hence, coefficient by
coefficient and with every `x` derivative retained,

```text
             J_(kn)^(d)(x)=K_(2n)^(d)(x)
             for all n,d>=0.                            (15)
```

This is the exact source-target map behind the dilation.  It is not a
numerical analogy and it loses only stable characters that cannot enter a
multiple-of-`k` coefficient.

The coefficient `alpha` is handled before decimation.  Rescaling a
quadratic source variable by `v=sqrt(alpha)u` transports the unit quadratic
identity to

```text
K_4 = alpha^2 (the J_0 block)
                 +alpha (the J_2 block).                (16)
```

Indeed the source-density coefficient of order `n` acquires the factor
`alpha^(-(n+1)/2)` under this change, so multiplying the unit identity by
`alpha^(5/2)` gives precisely `(16)`.  Combining `(15)` and `(16)` replaces
`(K_0,K_2,K_4)` by `(J_0,J_k,J_(2k))`.

## 2. Dilated identity for Q_6

For every two-form `(4)`,

```text
Lambda(J_(2k))
=alpha^2 [
   (16246280/531441)J_0(-1)-(4489/6561)J_0'(-1)
   -(5/81)J_0''(-1)-(64/81)J_0'(0)
   +(13390648/531441)J_0(1)-(6559/6561)J_0'(1)
   -(13/81)J_0''(1)
 ]
+alpha [
   (2012/2187)J_k(-1)+(5/27)J_k'(-1)
   -(2012/2187)J_k(1)-(13/27)J_k'(1)
 ].                                                        (17)
```

This is exactly THM-3642 equation (20), transported by `(15)--(16)`.
If the full source density is the constant `kappa`, then `J_0=kappa` and
all positive stable coefficients vanish.  Equation `(17)` becomes the
first contradiction in `(7)`.

## 3. Dilated identity for Q_*

The same transfer applied to THM-3642 equation (22) gives

```text
Lambda(J_(2k))
=alpha^2 [
   (2300/81)J_0(-1)-(1/9)J_0'(-1)-(5/81)J_0''(-1)
   +(3140/81)J_0(1)-(7/9)J_0'(1)-(13/81)J_0''(1)
 ]
+alpha [
   (4/9)J_k(-1)+(5/27)J_k'(-1)
   -(4/9)J_k(1)-(13/27)J_k'(1)
 ].                                                        (18)
```

For a nonzero constant source density, `(18)` becomes the second
contradiction in `(7)`.

## 4. Why k=1 is a genuine boundary

When `k=1`, the area term in `(10)` has no factor `t^(k-1)` forcing the
shifted `A_(kell+1)` character.  The decimated coefficient family then
contains the missing constant `A_0` direction, so `(12)` and the inherited
debt identity fail.  This is structural, not an artifact of the proof.

THM-3629 gives a global exact unit two-form and a local Darboux pair on this
affine boundary.  Independently, the exact companion enumerates the complete
degree-two two-form jet universe for `Q_6,k=1`: the lower retained rows have
rank `13`, while adjoining `Lambda(J_2)` raises the rank to `14`.  Thus no
analogue of `(17)` lies in that lower row space.

## 5. Scope and new frontier

The theorem closes

```text
Q in {Q_6,Q_*}, alpha!=0, monomial H(t)=alpha t^k, k>=2,
arbitrary target two-forms, hence arbitrary actual target pairs.          (19)
```

It does not cover

- a general polynomial `H` with two or more nonzero monomials;
- the affine or nonlinear `H'(0)!=0` mixed-pair frontier of THM-3629;
- another ordinary zero-debt polynomial from the THM-3641 hyperplanes;
- a THM-3641 nonordinary tangent collision;
- another compiler or an arbitrary polynomial map of `A^2`.

No Keller pair or counterexample is constructed.  `JC(2)` remains **OPEN**.
The sharp next probes are whether `(15)` has a triangular replacement for a
nonmonomial `H`, and whether the fourth-stable debt is a nonzero polynomial
on each ordinary zero-debt Hermite family rather than just at `Q_6,Q_*`.

## 6. Exact verification

The companion verifies the character exponent arithmetic symbolically and
checks `(17)--(18)` on the complete target two-form monomial jet universes
for `k=2,3,4,5`.  It includes the nonunit controls
`(Q_6,k,alpha)=(Q_6,3,2)` and `(Q_*,4,-3/5)`, and the affine rank-jump
control above.  The all-`k` quantifier is proved by `(8)--(16)`; the finite
enumerations are exact controls, not its replacement.

```bash
python3 04-computation/jc2_russell_cylinder_monomial_ramification_debt_dilation_thm3673.py
python3 -O 04-computation/jc2_russell_cylinder_monomial_ramification_debt_dilation_thm3673.py
```

Normal and optimized transcripts are byte-identical to the stored output.
The frozen companion reports `4202` active gates and contains no Python
assertion statements.
