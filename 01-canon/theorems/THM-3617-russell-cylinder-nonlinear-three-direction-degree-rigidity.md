---
id: THM-3617
title: "Russell-cylinder nonlinear three-direction degree rigidity"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  On every polynomial Russell-cylinder graph, the genuinely
  nonlinear target pair theta Y+R(B,C), V(B,C)S+Q(B,C) cannot have nonzero
  constant ordinary Jacobian when theta, [B^3]R, [B^2]V, and V(0,0) are
  nonzero and deg R,deg Q<=3, deg V<=2.  The x=0 trace is injective; three
  exact degree profiles violate Jung--van der Kulk divisibility; the sole
  cancellation h=-xq+n has [q]Jac=8 theta V(0,0).  This file is not a
  proved dependency before hostile audit and status promotion.
source: root / nonlinear_cylinder_shears follow-up, 2026-08-21
audit: >
  PENDING -- exact companion frozen provisionally; independent mathematical
  and replay audit required before promotion.
depends_on:
  - THM-3614-russell-cylinder-collision-free-full-linear-projection-rigidity
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3608-russell-cylinder-nonlinear-target-shear-rigidity
  - THM-3611-russell-cylinder-arm-separating-nonlinear-first-coordinate-rigidity
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
  - "Jung--van der Kulk: degree divisibility for polynomial automorphisms of the affine plane."
script: 04-computation/jc2_russell_cylinder_nonlinear_three_direction_degree_rigidity_thm3617.py
output: 05-knowledge/results/jc2_russell_cylinder_nonlinear_three_direction_degree_rigidity_thm3617.out
script_sha256: ee9d82cecb7e05408041b7443b8466d2f6e792c747c2daf20decc901267f6189
output_sha256: d5471aaed4655e64a84b381624eacbb6e8190a91b3483261d6d577319f417f06
hash_basis: raw LF bytes
---

# THM-3617 -- Russell-cylinder nonlinear three-direction degree rigidity

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**  The proof and exact controls are complete, but this file is
not a proved dependency until an independent hostile audit promotes it.

All rings and derivatives are over `C`.  Degrees of source polynomials mean
ordinary total degree in `(x,q)`; degrees of target polynomials mean ordinary
total degree in their displayed variables.

## 0. Statement

Retain the Russell graph functions

```text
D=1+x^2q,
c=xD(D+2),
b=(D-1)(D+2)^2,
e=q(D+3),                                               (1)

B=b+ch,
C=c,
Y=ce+(2b+4)h+ch^2,
S=((b+2)(e+3h^2)+ch(3e+h^2))/8                        (2)
```

for an arbitrary polynomial `h in C[x,q]`.  Let `T,U` be independent target
variables and take

```text
theta in C*,
R,Q in C[T,U],               deg R<=3, deg Q<=3,
V in C[T,U],                 deg V<=2,
lambda=[T^3]R != 0,
eta=[T^2]V != 0,             nu=V(0,0) != 0.           (3)
```

Define

```text
F=theta Y+R(B,C),
G=V(B,C)S+Q(B,C).                                      (4)
```

Then

```text
Jac_(x,q)(F,G) notin C*.                               (5)
```

The quantifiers in `(3)` are simultaneous.  Lower coefficients are arbitrary,
and no collision, symmetry, or genericity assumption is made on `h`.  The
coefficient `V(B,C)` of `S` is allowed to be genuinely nonlinear.

## 1. This family is beyond rank-two linear postcomposition

There is a useful closure of THM-3614.  If `(L,M)` is any rank-two linear
projection of `(B,C,Y,S)` and `H_1,H_2 in C[Z_1,Z_2]`, then the chain rule gives

```text
Jac(H_1(L,M),H_2(L,M))
 =Jac_(Z_1,Z_2)(H_1,H_2)(L,M) Jac(L,M).                 (6)
```

If the left side were a nonzero constant, both polynomial factors on the
right would be units of `C[x,q]`; in particular `Jac(L,M)` would be a nonzero
constant, contrary to THM-3614.  Thus every polynomial postcomposition through
a fixed rank-two linear target plane is already closed.

The family `(4)` does not secretly reduce to that corollary.  This must be
proved in the Russell quotient, not by treating all four ambient coordinates
as independent.  Write

```text
mathcal R=C[B,C,Y,S]/(CY-B(B+4)).                       (7)
```

Suppose, to the contrary, that `(F,G)=(H_1(L,M),H_2(L,M))` in `mathcal R`
for two linear forms `L,M`.  The functions `F,G` are algebraically
independent: on `C!=0`, the variables are `(B,C,S)`, `F` is nonconstant in
`(B,C)`, and `G` is affine in the free variable `S` with nonzero coefficient
`V`.  Hence `L,M` are algebraically independent and the substitution
`C[U,V]->mathcal R`, `(U,V)=(L,M)`, is injective.

Let `ell=L_S` and `m=M_S`.  If `ell=m=0`, every polynomial in `L,M` is
independent of `S`, contradicting the term `V(B,C)S` in `G`.  Otherwise,
differentiate `F=H_1(L,M)` with respect to `S` on `C!=0`.  It gives

```text
ell (H_1)_U + m (H_1)_V=0,
H_1=P(mU-ell V),
F=P(N),                 N=mL-ell M.                    (8)
```

Here `N` is an `S`-free linear form in `(B,C,Y)`.  In the normal form for
`C[B,C,Y]/(CY-B(B+4))`, a pure power `Y^j` is never reduced.  If the
`Y`-coefficient of `N` is zero, `P(N)` cannot contain `theta Y`.  If it is
nonzero and `deg P=j>=2`, `P(N)` has a nonzero pure `Y^j` term, again
impossible because `F=theta Y+R(B,C)` has `Y`-degree one.  Thus `deg P<=1`,
which cannot supply `lambda B^3`.  This proves that three target directions
are essential in the Russell quotient itself.

## 2. The source line is injective

Put `phi(q)=h(0,q)`.  On `x=0`,

```text
B=C=0,                 Y=4phi(q),
S=q+3phi(q)^2/4.                                      (9)
```

Therefore `(4)` restricts to

```text
q |-> (4theta phi(q)+R(0,0),
       nu(q+3phi(q)^2/4)+Q(0,0)).                     (10)
```

This map is injective: its first coordinate recovers `phi`, and its second
then recovers `q`, using `theta nu != 0`.  If `(4)` had nonzero constant
Jacobian, Gwozdziewicz's cited one-line theorem would make `(F,G)` a polynomial
automorphism of `A2`.  Jung--van der Kulk would then require the two nonconstant
ordinary degrees to be equal or the smaller to divide the larger.

## 3. Exact ordinary-degree profiles

Write

```text
c_7=x^5q^2,                  k=xq.                    (11)
```

### 3.1 Graph degree at least three

Let `d=deg h>=3`, with top homogeneous form `h_d`.  The Russell coordinates
have unique top forms

```text
B_(d+7)=c_7h_d,
Y_(2d+7)=c_7h_d^2,
S_(3d+7)=c_7h_d^3/8.                                  (12)
```

The unique top target monomials are `lambda B^3` in `F` and `eta B^2S` in
`G`; every other monomial allowed by `(3)` is strictly lower.  Hence

```text
F_(3d+21)=lambda c_7^3h_d^3,
G_(5d+21)=eta c_7^3h_d^5/8,
(deg F,deg G)=(3d+21,5d+21).                           (13)
```

Here `3d+21 < 5d+21 < 2(3d+21)`, so the smaller degree cannot divide the
larger.

### 3.2 Graph degree at most one and generic quadratic degree

If `deg h<=1`, then

```text
B_9=c_7k,                   S_13=c_7k^3/8,
F_27=lambda c_7^3k^3,       G_31=eta c_7^3k^5/8.       (14)
```

If `deg h=2`, write `h=h_2+h_1+n`.  The degree-nine and degree-thirteen
layers are

```text
B_9=c_7(k+h_2),             S_13=c_7(k+h_2)^3/8.       (15)
```

Thus, whenever `h_2!=-k`,

```text
F_27=lambda c_7^3(k+h_2)^3,
G_31=eta c_7^3(k+h_2)^5/8,
(deg F,deg G)=(27,31).                                  (16)
```

The degree pair `(27,31)` violates divisibility.

### 3.3 Cancelled quadratic degree

If `h_2=-k` but `h_1!=0`, the next layers are

```text
B_8=c_7h_1,                    S_10=c_7h_1^3/8,
F_24=lambda c_7^3h_1^3,       G_26=eta c_7^3h_1^5/8,
(deg F,deg G)=(24,26).                                  (17)
```

Again the smaller degree does not divide the larger.  Only

```text
h=-xq+n                                                     (18)
```

remains outside the degree argument.

## 4. The exceptional graph has an unavoidable trace coefficient

For `(18)`, direct differentiation of `(1)--(4)` at `x=0` gives

```text
F(0,q)=4theta n+R(0,0),
F_x(0,q)=theta(8q+3n^2)+3n R_T(0,0)+3R_U(0,0),
F_q(0,q)=0,
G_q(0,q)=nu.                                            (19)
```

Consequently

```text
Jac(F,G)(0,q)
 =nu(theta(8q+3n^2)+3n R_T(0,0)+3R_U(0,0)),
[q]Jac(F,G)=8theta nu != 0.                             (20)
```

This contradicts constancy and closes the last graph.  Together with
Sections 2--3, it proves the provisional claim `(5)` pending independent
hostile audit.

## 5. Preservation/loss ledger and sharp method boundaries

| item | exact content |
|---|---|
| source | every polynomial graph `w=h(x,q)` |
| target map | graph embedding into `(B,C,Y,S)`, followed by `(4)` |
| preserved | polynomiality, the labelled `x=0` trace, and ordinary source degree |
| first sidecar | one-line injectivity from `theta nu!=0` |
| second sidecar | Jung--van der Kulk degree divisibility |
| exceptional sidecar | `[q]Jac=8theta nu` |
| information not used | collisions, arm labels, formal puncture coordinates |
| still lost | higher target degrees, missing leading coefficients, arbitrary nonlinear pairs, and non-graph source planes |

Each nonvanishing hypothesis has an exact method-hostile boundary:

- if `theta=0`, the first line coordinate is constant, so the uniform line
  reconstruction used here disappears, as does `(20)`;
- if `nu=0`, the second line coordinate is constant, so the same uniform
  reconstruction and `(20)` again disappear;
- if `lambda=0`, the control `R=T^2`, `V=1+T^2`, `Q=0` has the divisible
  degree profile `(28,56)` at graph degree seven;
- if `eta=0`, the control `R=T^3`, `V=1`, `Q=T^3` has equal component
  degrees for every graph degree at least three;
- allowing `deg R=4` produces an equal-degree wall at graph degree seven
  for `R=T^4+T^3`, `V=1+T^2`, `Q=0`; and allowing `Q=T^6` produces the
  divisible profile `deg G=2 deg F` in the high-degree branch.

These are failures of this proof mechanism, not Keller pairs.  The theorem
does not claim the coefficient conditions or degree caps are necessary for
the no-go.  It also does not claim a result for implicit source planes or a
counterexample to JC(2).

## 6. Exact companion contract

The exact companion must verify without truth-bearing `assert` statements:

- the Russell identities and the line reconstruction `(9),(10)`;
- the quotient-safe one-linear-form obstruction `(7),(8)`, an ancillary
  ambient three-direction determinant, and the postcomposition chain rule
  `(6)`;
- the top forms `(12)--(17)` over explicit high-, low-, and quadratic graph
  universes, together with complete allowed-target monomial degree invoices;
- the divisibility inequalities on a stated integer grid;
- the exceptional formulas `(19),(20)` with arbitrary target coefficients;
- the `theta`, `nu`, `lambda`, and `eta` hostiles and the displayed
  higher-degree walls.

The companion checks the displayed algebra and finite support census.  The
all-degree inference uses the cited one-line and plane-automorphism theorems;
no finite degree box substitutes for those results.

Reproduction commands are

```text
python3 04-computation/jc2_russell_cylinder_nonlinear_three_direction_degree_rigidity_thm3617.py
python3 -O 04-computation/jc2_russell_cylinder_nonlinear_three_direction_degree_rigidity_thm3617.py
```

The normal and optimized transcripts must be byte-identical to each other and
to
`05-knowledge/results/jc2_russell_cylinder_nonlinear_three_direction_degree_rigidity_thm3617.out`.
