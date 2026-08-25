---
id: THM-3688
title: "Russell-cylinder exceptional-quartic actual-ring J1/J2 lift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over the
  irreducible quartic field of THM-3683, the uniform actual-ring J0=1 lift of
  THM-3687 extends through the full polynomial identities J1=J2=0.  The
  THM-3703 Apéry module compiles target membership to an 89-gap quotient: a
  177-square over the quartic field replaces the former 395-row raw solve.
  The 174 transverse equations are the high-degree part of one first-order
  operator U |-> U'+AU and reduce to exact polynomial division by A.  All
  four complex embeddings are covered.  No J3, J4, all-order, global-pair,
  Keller-map, or JC(2) conclusion is claimed.
source: jc_zero_debt_lift / exceptional-quartic gap-compiled coupled lift, 2026-08-22
audit: >
  INDEPENDENT HOSTILE AUDIT PASSED (jc_radial_escape_probe, 2026-08-23).
  The audit rederived irreducibility, squarefreeness, the exact discriminant
  and all four root embeddings; the J0/J1/J2 coefficient laws; every licensed
  division; actual-target-ring typing through the THM-3703 basis; and the
  precise J3/globalization boundary.  It repaired only quotient-space,
  matrix-solve, dependency-slug, and stopping-scope wording.  Independent
  normal and optimized executions line-normalize to the frozen output and
  the recorded script/output hashes.
depends_on:
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
  - THM-3687-russell-cylinder-exceptional-quartic-actual-j0-lift
  - THM-3703-russell-cylinder-exceptional-quartic-sagbi-module
related:
  - THM-3680-russell-cylinder-qdagger-coupled-stable-lift
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
  - THM-3737-russell-cylinder-exceptional-quartic-jacobian-image-hyperplane
  - THM-4039-exceptional-quartic-j3-lift-and-adjacent-gate-rigidity
  - THM-4043-exceptional-quartic-shifted-stable-identities-and-j6-lift
  - THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction
script: 04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j1_j2_lift_thm3688.py
output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_exact_j1_j2_lift_thm3688.out
script_sha256: 02cd67446b18b3863bc3665d48a6c5cccda81c394f94b754d2b90b1597c53ba6
output_sha256: 0c6bce1baef8dab98cd63219538c5e14aa9eae05903dc85fee647f776b514459
hash_basis: raw LF bytes
---

# THM-3688 -- all four exceptional folds pass the actual `J_1,J_2` gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The four
sixth-debt folds of THM-3683 are not stopped by either of the next two full
stable equations.  One exact certificate over their common quartic field
constructs actual target-ring coefficients through order three.

All rings below have characteristic zero.

## 1. Frozen field, ring, and coupled equations

Retain THM-3683/3687's field and exceptional fold

```text
F_6(T)=72783360T^4-77822208T^3-28419741T^2
                         +7849770T-1276420,
K=Q[alpha]/(F_6(alpha)),
q=Q_alpha(x).                                             (1)
```

Use the Russell-cylinder target ring and its restriction algebra

```text
R_K=K[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
B=(D-1)(D+2)^2,       C=xD(D+2),       E=q(D+3),
S=K[B,C,E] subset K[x].                                  (2)
```

Write `gamma` for restriction to `q=Q_alpha(x)` and `delta` for one normal
`q` derivative followed by restriction.  THM-3687 supplies actual
`F_1,G_1 in R_K` with

```text
C'G_1-F_1E'=1.                                           (3)
```

As in the independently audited expansion of THM-3680, introduce actual
target coefficients `F_2,G_2,F_3,G_3`.  Their full source-polynomial stable
equations are

```text
J_1=K_1+2L_0(F_2,G_2),

J_2=K_2+M(F_2,G_2)+3L_0(F_3,G_3),                       (4)

L_0(P,Q)=-E'P+C'Q,
M(P,Q)=P'G_1-2PG_1'+2F_1'Q-F_1Q',                       (5)

K_1=2C'delta(E)+F_1'G_1-F_1G_1'-2delta(C)E',

K_2=3C'delta(G_1)+2F_1'delta(E)+delta(C)'G_1
       -F_1delta(E)'-2delta(C)G_1'-3delta(F_1)E'.        (6)
```

The exact degrees are

```text
deg K_1=392,                    deg K_2=207.             (7)
```

No retained-point or truncated-jet equality substitutes for `(4)`: the
companion checks the complete polynomials in `K[x]`.

## 2. The 89-gap quotient compiler

THM-3703 gives a monic free presentation

```text
S=direct_sum_(r=0)^17 K[Z]p_r,        deg Z=18,          (8)
```

with one monic basis element in every degree of

```text
Gamma=<18,21,30,71,83,124>.
```

There are exactly `89` missing degrees, the largest is `169`, and every
degree at least `170` occurs.  Put

```text
S_(<=375)=S intersect K[x]_(<=375).
```

Descending against the monic basis defines the exact `K`-linear quotient
remainder

```text
rho:K[x]_(<=375) -> K[x]_(<=375)/S_(<=375),
                                      dim_K im(rho)=89. (9)
```

The Bezout identity `(3)` turns every polynomial solution of `L_0(F,G)=h`
into

```text
F=F_1h+C'U,                    G=G_1h+E'U.              (10)
```

Conversely every solution has this form.  Requiring the pair in `(10)` to
belong to `S^2` is therefore the gap-membership problem

```text
Phi(U)=(rho(C'U),rho(E'U)).                              (11)
```

For `deg U<=355`, `(11)` has `356` field variables and `178=2*89` quotient
coordinates.  The first `177` monomial variables and first `177` quotient
coordinates form an exact invertible square over `K`; its expansion in the
power basis `(1,alpha,alpha^2,alpha^3)` is only

```text
177 by 177 over K,                  708 by 708 over Q.   (12)
```

The sole omitted coordinate is the `G`-copy of the last gap, degree `169`.
The companion does not infer `(12)` from modular failure: it rebuilds the
square over `K`, solves it exactly over `Q`, reconstructs field elements, and
checks the final omitted coordinate and all target-membership remainders.

For `h=-K_1/2`, this produces one base pair in `S^2`.  To parametrize the
transverse correction, use free multiplier coordinates

```text
deg U=182,183,...,355,                 174 directions.  (13)
```

The five immediately preceding free coordinates `177,...,181` are omitted
from the transverse square.  This is the structural coordinate explanation
of the old split-prime selector's `179=5+174` count.  The proof does not need
to assert that each provisional coordinate direction is separately an actual
target pair: after solving the transverse system it gates the one final
linear combination in all `178` membership coordinates.

## 3. The transverse operator is `D+A`

On a homogeneous pair `(C'U,E'U)`, direct expansion of `(5)` and `(3)` gives

```text
M(C'U,E'U)=U'+AU,                                         (14)

A=C''G_1-2C'G_1'+2F_1'E'-F_1E''.                         (15)
```

Indeed the coefficient of `U'` is precisely `C'G_1-F_1E'=1`.  This is the
mechanism behind the former numerical response rank.  Exact leading terms
give

```text
deg A=214,
lc(A)=162 lc(E)lc(F_1) !=0.                              (16)
```

Thus the directions `(13)` meet the high `J_2` rows `396,...,569`
triangularly, with the same nonzero diagonal `(16)`.  More strongly, no
separate `174 by 174` transverse matrix solve is needed.  If
`(F_2^0,G_2^0)` is the base `J_1` solution, divide

```text
-(K_2+M(F_2^0,G_2^0)) by A.                             (17)
```

The quotient coefficients in degrees `182,...,355` are exactly the `174`
parameters.  The nominal rows `570,571,572` say exactly that `(17)` has no
quotient terms above degree `355`; the exact division passes this gate.  The
resulting actual pair satisfies

```text
K_2+M(F_2,G_2) in K[x]_(<=395).                         (18)
```

Applying the same gap compiler to `h=-(K_2+M(F_2,G_2))/3` produces actual
`F_3,G_3 in S`.  Its omitted gap-`169` coordinate is exactly zero.  Full
substitution then proves

```text
J_0=1,                         J_1=0,       J_2=0        (19)
```

as identities in `K[x]`.

## 4. Actual target-ring certificate

Every canonical basis element `Z^kp_r` in `(8)` was constructed in THM-3703
using target-ring addition, multiplication, constant shifts, and units of
`K`.  Coefficients indexed by that basis therefore type actual elements of
`R_K`; they are not arbitrary source-polynomial solutions.

The frozen certificate is

```text
       target support  restriction degree  canonical-target hash
F_2         287               375  5873ad24abe36b3b88c91f683dbf2915714fcbf9ac4586a7c5830a8d2db95058
G_2         284               372  1d793954556762bf9b456c9167f895c7d2c214ec6491c372b5cf43ae3b8275f3
F_3         108               196  4f3cfaf8f6ed40dad74d542d437dc6d429ddf860cdb154922ffb4914d75116c8
G_3         287               375  e2e0059c3e55f86f5b688e982e7a35958e5df12f37c98cb7efc2e6ded237cb84. (20)
```

The companion and transcript contain the full hashes.  Each field coefficient
is serialized in the ordered power basis `(1,alpha,alpha^2,alpha^3)`.

Because `F_6` is irreducible, every one of its four complex roots defines an
embedding `K -> C`.  Base-changing `(19)` along any embedding proves the same
actual `J_0,J_1,J_2` lift for all four exceptional folds.

## 5. Strict stopping boundary

This is a substantial continuation result, not a planar Jacobian
counterexample.  It does not prove

- `J_3=0`, `J_4=0`, or any later stable equation;
- that the finite stable series algebraizes to global target elements;
- a positive global pair on the Russell cylinder;
- a polynomial Keller map or its noninjectivity; or
- a counterexample to `JC(2)`.

THM-3629 becomes applicable only after a positive **global** pair is
constructed.  Downstream THM-3737 reduces the `J_3` stage to the scalar gate
`Lambda(D_3)=0`, and THM-4039 proves that gate and obtains stagewise actual
coefficients through `J_4`; THM-4043 then uses shifted stable identities to
reach `J_6`.  Those later theorems do not retroactively add a `J_3` claim to
the frozen certificate proved here.  THM-4046 reaches `J_7` and proves the
first unavoidable scalar debt at `J_8`.

## 6. Reproduction

From the repository root:

```bash
python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j1_j2_lift_thm3688.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j1_j2_lift_thm3688.py
```

Normal and optimized exact runs must byte-match the stored LF transcript.
The script contains no Python `assert` statements.  **QED.**
