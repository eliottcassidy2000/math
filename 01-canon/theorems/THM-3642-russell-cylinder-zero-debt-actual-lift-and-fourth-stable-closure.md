---
id: THM-3642
title: "Russell-cylinder zero-debt actual lift and fourth-stable closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the displayed
  degree-six zero-second-debt fold, a
  deterministic exact target-ring certificate reaches J_0=1 and
  J_1=J_2=0.  For that fold and a symmetric-curvature degree-seven control,
  universal endpoint identities force respectively
  lambda(J_4)=365888/6561 and 5440/81 whenever J_0=1 and J_2=0.  The
  fourth-order identities quantify over arbitrary target two-forms and do
  not require J_1 or J_3.  The fixed compiler and two fixed folds are the
  exact scope; no all-order, Keller-pair, or JC(2) claim is made.
source: root/zero-debt-next-order/2026-08-21
depends_on:
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
related:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness
  - THM-3639-russell-cylinder-all-retained-cells-universal-second-stable-debt
script: 04-computation/jc2_russell_cylinder_zero_debt_actual_lift_fourth_stable_closure_thm3642.py
output: 05-knowledge/results/jc2_russell_cylinder_zero_debt_actual_lift_fourth_stable_closure_thm3642.out
script_sha256: 348accd3bd488a2ca146d18000c74a7613b6801d2c3b1c22c634d5816968c681
output_sha256: 5d1e7e975335c60a926a0e45f282deb9bb9f37c4a1dbe0b53bb81290b2438d81
hash_basis: raw LF bytes
audit: >
  PASS -- independent hostile audit reconstructed the exact Q6 target
  representatives, derivative hashes, modular pivot skeleton, selected
  618-by-618 rational solve, and full J1=J2 residual identities.  A separate
  row-space reconstruction recovered every coefficient in both universal J4
  identities and proved the 105-element two-form jet universe complete,
  including the invisible degree-five boundary.  Normal, optimized, and
  stored transcripts are byte-identical at 912 active gates; AST0, hashes,
  docs, diff, and whitespace checks pass.
---

# THM-3642 -- Russell-cylinder zero-debt actual lift and fourth-stable closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The exact
certificate and both universal endpoint identities were independently
reconstructed on the frozen bytes.  Their restricted scope is part of the
theorem.

All rings and points are over `C`; the certificate itself is rational.

## 0. Setup and theorem

Use

```text
R=C[b,c,e]/(c^2 e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3).       (1)
```

Put `y=c/3`, `z=e+3`, and consider the quadratic fold

```text
(x,t) |-> (x,Q(x)+t^2,w=t).                            (2)
```

The two collision polynomials are

```text
Q_6=Q_1-(259/36)x^2(x^2-1)^2,

Q_1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,                (3)

Q_*= -x^7-(27/4)x^6+3x^5+18x^4-3x^3
       -(27/2)x^2+x-3/4.                               (4)
```

Both have retained values, slopes, and ordinary-triple tangents

```text
Q(-1,0,1)=(-3,-3/4,-3),
Q'(-1,0,1)=(-9/2,1,9/2),

(y',z')=((1,-9),(1,4),(1,9)).                          (5)
```

Their curvature triples are respectively

```text
Q_6''=(-451/18,-251/9,-163/18),
Q_*''=(-27/2,-27,-27/2),                               (6)
```

and both lie on the zero-second-debt wall

```text
5Q''(-1)+13Q''(1)+243=0.                               (7)
```

For a target pair `F^#,G^# in R[w]`, write its pulled-back source Jacobian
as

```text
Jac_x,t(F^#,G^#)=sum_(n>=0)t^n J_n(x).                  (8)
```

The theorem has two parts.

1. For `Q_6`, there is one deterministic actual target-ring choice through
   target stable order three with

   ```text
   J_0=1,                 J_1=0,                 J_2=0  (9)
   ```

   as exact polynomials in `Q[x]`.
2. For every target pair, indeed for every target two-form, the identities
   in Sections 2 and 3 hold.  Consequently `(J_0,J_2)=(1,0)` forces

   ```text
   Q_6: lambda(J_4)=365888/6561,
   Q_*: lambda(J_4)=5440/81,                            (10)

   lambda(P)=(5P(-1)-18P(0)+13P(1))/18.
   ```

Thus the exact `Q_6` lift reaches the second-stable equality wall but cannot
extend to a constant source Jacobian.  The conclusion is only about the two
fixed folds `(3)`--`(4)` and compiler `(1)`.

## 1. Deterministic actual `Q_6` lift

Take the zero-stable restrictions

```text
U=c,                         V=e.                       (11)
```

List raw target monomials `b^i c^j e^k` in nested `(i,j,k)` order, using the
restriction-degree weights

```text
(deg gamma(b),deg gamma(c),deg gamma(e))=(24,17,14).    (12)
```

These are **restriction-degree cutoffs**, not target total degrees, and no
minimality is claimed.  At cutoff `126` there are `105` monomials.  The
`143 x 210` rational system with columns

```text
[c' gamma(m)]_m,                  [-e' gamma(m)]_m      (13)
```

has rank `142`; setting its RREF free variables to zero selects actual
target representatives `F_1,G_1` with

```text
c' gamma(G_1)-gamma(F_1)e'=1.                           (14)
```

Their restriction data are

```text
                       degree   coefficient hash
gamma(F_1)               126    648b6d09...bbbd75
gamma(G_1)               123    e9363ac40...34587
delta(F_1)               120    8e44affbd...1eb22
delta(G_1)               117    2fe85085f...05dc6.      (15)
```

The actual coefficient vectors, not merely their restrictions, are frozen.
Their hashes bind every nonzero coefficient to its target monomial as
semicolon-separated `i,j,k:numerator/denominator` strings:

```text
F_1: support 71, hash 88c69cb9...e4da77,
G_1: support 68, hash a21b37f6...57267.                 (16)
```

For cutoff `240` there are `561` raw monomials.  The coupled `J_1,J_2`
operator has `623 x 2244` entries over `F_1000003` and modular operator and
augmented ranks both `618`.  This modular calculation is only a pivot
skeleton.  It selects column and row hashes

```text
columns 38236094...d91b1b,          rows 46fcca02...04492d. (17)
```

Rebuilding that selected `618 x 618` square over `Q`, solving it exactly,
and substituting into the **full rational polynomial identities** proves
`J_1=J_2=0`.  No claim is made that the full rational operator has rank
`618`.  The chosen target coefficient vectors are

```text
F_2 support 185, hash efa9ce67...0c953,
G_2 support 185, hash 4e69d81e...2b449,
F_3 support 185, hash 1a2ee587...3b756,
G_3 support  63, hash 6751af82...8a267.                 (18)
```

Their restriction degrees and hashes are

```text
gamma(F_2): 235, b9166bc8...e7d09,
gamma(G_2): 232, c85105a3...f347,
gamma(F_3): 240, f8ee4cfd...c495c,
gamma(G_3): 238, 64a56487...3fb11.                     (19)
```

The displayed hashes in `(15)`--`(19)` are abbreviated prose labels; the
companion prints and gates their full SHA-256 values.

## 2. Universal `Q_6` fourth-stable identity

For arbitrary target data, with primes denoting `x` derivatives,

```text
lambda(J_4)=
  (16246280/531441)J_0(-1)-(4489/6561)J_0'(-1)
  -(5/81)J_0''(-1)-(64/81)J_0'(0)
  +(13390648/531441)J_0(1)-(6559/6561)J_0'(1)
  -(13/81)J_0''(1)
  +(2012/2187)J_2(-1)+(5/27)J_2'(-1)
  -(2012/2187)J_2(1)-(13/27)J_2'(1).                  (20)
```

Setting `J_0=1,J_2=0` gives

```text
lambda(J_4)=365888/6561 !=0.                            (21)
```

No `J_1` or `J_3` hypothesis enters `(20)`.

## 3. Universal `Q_*` fourth-stable identity

For the degree-seven control,

```text
lambda(J_4)=
  (2300/81)J_0(-1)-(1/9)J_0'(-1)-(5/81)J_0''(-1)
  +(3140/81)J_0(1)-(7/9)J_0'(1)-(13/81)J_0''(1)
  +(4/9)J_2(-1)+(5/27)J_2'(-1)
  -(4/9)J_2(1)-(13/27)J_2'(1).                        (22)
```

Hence

```text
J_0=1, J_2=0             ==>       lambda(J_4)=5440/81 !=0. (23)
```

Again `J_1,J_3` are absent.

## 4. Why the identities cover actual pairs

In local target coordinates `(y,z,w)`, write an arbitrary two-form as

```text
omega=A dy^dz+B dy^dw+C dz^dw.                         (24)
```

On the folded source its Jacobian density is

```text
A(Y_x Z_t-Y_t Z_x)+B Y_x+C Z_x.                        (25)
```

Equations `(20)` and `(22)` are linear in `(A,B,C)` and their target jets.
The companion checks all `35` coefficient monomials of total target degree
at most four in each of the three slots, for `105` basis two-forms.  Higher
target degree has source order too large to enter `J_4` values, `J_2` first
derivatives, or `J_0` second derivatives.  No closedness, decomposability,
rank assumption, or division is injected.  Therefore the identities hold
for arbitrary two-forms and in particular for `dF^# wedge dG^#` from every
actual target-ring pair.

## 5. Actual local controls and the stopping boundary

The lower retained jets in `(20)`--`(22)` are nonvacuous.  Put `G=w` and

```text
F=F_0(y,z)+w^2F_2(y,z).                                 (26)
```

The exact finite polynomials `F_0,F_2` are printed in the companion.  For
each fold they give, at every retained branch,

```text
(J_0,J_0',J_0'')=(1,0,0),        (J_2,J_2')=(0,0),
J_1=J_3=0.                                             (27)
```

Their fourth-stable value vectors are

```text
Q_6: (97158947953,-55841390863,-37631662223)/997928100,

Q_*: (7993879/91260,-29539213/456300,-13841293/456300), (28)
```

whose `lambda` values are exactly `(21)` and `(23)`.  These are actual
finite local target-ring controls.  They are not additional global
`J_0=1,J_2=0` polynomial pairs; the global `Q_6` nonvacuity comes from
Section 1.

## 6. Scope and non-consequences

This result fixes compiler `(1)`, the quadratic fold `(2)`, and
the two polynomials `(3)`--`(4)`.  It does not cover another collision
polynomial, a nonquadratic fold, or arbitrary maps of `A^2`.  It neither
constructs nor excludes every Keller pair and gives no proof of `JC(2)`.
The two-dimensional Jacobian conjecture remains **OPEN**.

The exact transcript reproduces under both normal Python and `python -O`,
matches byte-for-byte, and reports zero assertion AST nodes.  The independent
audit reconstructed both the deterministic actual lift and the universal
fourth-stable identities before promotion.
