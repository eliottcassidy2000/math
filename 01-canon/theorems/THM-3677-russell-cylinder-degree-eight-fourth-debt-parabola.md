---
id: THM-3677
title: "Russell-cylinder degree-eight fourth-debt parabola"
status: >
  PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.  The complete
  degree-at-most-eight principal Hermite family with the Q_6 retained values,
  slopes, and zero second debt is the affine plane Q_(p,r)=Q_6+pR_1+rR_2.
  Every member has a universal arbitrary-two-form fourth-stable identity;
  under J_0=1,J_2=0 its debt is
  64(729p-42120r^2+15192r+5717)/6561.  Thus every point off one rational
  parabola is closed, for every nonzero critical vertical displacement after
  character dilation and formal conjugacy.  On the parabola the retained
  fourth debt vanishes; this supplies no lift or Keller pair.  The displayed
  degree-eight Q_dagger is the sharp next candidate.  JC(2) remains open.
source: kps-s194 / ordinary zero-debt Hermite-family continuation, 2026-08-21
depends_on:
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
  - THM-3642-russell-cylinder-zero-debt-actual-lift-and-fourth-stable-closure
  - THM-3673-russell-cylinder-monomial-ramification-debt-dilation
  - THM-3675-russell-cylinder-critical-fold-formal-conjugacy-closure
script: 04-computation/jc2_russell_cylinder_degree8_fourth_debt_parabola_thm3677.py
output: 05-knowledge/results/jc2_russell_cylinder_degree8_fourth_debt_parabola_thm3677.out
script_sha256: 9d96fe3bac02121baa2a9a4746ba055893f4936023538ba339c6f82c42cb601d
output_sha256: 9c1e7280691a7af557eac10473ebe549282474ce8cbf72108acf9c1f8bf3834b
hash_basis: raw LF bytes
---

# THM-3677 -- the fourth debt cuts out one rational parabola

**PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**  The two
isolated fourth-debt calculations in THM-3642 extend to the complete
degree-at-most-eight ordinary zero-second-debt family.  The fourth debt is a
quadratic polynomial, not an unexplained pointwise constant.  Its zero set is
one explicit rational parabola.

All rings, jets, and two-forms are over `C`.

## 0. The complete Hermite plane

Put

```text
P=x^2(x^2-1)^2,
R_1=P(1-x^2),                  R_2=P(4-9x),

Q_(p,r)=Q_6+pR_1+rR_2,                                  (1)
```

where

```text
Q_6=Q_1-(259/36)P,
Q_1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4.                 (2)
```

Every member has

```text
Q(-1,0,1)=(-3,-3/4,-3),
Q'(-1,0,1)=(-9/2,1,9/2),
5Q''(-1)+13Q''(1)+243=0.                               (3)
```

Conversely, every polynomial of degree at most eight satisfying `(3)` is
exactly one `Q_(p,r)`.  Indeed the six value-and-slope conditions force the
difference from `Q_6` to be

```text
P(a+bx+cx^2).                                           (4)
```

The homogeneous second-debt row is `9a+4b+9c=0`, and `(R_1,R_2)` form a
basis of that two-dimensional kernel.  Equivalently, the companion verifies
rank seven for the seven confluent-Hermite rows on polynomials of degree at
most eight and verifies that `R_1,R_2` are an independent nullspace basis.

All these folds have the same ordinary retained tangent triple as `Q_6`.
In particular, this theorem does not enter a tangent-collision boundary.

## 1. Universal fourth-stable identity

Use the compiler and quadratic fold

```text
D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),
y=c/3,                z=e+3,

q=Q_(p,r)(x)+t^2,     w=t.                              (5)
```

For an arbitrary target two-form, write its source density as

```text
j(x,t)=sum_(n>=0)t^n J_n(x),
Lambda(S)=(5S(-1)-18S(0)+13S(1))/18.                  (6)
```

Define

```text
A_- =8(1819584pr+404352p-89092224r^2
        -13115349r+8006926)/3326427,

B_- =(177840r-26159)/28431,

A_0 =-8(1819584pr-2552472p+81746496r^2
         -74734101r-15181226)/3326427,

B_0 =64(9r-1)/81,
B_+ =-13(144r+65)/2187,

C_- =-8(2340r-503)/3159,
C_0 = 8(2340r-503)/3159.                               (7)
```

Then every arbitrary target two-form obeys

```text
Lambda(J_4)
 = A_- J_0(-1)+B_- J_0'(-1)-(5/81)J_0''(-1)
   +A_0 J_0(0)+B_0 J_0'(0)
   +B_+ J_0'(1)-(13/81)J_0''(1)

   +C_- J_2(-1)+(5/27)J_2'(-1)
   +C_0 J_2(0)-(13/27)J_2'(1).                         (8)
```

No closedness, decomposability, `J_1`, or `J_3` equation is used.  The
companion verifies `(8)` symbolically in the parameters `p,r` on the complete
`105`-monomial target two-form jet universe.  Hence it holds in particular
for `dF^# wedge dG^#` from every actual target pair.

## 2. The fourth-debt polynomial

If `J_0=1` and `J_2=0` as polynomial identities, `(8)` becomes

```text
Lambda(J_4)=D_4(p,r),

D_4=A_-+A_0
   =64(729p-42120r^2+15192r+5717)/6561.                (9)
```

This recovers both THM-3642 controls:

```text
Q_6=Q_(0,0):        D_4=365888/6561,

Q_*=Q_(0,1/9):      D_4=5440/81.                       (10)
```

The exact zero set is the rational parabola

```text
p=(520/9)r^2-(1688/81)r-5717/729.                      (11)
```

Thus every point off `(11)` has a nonzero arbitrary-two-form obstruction at
stable order four.  On `(11)`, equation `(8)` becomes a zero invoice; it does
not construct a two-form or solve any lower global identity.

## 3. Closure for every critical vertical displacement off the parabola

The character-decimation proof of THM-3673 is independent of the particular
collision polynomial.  It transports `(8)` from stable orders `(0,2,4)` to
`(0,k,2k)` under `q=Q+alpha t^k`, multiplying the `J_0` block by `alpha^2`
and the `J_k` block by `alpha`.

Moreover, the two `J_k` value coefficients in `(8)` sum to zero:

```text
C_-+C_0=0.                                              (12)
```

The formal-conjugacy argument of THM-3675 therefore applies verbatim.  For
any nonzero `H in t^2 C[t]`, a hypothetical constant source density becomes
an `x`-constant series after monomialization.  `Lambda` and all derivatives
kill that series, `(12)` kills the middle block, and the remaining equation is

```text
0=alpha^2 kappa D_4(p,r).                               (13)
```

Consequently

```text
D_4(p,r)!=0
  ==> no arbitrary target two-form has nonzero constant pullback
      for any 0!=H in t^2 C[t].                         (14)
```

This closes every degree-at-most-eight ordinary zero-second-debt critical
fold except the parabola `(11)`.

## 4. The first zero-fourth-debt candidate

At `r=0`, equation `(11)` gives `p=-5717/729`.  The corresponding polynomial
is

```text
Q_dagger=
 (22868x^8-89583x^6+2916x^5+123684x^4-5832x^3
   -63530x^2+2916x-2187)/2916.                         (15)
```

It has the fixed retained values and slopes and satisfies both retained debt
conditions

```text
D_2=0,                     D_4=0.                       (16)
```

This is a **candidate**, not a survivor theorem.  Nothing here proves that
`Q_dagger` has an actual target-ring lift with `J_0=1,J_1=J_2=0`, let alone
that such a lift reaches `J_4=0`.  The decisive next exact test is the
THM-3642 actual-lift system with `Q_6` replaced by `Q_dagger`.

## 5. Scope

The theorem is complete only for the degree-at-most-eight principal slope
packet `(3)`.  It does not cover

- higher-degree ordinary zero-second-debt polynomials;
- other projective slope charts in THM-3641;
- nonordinary tangent collisions;
- `H'(0)!=0` mixed pairs;
- another compiler or arbitrary planar polynomial maps.

On the parabola `(11)`, both a positive lift and the next obstruction remain
**OPEN**.  No Keller pair is constructed, and `JC(2)` remains **OPEN**.

## 6. Exact verification

```bash
python3 04-computation/jc2_russell_cylinder_degree8_fourth_debt_parabola_thm3677.py
python3 -O 04-computation/jc2_russell_cylinder_degree8_fourth_debt_parabola_thm3677.py
```

Normal and optimized transcripts are byte-identical to the stored output.
The exact companion reports `122` active gates and has no Python assertion
statements.
