---
id: THM-2927
title: "General-width flagged Macaulay leading coefficient"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  For every width M>=3 and every 0<a<b<M, the fixed
  20Q+10C+6F Macaulay determinant from THM-2921 has exact degree
  58M-36 and strictly positive leading coefficient.  Its pure-power
  chart factors as
  3*u0^14*v0^4*det2(u,v)^2*det3(u,v,w)^24.
  The two flag minors are strictly positive on factorial response rows
  by Cauchy--Binet and generalized-Vandermonde total positivity.
  Consequently the fixed chart is eventually positive at every width
  and support.  Finite-depth nonvanishing and arbitrary-width SFC(4)
  remain open.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
related:
  - THM-2830-consecutive-three-slot-factorial-moment-atomic-orientation
  - THM-2909-nonconsecutive-response-row-tp3-arc
  - THM-2921-diameter-four-nonconsecutive-macaulay-newton-closure
  - THM-2922-diameter-five-macaulay-newton-atlas
  - THM-2924-diameter-six-macaulay-newton-atlas
script: 04-computation/gmc_general_width_flagged_macaulay_leading_thm2927.py
output: 05-knowledge/results/gmc_general_width_flagged_macaulay_leading_thm2927.out
script_sha256: a473ffeb2725d333226f80274034f1ddc0cdd6e408a87067545b8972ff98a5ca
output_sha256: d701ed613a95761c66ad92ff6f6a39dc2a0d4382f846304f665af56c3d067d59
constructor_dependency_sha256: 42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64
hash_basis: LF-normalized bytes
---

# THM-2927 -- general-width flagged Macaulay leading coefficient

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Statement

Fix integers

```text
M>=3,                         0<a<b<M.                 (1)
```

For the factorial functional `L(s^j)=j!`, eliminate the mean on the
translated support `(n,n+a,n+b,n+M)` and scale the quadratic, cubic and
quartic moment forms by the PROVED THM-2925 denominators.  Call the
forms `Q~,C~,F~`.

Use the degree-seven Macaulay chart with the fixed row set

```text
0,...,19;
21,...,29,35;
36,...,41,                                             (2)
```

that is, `20Q+10C+6F`, and write its determinant as

```text
P_(M,a,b)(n).                                          (3)
```

Then

```text
deg P_(M,a,b)=58M-36,                                  (4)
```

and its leading coefficient is strictly positive.

More precisely, put

```text
R_j(d)=j^M-j^d,
r_j=(R_j(0),R_j(a),R_j(b)),               j=2,3,4.    (5)
```

The leading coefficient in `(4)` is

```text
3 R_2(0)^14 R_3(0)^4
  det[r_2|_(0,a);r_3|_(0,a)]^2
  det[r_2;r_3;r_4]^24,                                (6)
```

which is positive.

Consequently, for every fixed `(M,a,b)`, the same Macaulay chart is
positive for all sufficiently large integer depths `n`.  This last
threshold is not made uniform here.

## 2. Top scaled forms

For tensor order `j`, the normalized factorial entry at offsets
`d_1,...,d_j` is

```text
T_j(d_1,...,d_j;n)
 =(jn+1)_(sum_i d_i)/prod_i(n+1)_(d_i).                (7)
```

At infinity,

```text
T_j(d_1,...,d_j;n) -> j^(sum_i d_i).                  (8)
```

The mean-difference inclusion sum therefore has limiting form

```text
[sum_(d in {0,a,b}) x_d(j^d-j^M)]^j
 =(-1)^j(r_j dot x)^j.                                (9)
```

THM-2925 proves that the scaled order-`j` form has degree at most
`(j-1)M-1`, with monic scaling denominator.  Hence `(9)` is its
top-degree coefficient form.

The determinant uses ten cubic rows, so the ten minus signs in `(9)`
cancel.  Its candidate top coefficient is therefore the selected
Macaulay determinant of

```text
(r_2 dot x)^2,             (r_3 dot x)^3,
(r_4 dot x)^4.                                      (10)
```

## 3. The flagged pure-power identity

Let `u,v,w` be three row vectors in a characteristic-zero field.  Form
the degree-seven Macaulay matrix for

```text
(u dot x)^2,                 (v dot x)^3,
(w dot x)^4,                                        (11)
```

and select the rows `(2)`.  Its determinant is

```text
Delta(u,v,w)
 =3 u_0^14 v_0^4
   (u_0v_1-u_1v_0)^2
   det[u;v;w]^24.                                    (12)
```

Here is an exact flag-elimination proof.  The determinant has
multidegree

```text
deg_u Delta=40,                  deg_v Delta=30,
deg_w Delta=24.                                         (13)
```

The universal degree-`(2,3,4)` resultant divides every maximal
degree-seven Macaulay minor.  On pure powers it is

```text
Res((u dot x)^2,(v dot x)^3,(w dot x)^4)
 =det[u;v;w]^(2*3*4).                                  (14)
```

After removing `(14)`, the residual chart factor is independent of
`w` and has bidegree `(16,6)` in `(u,v)`.

Order the monomials first by `x_0` degree and then by `x_1` degree, as
in the constructor.  Block elimination of the residual flag map gives

| flag boundary | residual multiplicity |
|---|---:|
| `u_0=0` | `14` |
| `v_0=0` | `4` |
| `u_0v_1-u_1v_0=0` | `2` |

The three factors already exhaust bidegree `(16,6)`.  On the dense flag
chart the residual determinant is therefore

```text
c u_0^14 v_0^4(u_0v_1-u_1v_0)^2.                     (15)
```

For completeness, the block multiplicities can also be read from the
associated graded chart: at generic points of the three flag divisors,
the selected matrix ranks are respectively `29`, `34`, and `35`; the
missing pure-power pivots enter quadratically.  At

```text
u=(1,b,0),                    v=(d,1,0),
w=(0,0,1),                                             (16)
```

the remaining triangular block has determinant

```text
3d^4(bd-1)^26.                                        (17)
```

Taking `b=0,d=1` gives `c=3`.  Because both sides are polynomial,
the identity extends from the dense flag chart to all `u,v,w`, proving
`(12)`.

This explains why the stored 36-by-36 minor is stable: it is a Borel
flag semi-invariant, not an arbitrary numerical chart.

## 4. Total positivity of the response flags

Rewrite `(5)` as

```text
R_j(d)=(j-1) sum_(k=d)^(M-1) j^k.                     (18)
```

Let

```text
A_(j,k)=(j-1)j^k,                  j in {2,3,4},
H_(k,d)=1_(k>=d).                                      (19)
```

Then the response matrix is the product `AH`.  Cauchy--Binet gives

```text
det[r_2|_(0,a);r_3|_(0,a)]
 =sum_(k_1<a<=k_2)
   det[A_(2,3),(k_1,k_2)],                             (20)
```

and

```text
det[r_2;r_3;r_4]
 =sum_(k_1<a<=k_2<b<=k_3)
   det[A_(2,3,4),(k_1,k_2,k_3)].                      (21)
```

Every summation range is nonempty by `(1)`.  Each determinant on the
right is strictly positive: after removing the positive row factors
`j-1`, it is a generalized Vandermonde determinant at

```text
2<3<4,                       k_1<k_2<k_3,              (22)
```

equivalently an ordinary Vandermonde times a positive Schur polynomial.
Thus both flag minors in `(6)` are positive.  Also `R_2(0),R_3(0)>0`.
Equation `(12)` now proves `(6)`, so the upper degree bound from
THM-2925 is attained and `(4)` follows.

## 5. Consequences and sharp boundary

The result removes two possible all-width obstructions:

1. the fixed chart never loses its top degree through cancellation; and
2. its eventual sign never reverses.

It does **not** rule out a finite positive-depth root of
`P_(M,a,b)`.  Therefore it does not prove one-chart persistence at every
depth, arbitrary-width SFC(4), SFC(5), or arbitrary-radial GMC(2).
Widths four, five and six remain separately closed by THM-2921,
THM-2922 and THM-2924.  The next decisive target is positivity of the
lower deformation coefficients, or a finite Pluecker atlas that
switches charts across any finite-depth wall.

The strict inequalities `0<a<b<M` are load-bearing.  At a repeated
offset, one of the generalized Vandermonde/step-incidence minors
vanishes, exactly matching support collapse.

## 6. Exact verification

The companion:

- symbolically evaluates the two-parameter flag slice `(16)`;
- checks `(12)` at `40` deterministic generic integer triples;
- reproduces the three boundary ranks and one full-rank control;
- verifies the Cauchy--Binet identities, positivity, and `(12)` on all
  `220` supports with `3<=M<=12`; and
- records `4,004` positive TP2 and `5,005` positive TP3 summands.

These checks audit the universal algebraic proof; they do not replace
its quantifiers.  Run

```text
python 04-computation/gmc_general_width_flagged_macaulay_leading_thm2927.py
python -O 04-computation/gmc_general_width_flagged_macaulay_leading_thm2927.py
```

Normal and optimized executions byte-match the stored output with the
declared LF-normalized hashes.

**QED (candidate pending independent hostile audit).**
