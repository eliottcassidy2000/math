---
id: THM-2943
title: "Width-seven/eight two-chart Macaulay resultant closure"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  For every n>=0, first-window SFC(4) holds on every
  translated four-slot support of exact width seven or eight.  The
  original and stable-mutated 20Q+10C+6F Macaulay determinants have a
  raw polynomial gcd which is the genuine resultant times the common
  positive flag content and selected-basis factors.  In all 36 new
  width-seven/eight families this raw gcd has nonnegative coefficients
  and positive constant term, while the residual chart cofactors are
  coprime.  Width six is replayed only as an inherited control.  No
  width-nine, arbitrary-width, shifted-window, or SFC(5) claim is made.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
related:
  - THM-2849-four-slot-first-window-macaulay-box
  - THM-2924-diameter-six-macaulay-newton-atlas
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
  - THM-2944-width-nine-ten-two-chart-macaulay-resultant-closure
script: 04-computation/gmc_width_seven_eight_two_chart_resultant_closure_thm2943.py
output: 05-knowledge/results/gmc_width_seven_eight_two_chart_resultant_closure_thm2943.out
script_sha256: d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba
output_sha256: 99b262e58a7e338294ad473629fdd057e415779c5e12ec0706eff5ed2d855049
factorization_dependency_sha256: 85ce40de2aa777b8af091dfec934de68beeb42676bea3b48dbf815990a51e0e9
hash_basis: LF-normalized bytes
---

# THM-2943 -- width-seven/eight two-chart Macaulay resultant closure

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^j)=j!.          (1)
```

Fix

```text
M in {7,8},              0<a<b<M,              n>=0.   (2)
```

For every nonzero polynomial

```text
H=c_0 s^n+c_1 s^(n+a)+c_2 s^(n+b)+c_3 s^(n+M),       (3)
```

at least one of

```text
L(H),                  L(H^2),                  L(H^3), L(H^4) (4)
```

is nonzero.  Equivalently, first-window SFC(4) holds on all

```text
C(6,2)+C(7,2)=15+21=36                                (5)
```

translated four-slot supports of exact width seven or eight.

## 2. Mean elimination and the common Macaulay model

Put

```text
f_j=s^j/j!,                         L(f_j)=1.           (6)
```

On the mean-zero hyperplane write

```text
H=x(f_n-f_(n+M))
 +y(f_(n+a)-f_(n+M))
 +z(f_(n+b)-f_(n+M)).                                (7)
```

For moment order `m`, PROVED THM-2925 supplies the exact positive
denominator

```text
D_(m,M)
 =[prod_(j=1)^(M-1)(n+j)^(m-1)](n+M)^(m-2).          (8)
```

After clearing `(8)`, let `Q,C,F` be the resulting homogeneous forms
of degrees `2,3,4` in `x,y,z`.  The scaling is nonzero for every
integer `n>=0`, so it preserves their common projective zero set.

Use the degree-seven Macaulay map

```text
Phi_n:S_5 direct_sum S_4 direct_sum S_3 -> S_7,
Phi_n(A,B,D)=A Q+B C+D F.                             (9)
```

The inherited stored transpose has `46` rows and `36` columns.  Both
charts use the quadratic and cubic rows

```text
Q rows: 0,...,19,
C rows: 21,...,29,35.                                (10)
```

The original quartic rows are

```text
J0 local  ={0,1,2,3,4,5},
J0 global ={36,37,38,39,40,41},                      (11)
```

while the stable two-row mutation is

```text
J1 local  ={0,3,4,5,6,7},
J1 global ={36,39,40,41,42,43}.                      (12)
```

Call the two selected determinants `Delta_0(n),Delta_1(n)`.
THM-2925 gives the uniform degree invoice

```text
deg Delta_i <=58M-36.                                 (13)
```

## 3. The two universal flag factors

Use the coefficient notation of PROVED THM-2942.  Its universal
factorization gives

```text
Delta_0=q200^6*c300*K*R,
Delta_1=q200^5*c300*P_alt*R,                          (14)
```

where

```text
R=Res(Q,C,F)                                          (15)
```

is the genuine ternary resultant,

```text
K
 =c120*q200^2-c210*q110*q200
  -c300*q020*q200+c300*q110^2,                       (16)
```

and

```text
P_alt
 = c012*q020*q200^2
  -c021*q011*q200^2
  -c210*q002*q020*q200
  +c210*q011^2*q200
  +c300*q002*q020*q110
  -c300*q011^2*q110.                                 (17)
```

These are global polynomial identities.  They remain valid on the
boundary where the selected quotient basis used to derive them
degenerates.

## 4. The two distinct exact gates

For one normalized support `(M,a,b)`, put

```text
g=gcd(K,P_alt) in Z[n]                                (18)
```

with positive leading coefficient, and write

```text
K=g K_0,                         P_alt=g P_0.          (19)
```

The exact replay proves the **flag-separation gate**

```text
gcd(K_0,P_0)=1,
gcd(q200*K_0,P_0)=1.                                 (20)
```

The second identity in `(20)` retains the extra original-chart power
of `q200`; omitting it would leave a possible common chart wall.

Now normalize the raw polynomial gcd

```text
G_raw=gcd(Delta_0,Delta_1)                            (21)
```

to have positive leading coefficient.  Equations `(14),(19),(20)`
and direct exact division give the association in `Q[n]`

```text
G_raw ~ q200^5*c300*R*g.                              (22)
```

The positive rational associate in `(22)` is recorded exactly for
every family.  The replay then proves the independent
**common-resultant-content gate**

```text
G_raw(0)>0,
every coefficient of G_raw is nonnegative.           (23)
```

Consequently

```text
G_raw(n)>0                         for every n>=0.     (24)
```

The separation between `(20)` and `(23)` is load-bearing.  THM-2942
already proves that the variable flag coordinates are separated
through width twenty, but that alone does not rule out the common
resultant wall.  Here `(22),(23)` prove that the entire common factor,
including `R`, is nonzero on the integer depth ray.

Indeed, if `R(n_0)=0` for some integer `n_0>=0`, then `(22)` would give
`G_raw(n_0)=0`, contradicting `(24)`.  Therefore

```text
Res(Q,C,F)(n)!=0                    for every n>=0.    (25)
```

## 5. Projective consequence

By the defining property of the ternary resultant, `(25)` says that
`Q,C,F` have no common point in `P^2`.  Restoring the eliminated mean
in `(7)`, a nonzero polynomial `(3)` therefore cannot satisfy all four
equalities in `(4)`.  This proves the statement.

Equivalently, at every depth at least one of the two Macaulay charts is
nonzero.  This is stronger than merely changing a variable Pluecker
coordinate: the common selected-basis/resultant content has itself
been shown positive on the depth ray.

## 6. The exact 46-family replay

The companion deliberately includes exact width six as an independent
control:

```text
width 6: C(5,2)=10 families,
width 7: C(6,2)=15 families,
width 8: C(7,2)=21 families,
total:              46 families.                     (26)
```

Width six is already closed by PROVED THM-2924 and is not a new
consequence here.

For each support, the script interpolates both raw determinants from
the complete prefix required by `(13)`, then:

1. divides them by the exact universal factors in `(14)` and verifies
   that both quotients are the same resultant polynomial;
2. verifies `(20)` and the associated raw-gcd identity `(22)`;
3. checks every coefficient and the constant term in `(23)`;
4. compares both interpolants with direct `36`-by-`36` determinants at
   the three outside-grid depths

   ```text
   58M-35, 58M-34, 58M-33,                            (27)
   ```

   giving `46*3=138` paired direct controls; and
5. checks one exact three-term Grassmann--Pluecker exchange per
   support.  With `I={0,3,4,5}` and quartet `{1,2,6,7}`, the checked
   identity is

   ```text
   Delta_(I12) Delta_(I67)
   -Delta_(I16) Delta_(I27)
   +Delta_(I17) Delta_(I26)=0.                        (28)
   ```

The original reduced chart cofactor has both positive and negative
coefficients in

```text
width 6: 7,
width 7: 12,
width 8: 18,
total:   37                                             (29)
```

families.  Thus single-chart coefficientwise positivity genuinely
stops in most of this bank.  Equation `(20)` is the exact reason the
stable mutation removes every *common* chart obstruction; `(29)` does
not assert that each mixed polynomial has an integral root.

The complete family record has digest

```text
0449a7d9f1799da7af7be9e8a41e84d6d79c1aa3b68aa4c13bfd5eb09a342fa4.
                                                               (30)
```

## 7. Scope and stopping boundary

The proved-candidate scope is exactly:

```text
slots:        four distinct exponents;
widths:       exactly seven or eight;
window:       moments one through four;
depth:        every integer n>=0;
certificate:  original/stable two-chart resultant atlas.       (31)
```

The replay contains width six only as a control.  Width nine and ten
have a separately reserved follow-on, THM-2944, and are not silently
included here.  No arbitrary-width SFC(4), shifted moment window,
SFC(5), arbitrary-radial GMC(2), or Jacobian-conjecture conclusion is
claimed.

The finite result also does not establish a closed formula for
`deg g`; THM-2942 records a genuine width-eleven hostile to the naive
degree pattern.  The proof uses only the exact cellwise properties
`(20),(23)`.

## 8. Exact verification

The companion hash-pins the independently audited THM-2942
factorization artifact.  It uses only explicit `require` checks, so
optimized Python retains every truth-bearing test.  Run

```text
python 04-computation/gmc_width_seven_eight_two_chart_resultant_closure_thm2943.py
python -O 04-computation/gmc_width_seven_eight_two_chart_resultant_closure_thm2943.py
```

Normal and optimized executions must byte-match the stored output and
the declared LF-normalized hashes.

**Proof-complete candidate.**
