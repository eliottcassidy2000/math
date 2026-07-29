---
id: THM-2854
title: "Four-slot shifted-window Macaulay box"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY HOSTILE-AUDITED.
  For every four-slot support with top exponent at most eight, each of the
  shifted windows k=1 and k=2 has no nonzero common zero of its four
  factorial moment forms.  The sharp Macaulay maps in degrees 11 and 15
  have full target ranks 364 and 816 modulo 1000003 in all 252 cells.  A
  nonzero modular maximal minor lifts to characteristic zero.  The result
  proves neither larger supports nor arbitrary shifted SFC(4).
source: root/gmc-four-slot-shifted-macaulay-2026-07-28
depends_on: []
related:
  - THM-2836-sfc3-arbitrary-support-shifted-window-census
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2849-four-slot-first-window-macaulay-box
script: 04-computation/sfc4_shifted_windows_macaulay_box_thm2854.py
output: 05-knowledge/results/sfc4_shifted_windows_macaulay_box_thm2854.out
script_sha256: f5f98b19202b603309c56a30e8179adc2d68859f74a526b6e70e5732b59df699
output_sha256: b97309e98af2bf013636f310a039f84f98c3213aaffa426679654bd0e3360b43
hash_basis: LF-normalized bytes
---

# THM-2854 -- four-slot shifted-window Macaulay box

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY
HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L:C[s] -> C,                    L(s^n)=n!,
f_a=s^a/a!,                     L(f_a)=1.                (1)
```

For every support

```text
0<=a_0<a_1<a_2<a_3<=8,                                  (2)
```

every `k in {1,2}`, and every

```text
H=c_0f_a0+c_1f_a1+c_2f_a2+c_3f_a3,                      (3)
```

the four equations

```text
L(H^(k+1))=L(H^(k+2))=L(H^(k+3))=L(H^(k+4))=0           (4)
```

force `H=0`.  The same conclusion holds for coefficients in the
unnormalized monomial basis, because multiplication of each coordinate by
the nonzero scalar `a_i!` is a diagonal bijection.

There are exactly

```text
binom(9,4)*2=126*2=252                                  (5)
```

cells in the stated universe.  The proof below supplies an exact
characteristic-zero certificate for every cell.

## 2. The four moment forms

For an order `m`, write `P_m(c)=L(H^m)`.  It is a homogeneous form of
degree `m` in `R=Q[c_0,c_1,c_2,c_3]`, with coefficient formula

```text
P_m(c)=
 sum_(e_0+...+e_3=m)
 [m!/prod_i e_i!]
 [(sum_i e_i a_i)!/prod_i(a_i!)^e_i]
 prod_i c_i^e_i.                                        (6)
```

The companion constructs each form in two separate ways modulo

```text
p=1000003:                                               (7)
```

1. enumerate the exponent counts in `(6)` and insert the multinomial
   coefficient;
2. enumerate all `4^m` ordered tensor words and accumulate their monomial
   counts.

The resulting coefficient dictionaries agree for every order `2<=m<=6`
on all 126 supports.  A third path constructs the forms over `Q` and
compares their values with direct univariate polynomial powering on three
separated supports, two coefficient vectors, and five orders, giving 30
exact controls.

The script verifies by trial division that `p` is prime.  Moreover

```text
p>6*8=48,                                                (8)
```

so none of the factorials used in the rational coefficients meets `p`.
In particular every denominator in `(6)` is invertible modulo `p`.

## 3. Sharp Macaulay certificates

For fixed `k`, put

```text
(d_1,d_2,d_3,d_4)=(k+1,k+2,k+3,k+4),
D=sum_i d_i-3=4k+7.                                     (9)
```

Consider the degree-`D` multiplication map

```text
Phi_k:
 direct_sum_(i=1)^4 R_(D-d_i) -> R_D,
 (A_1,...,A_4) |-> sum_i A_i P_(k+i).                  (10)
```

The degree in `(9)` is sharp for a four-form complete intersection: the
Artinian quotient's last possible degree is
`sum_i(d_i-1)=D-1`.  In monomial bases, the two matrix shapes are

```text
k=1: D=11,      589 rows by 364 columns;
k=2: D=15,     1325 rows by 816 columns.                (11)
```

For every support in `(2)`, exact sparse elimination over `F_p` gives

```text
rank(Phi_1)=364,                 rank(Phi_2)=816.        (12)
```

Thus in every cell a maximal minor is nonzero modulo `p`.  Clear the
rational coefficient denominators by a number coprime to `p`.  The
corresponding integer minor remains nonzero modulo `p`, hence is a nonzero
integer.  Therefore the same map has full target rank over `Q`, and

```text
(P_(k+1),P_(k+2),P_(k+3),P_(k+4))_D=R_D.               (13)
```

If a common projective zero `[c_0:c_1:c_2:c_3]` existed over `C`, every
element of the ideal component in `(13)` would vanish there.  Some
coordinate `c_j` is nonzero, but `c_j^D` belongs to `R_D`, a
contradiction.  Hence the projective common-zero set is empty, proving
`(4)`.

## 4. Positive and hostile controls

The exact census is surrounded by controls that distinguish full rank
from a merely large rank.

- For each `k=0,1,2`, the monomial complete intersection
  `(c_0^(k+1),c_1^(k+2),c_2^(k+3),c_3^(k+4))` has full rank in the sharp
  degree.
- A synthetic system with all four forms vanishing at
  `[0:0:0:1]` has ranks `114,340,756`, strictly below target dimensions
  `120,364,816`.
- On the actual support `(0,2,5,8)`, deleting the highest moment gives
  exactly the same three deficient ranks.  The respective deficits

  ```text
  6,24,60=(k+1)(k+2)(k+3)                              (14)
  ```

  match the expected Bézout lengths of the first-three-form complete
  intersections in the tested degree.  Thus the fourth form is
  load-bearing in this control rather than an eliminator artifact.
- A pivot engine using the opposite column order agrees on the supports
  `(0,1,2,3)`, `(0,2,5,8)`, and `(5,6,7,8)` in both shifted windows.
- A SHA-256 digest binds every support, every coefficient vector, both
  window labels, and all 252 ranks:

  ```text
  21c3c85cf89184651e302965bafef49274e22f94545414c463087baddea0d2a9.
  ```

## 5. Evidence and independent audit

Run

```text
python 04-computation/sfc4_shifted_windows_macaulay_box_thm2854.py
python -O 04-computation/sfc4_shifted_windows_macaulay_box_thm2854.py
```

The standard-library companion uses explicit exceptions rather than
Python `assert` for every truth-bearing gate.  Normal, optimized, and
stored runs agree byte for byte.

The independent hostile audit rederived formula `(6)`, the sharp degree
and matrix dimensions in `(9)--(11)`, the modular-minor lift, and the
projective-emptiness implication.  It inspected the two coefficient
constructors, opposite-pivot controls, synthetic common-point system, and
the exact deficit identity `(14)`, and replayed the final stored
certificate.  No remaining mathematical or evidence defect was found.

## 6. Boundary

This theorem proves exactly

```text
four slots,     0<=a_0<a_1<a_2<a_3<=8,     k in {1,2}. (15)
```

It proves no support with top exponent above eight, no window `k>=3`,
no arbitrary shifted SFC(4), no SFC(5), and no new unbounded GMC case.
THM-2849 gives the complementary bounded first-window result through top
exponent fifteen.  Combining two finite boxes does not make a uniform
theorem.
