---
id: THM-2849
title: "Four-slot first-window Macaulay closure through exponent fifteen"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  four-slot support with top exponent at most fifteen,
  the first four factorial moments have no common nonzero coefficient
  vector.  After eliminating the mean, two independent exact constructors
  produce the quadratic, cubic, and quartic forms; a degree-seven Macaulay
  matrix has rank 36 modulo 1000003 in all 1,820 cells.  A nonzero modular
  maximal minor lifts to characteristic zero and proves projective
  emptiness.  The bounded result neither proves arbitrary SFC(4) nor any
  shifted moment window.
source: root/gmc-four-slot-macaulay-box-2026-07-28
depends_on: []
related:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2836-bounded-three-slot-strong-factorial-macaulay-census
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
script: 04-computation/sfc4_first_window_macaulay_box_thm2849.py
output: 05-knowledge/results/sfc4_first_window_macaulay_box_thm2849.out
script_sha256: d32bd79201bbb0ae3db9904b27869f7e112084777d1cc7d33b28a37d46734a61
output_sha256: 98afca5d00f25ff2fda409a58df96973e04e8f30e27d90764525a28470a1d528
hash_basis: LF-normalized bytes
---

# THM-2849 -- four-slot first-window Macaulay closure through exponent fifteen

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^n)=n!.          (1)
```

For every support

```text
0<=a_0<a_1<a_2<a_3<=15
```

and every

```text
H=c_0s^a_0+c_1s^a_1+c_2s^a_2+c_3s^a_3,               (2)
```

the equations

```text
L(H)=L(H^2)=L(H^3)=L(H^4)=0                           (3)
```

force `H=0`.  Equivalently, SFC(4) holds in the first moment window on all
`binom(16,4)=1820` supports with top exponent at most fifteen.

This is an exact finite theorem, not a statistical census.  Its proof is a
nonzero-minor certificate over a finite field which lifts to
characteristic zero.

## 2. Mean elimination and exact moment forms

Normalize

```text
f_a=s^a/a!,                         L(f_a)=1.           (4)
```

Rescaling the four coefficients by nonzero factorials does not change
projective nullity.  The equation `L(H)=0` therefore permits the unique
parametrization

```text
H=x(f_a0-f_a3)+y(f_a1-f_a3)+z(f_a2-f_a3).             (5)
```

Put

```text
Q=L(H^2),                 C=L(H^3),                 F=L(H^4).
```

These are homogeneous forms of degrees `2,3,4` in `R=Q[x,y,z]`.  Before
the substitution `c_3=-c_0-c_1-c_2`, their coefficients are given
explicitly by

```text
M_m(c)=
 sum_(e_0+...+e_3=m)
 [m!/prod_i e_i!]
 [(sum_i e_i a_i)!/prod_i(a_i!)^e_i]
 prod_i c_i^e_i.                                         (6)
```

The companion reconstructs each of `Q,C,F` in two independent ways:

1. substitute `c_3=-x-y-z` into `(6)` and expand its multinomials;
2. expand every ordered tensor
   `L(prod_j(f_a_i_j-f_a3))` through its `2^m` signed factorial terms.

The coefficient dictionaries agree exactly in every one of the 1,820
cells.

## 3. The degree-seven Macaulay certificate

Let `R_d` denote the degree-`d` part of `R`.  For each support form the
linear map

```text
Phi:R_5 direct_sum R_4 direct_sum R_3 -> R_7,
Phi(A,B,D)=AQ+BC+DF.                                    (7)
```

In monomial bases the script stores the transpose of the conventional
matrix for `Phi`, with

```text
dim R_5+dim R_4+dim R_3=21+15+10=46 rows,
dim R_7=36 columns.                                      (8)
```

The exact companion reduces the rational matrix modulo

```text
p=1000003>60.                                            (9)
```

It verifies primality of `p`, so every factorial denominator occurring in
the box is a unit.  For all 1,820 supports,

```text
rank_Fp Phi=36.                                         (10)
```

Hence some `36 by 36` minor is nonzero modulo `p`.  After clearing the
factorial denominators by a number prime to `p`, the corresponding integer
minor is not divisible by `p`, and in particular is nonzero.  Therefore

```text
rank_Q Phi=36,                 (Q,C,F)_7=R_7.           (11)
```

If a projective common zero `[x:y:z]` existed over `C`, every element of
`(Q,C,F)_7` would vanish there.  But at least one coordinate is nonzero,
and its seventh power belongs to `R_7`, contradicting `(11)`.  Thus the
projective variety is empty and `(3)` has no nonzero solution.

## 4. Exact controls

The companion checks more than the headline ranks.

- It enumerates the universe independently as `binom(16,4)=1820`.
- Both moment-form constructors agree cell by cell.
- A NumPy integer modular eliminator and a separate pure-Python eliminator
  agree on the consecutive cells `(0,1,2,3)`, `(1,2,3,4)`, the far cell
  `(0,5,10,15)`, and the deletion control.
- Deleting the quartic rows at `(1,2,3,4)` leaves rank `30<36`.
  Thus the test does not accidentally certify the false assertion that
  the quadratic and cubic alone have no common projective point.
- A SHA-256 digest binds every support, all three exact coefficient
  vectors, and its rank to the stored transcript.

Normal, optimized, and stored runs agree byte for byte; all truth-bearing
gates use explicit exceptions rather than Python `assert`.  An independent
hostile audit replayed both modes and the stored transcript, rederived the
modular-minor lift, and checked that full degree-seven span implies
projective emptiness.  It found no mathematical defect.

## 5. Boundary

This theorem proves only

```text
number of slots=4,               first window k=0,
maximum exponent<=15.                                    (12)
```

It does not prove unbounded SFC(4), any shifted window, SFC(5), or a
uniform symbolic multipole-separation inequality.  THM-2846 shows why the
quartic is genuinely load-bearing: on the support `{1,2,3,4}` a positive
two-cone plane already has a nonzero common quadratic/cubic null line.
THM-2848 isolates the unbounded geometric target as factorial
cubic--quartic multipole separation.  The finite Macaulay box is rigorous
evidence for that target, not a replacement for it.

The bounded conclusion is therefore proved with the scope in `(12)`.
