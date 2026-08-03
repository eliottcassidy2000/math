---
id: THM-3311
title: "Alternating-sextic HFC(3) exclusion"
status: >
  PROVED + VERIFIED-EXACT.  Every homogeneous degree-six polynomial in the
  S3-sign eigenspace is uniquely
  V(A e1^3+B e1 e2+C e3).  Its odd simplex moments vanish by symmetry, but
  the exact even-moment equations I2=I4=I6=0 have no nonzero complex
  projective solution.  The A=1 ideal is the unit ideal over Q; an independent
  resultant certificate has degrees 8 and 12 and coprime reductions modulo
  17.  The A=0 line is excluded separately.  The pair I2=I4=0 does retain a
  genuine irreducible degree-eight survivor, so I6 is load-bearing.  This
  excludes one homogeneous symmetry cell only and does not prove HFC(3) or
  FC(3).
source: root/creative-passports/2026-08-03
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
related:
  - THM-3304-fourier-dirichlet-kernel-and-alternating-quintic-hfc3-exclusion
  - THM-3310-degree-four-cyclic-eigenspace-on-the-triangle
script: 04-computation/factorial_hfc3_alternating_sextic_thm3311.py
output: 05-knowledge/results/factorial_hfc3_alternating_sextic_thm3311.out
script_sha256: 743ae0a937bb18965ee2d6469235acd47e4fcd86bd5ddb52bf10a0efbc31b75e
output_sha256: 07a42600cb8ed306adfe5be16dbf28ecb4b199b5575b17aedf16303857541cdb
hash_basis: LF-normalized bytes
---

# THM-3311 -- alternating-sextic HFC(3) exclusion

**PROVED + VERIFIED-EXACT.**

Let `Delta_2` be the barycentric triangle `x+y+z=1`, equipped with the standard
simplex coordinate measure `dA=dx dy` after `z=1-x-y` (a harmless nonzero
normalization does not affect zeros).  Write

```text
V=(x-y)(y-z)(z-x),
e1=x+y+z,       e2=xy+yz+zx,       e3=xyz.                (1)
```

This theorem closes the next `S_3`-sign cell after the alternating quintic of
THM-3304.  It deliberately stays inside the homogeneous subclass identified
by THM-3018 and MISTAKE-350.

## 1. Normal form and automatic moments

Every alternating polynomial is divisible by `V`.  After division, a
homogeneous alternating sextic leaves a symmetric cubic, whose basis is
`e1^3,e1 e2,e3`.  Thus the cell is exactly

```text
G=V(A e1^3+B e1 e2+C e3),       [A:B:C] in P^2(C).        (2)
```

A transposition preserves the triangle measure and reverses `G`, so every odd
moment of `G` is zero.  Put `I_m=int_Delta G^m dA`.  Two independent exact
integration routes give primitive numerator forms `Q_2,Q_4,Q_6` with

```text
I_2=Q_2/25225200,
I_4=Q_4/133989577722000,
I_6=Q_6/6909970509968180054400,                            (3)
```

and term counts `6,15,28`.  In particular,

```text
Q_2=15015 A^2+6006 AB+182 AC+637 B^2+42 BC+C^2.           (4)
```

The companion freezes every coefficient of `Q_4,Q_6`; their hashes are

```text
Q2  f9abb5b86dc5f66a0f7e3264f60972178565c0a8c3d7db52a0626072d0449c70
Q4  4d5c643af1ce370e7642b2f89356ba8e275382c925c21fc03507af610125bf99
Q6  bf3cf3c34f410f1ca6731b8be828bb46efff08402c5d6510aa647c3e9214fec4.
```

One route expands homogeneous barycentric monomials and uses the factorial
Dirichlet integral.  The other substitutes `z=1-x-y` first and integrates a
sparse polynomial over `0<=y<=1-x, 0<=x<=1`.  They agree coefficientwise.

## 2. The affine chart is empty

On `A=1`, regard the three forms as polynomials in `B,C`.  Exact Groebner
elimination over `Q` gives

```text
<Q_2(1,B,C),Q_4(1,B,C),Q_6(1,B,C)> = <1>.                (5)
```

There is also a compact independent certificate.  Eliminate `B` from
`(Q_2,Q_4)` and `(Q_2,Q_6)` and take primitive resultants
`R_24(C),R_26(C)`.  Their degrees are `8,12`.  Modulo `17` they retain those
degrees, have leading coefficients `12,9`, and satisfy

```text
gcd(R_24,R_26)=1 in F_17[C].                              (6)
```

Full-degree reduction makes `(6)` a characteristic-zero coprimality
certificate: a common primitive factor over `Q` would retain a nonconstant
leading term modulo `17` and divide both reductions.

## 3. The projective boundary is empty

On `A=0,B=1`, exact univariate elimination already gives

```text
Res_C(Q_2,Q_4)=14827264400 != 0.                          (7)
```

The only point missed by that chart is `[A:B:C]=[0:0:1]`; its three primitive
coefficients are `(1,1,12)`.  Equations `(5)--(7)` therefore cover all of
`P^2(C)` and prove

```text
I_2=I_4=I_6=0       implies       A=B=C=0.                (8)
```

Hence no nonzero alternating sextic can satisfy all positive-power factorial
moments.

## 4. First failure boundary

The sixth moment is not decorative.  On `A=1`, `R_24` is irreducible of
degree eight over `Q`.  The eliminated `B`-degrees of `Q_2,Q_4` are constant,
so its roots represent genuine common zeros of `I_2,I_4` over an algebraic
extension.  Adding `I_6` changes this nonempty degree-eight survivor to the
unit ideal `(5)`.

Thus the sharp result of this computation is not a positivity claim: complex
cancellation survives two even moments, and the first load-bearing extra
condition in this cell is `I_6`.

## 5. Scope and reproduction

Proved: the normal form `(2)`, the two-route moment forms, projective
emptiness `(8)`, and the stated two-moment survivor.  Open: every other
degree-six representation cell, arbitrary homogeneous `HFC(3)`, and full
nonhomogeneous `FC(3)`.

Run

```text
python 04-computation/factorial_hfc3_alternating_sextic_thm3311.py
python -O 04-computation/factorial_hfc3_alternating_sextic_thm3311.py
```

and compare LF-normalized bytes with the declared output.  Normal and
optimized outputs are byte-identical.  The script contains no assertion node
or floating literal and uses exact integer, rational, modular, and Groebner
arithmetic only.

**QED.**
