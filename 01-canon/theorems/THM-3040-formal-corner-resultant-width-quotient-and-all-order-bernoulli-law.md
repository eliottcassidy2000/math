---
id: THM-3040
title: "Formal factorial-resultant corner width quotient and the all-order Bernoulli law"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; independent immutable audit pending.  For
  every fixed three-lower-slot support and every width M above it, the intrinsic
  pure-resultant formal U=V=0 corner germ has exact normalized width quotient
  (1+Mt)^26(1+(M+1)t)^20.  More generally, the k-slot factorial system has
  exponents k!(H_k-1) and k!(k+1-2H_k).  Consequently THM-3030's C-corner
  Bernoulli--Faulhaber law holds at every jet order, without constructing the
  global P_j.  This is a coefficientwise t-adic corner, not a physical width,
  raw selected-chart identity, or wall-stripped-core statement.
source: kind-pasteur-2026-08-01-S?
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2946-full-macaulay-maximal-minor-gcd-and-chart-free-resultant
related:
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-3030-eighth-resultant-jet-and-corner-constant-closed-forms
script: 04-computation/gmc_formal_corner_resultant_width_quotient_thm3040.py
output: 05-knowledge/results/gmc_formal_corner_resultant_width_quotient_thm3040.out
script_sha256: 1428cb92776ce9b27839c8fe9fefb3d26c00cf6a3e895e052fa5678bfdc26d0e
output_sha256: b17fbaecd97e810fe938d6717fc5b21b3a6b74fc51c2947d9aeea49d2ccad404
hash_basis: LF-normalized bytes
---

# THM-3040 -- the formal corner quotient

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.**

THM-3030 found, through the first eight exact resultant log jets, the width
difference

```text
Delta_M((-1)^(j-1) j C_j)=26M^j+20(M+1)^j.           (1)
```

It correctly left the all-order continuation open.  The missing mechanism is
not a larger interpolation.  It is the multihomogeneous character of the
pure resultant after taking the coefficientwise formal corner.

## 1. Statement

Fix three distinct lower offsets

```text
0 <= a < b < c < M,
```

independent of `M`, and form the factorial quadratic, cubic, and quartic of
THM-2925 for the four-slot support `(a,b,c,M)`.  Put `t=1/n`.  Let
`R_M^C(t)` be the normalized **intrinsic pure-resultant** germ obtained by
taking the coefficientwise formal corner

```text
U=2^M=0,       V=3^M=0,       4^M=U^2=0.             (2)
```

Then `R_M^C(0) != 0`, and for every `M>c`, with one fixed standard resultant
sign convention,

```text
R_(M+1)^C(t) / R_M^C(t)
   =(1+Mt)^26 (1+(M+1)t)^20.                          (3)
```

For the first-gap support `(0,1,2,M)`, define

```text
C_j(M)=[t^j] log(R_M^C(t)/R_M^C(0)),
p_j^C(M)=(-1)^(j-1) j C_j(M).                         (4)
```

This is exactly THM-3030's `C_j=Q_j(M,0,0)` wherever its frozen global jet
`Q_j` has been constructed.  At every order `j>=1`, whether or not the global
three-variable `P_j` has been constructed,

```text
p_j^C(M+1)-p_j^C(M)=26M^j+20(M+1)^j,                 (5)

p_j^C(M)=46 sum_(s=1)^(M-1) s^j+20M^j+K_j,           (6)
```

where `K_j` is independent of `M` (and may depend on the fixed lower triple).

## 2. The coefficientwise corner kills every terminal face

For tensor order `r in {2,3,4}`, offsets `d_1,...,d_r`, and
`S=sum_i d_i`, THM-2925 uses

```text
T_r(d;n)=(rn+1)_S / product_i (n+1)_(d_i).            (7)
```

As a formal `t`-series,

```text
T_r(d;1/t)
 =r^S product_(q=1)^S(1+qt/r)
      / product_i product_(q=1)^(d_i)(1+qt).          (8)
```

Equivalently, with `P_k(N)=sum_(q=1)^N q^k`,

```text
T_r(d;1/t)=r^S exp sum_(k>=1) ((-1)^(k-1)t^k/k)
             (r^(-k)P_k(S)-sum_i P_k(d_i)).           (9)
```

Thus each fixed `t`-coefficient is `r^S` times a polynomial in the offsets.
If a nonempty inclusion face replaces `q>0` lower offsets by the terminal
offset `M`, then

```text
r^S=r^B (r^M)^q.                                      (10)
```

It lies in `(U)` for `r=2`, in `(V)` for `r=3`, and in `(U^2)` for `r=4`.
The corner (2) therefore kills every nonempty terminal face and retains only
the empty face.  The three surviving unscaled forms, denoted
`G_2,G_3,G_4`, use only `(a,b,c)` and are independent of `M`.

This assertion is deliberately **coefficientwise in the `t`-adic
completion**.  It is not a claim that a finite rational function containing a
terminal offset is globally divisible by `r^M` after arbitrary operations at
infinity.  Resultant and logarithm coefficients use only finitely many jets,
so the coefficientwise corner is the lawful object.

## 3. The corner logarithm exists

At `t=0`, multinomial expansion gives

```text
G_r(0;x)=(r^a x_0+r^b x_1+r^c x_2)^r.                (11)
```

The generalized Vandermonde matrix `(r^d)`, with rows `r=2,3,4` and columns
`d=a,b,c`, has nonzero determinant (strict total positivity).  Hence

```text
Res(G_2,G_3,G_4)|_(t=0)
 =det((r^d))^(2*3*4) !=0.                             (12)
```

For `(a,b,c)=(0,1,2)`, the determinant is `2`, so the constant is `2^24`.
Therefore the low resultant is an `M`-independent unit series.  Corner
evaluation commutes with the pure resultant and its formal logarithm.

## 4. Resultant character and the exact quotient

Normalize THM-2925's denominator clearers by their leading powers of `n`:

```text
d_(r,M)(t)
 =product_(s=1)^(M-1)(1+st)^(r-1) (1+Mt)^(r-2).      (13)
```

Section 2 gives the corner forms exactly as `d_(r,M)G_r`.  The ternary
resultant of degrees `(2,3,4)` has coefficient multidegrees

```text
(3*4,2*4,2*3)=(12,8,6).                              (14)
```

Consequently, writing `R_low(t)=Res(G_2,G_3,G_4)`,

```text
R_M^C(t)
 =R_low(t) d_(2,M)^12 d_(3,M)^8 d_(4,M)^6
 =R_low(t)
    product_(s=1)^(M-1)(1+st)^46 (1+Mt)^20.           (15)
```

The two exponents are the pairings

```text
(12,8,6).(1,2,3)=46,       (12,8,6).(0,1,2)=20.       (16)
```

Increasing the width promotes the old terminal factor from weight `20` to
interior weight `46`, a net exponent `26`, and attaches the new terminal with
weight `20`.  Equation (3) follows.  Equivalently, before `t` normalization,

```text
R_M^C(n)=R_low(n)
  (Gamma(n+M)/Gamma(n+1))^46 (n+M)^20.                (17)
```

No Barnes function is needed: Faulhaber is the logarithmic ghost of this one
Gamma flag.

## 5. All-order Bernoulli slots

Taking `[t^j] log` in (3) proves (5), and telescoping proves (6).  Faulhaber's
formula now gives, for every `m>=1` and every `j` with `2m-1<j`,

```text
[M^(j-2m+1)] C_j
 =(-1)^(j+1) 46 B_(2m)/(2m)!
    (j-1)(j-2)...(j-2m+2).                           (18)
```

Since `sign(B_(2m))=(-1)^(m+1)`, this is THM-3030's corrected
`(-1)^(j+m)` sign law, now at every order, with

```text
c_m^C=46 |B_(2m)|/(2m)!       for every m>=1.         (19)
```

All even nonterminal slots vanish because `B_(2m+1)=0` for `m>=1`.  The two
top slots are

```text
[M^(j+1)]C_j=(-1)^(j+1)46/(j(j+1)),
[M^j]C_j=(-1)^j 3/j.                                  (20)
```

The preregistered tests of THM-3030 are therefore theorems:

```text
c_5^C=23/23950080,       [M]C_10=-23/66,
c_6^C=15893/653837184000,[M]C_12=15893/16380.          (21)
```

The factor `691` from `B_12=-691/2730` really is the first break of the
observed reduced numerator `23`.

## 6. General `k`-slot boundary character

The argument is not special to four slots.  Let `k>=3`, fix `k-1` distinct
lower offsets, and take the `k-1` homogeneous factorial forms of degrees
`r=2,...,k` in `k-1` variables.  Give the terminal exponential `r^M` its own
formal symbol and take the augmentation corner that kills every positive
terminal monomial.  The coefficientwise argument of section 2 leaves one
fixed low form at each order.

The resultant multidegree of the order-`r` form is

```text
mu_r=product_(2<=s<=k,s!=r) s = k!/r.                 (22)
```

Therefore the universal quotient is

```text
R_(M+1)^C/R_M^C
 =(1+Mt)^A_k (1+(M+1)t)^B_k,                          (23)

A_k=sum_(r=2)^k mu_r       =k!(H_k-1),
B_k=sum_(r=2)^k (r-2)mu_r  =k!(k+1-2H_k).             (24)
```

The interior density and normalized resultant degree are

```text
I_k=A_k+B_k=k!(k-H_k),
deg R_M^C=I_k M-A_k.                                  (25)
```

For consecutive lower offsets `0,...,k-2`, the low response determinant is

```text
det(r^a)_(r=2..k,a=0..k-2)=product_(j=1)^(k-2) j!,    (26)
```

so the low resultant constant is the `k!`-th power of (26), in particular a
unit.  General fixed distinct lower offsets use the same strict generalized-
Vandermonde argument as section 3.

The first cases are

| slots `k` | `A_k` | `B_k` | `I_k` |
|---:|---:|---:|---:|
| 3 | 5 | 2 | 7 |
| 4 | 26 | 20 | 46 |
| 5 | 154 | 172 | 326 |

Writing `p_j=(-1)^(j-1)jC_j` gives universally

```text
p_j(M)=I_k sum_(s=1)^(M-1)s^j+B_k M^j+K_j,            (27)
```

so the all-order odd-slot magnitude is

```text
I_k |B_(2m)|/(2m)!,                                   (28)
```

and every even nonterminal slot vanishes.  This is an abstract factorial-
resultant theorem.  It does not assert a global higher-SFC closure, a useful
sign for its wall-stripped core, or existence of a single physical chart at
all widths.

## 7. What was and was not used

This proof uses only THM-2925's all-width denominator law and the intrinsic
resultant multihomogeneity isolated by THM-2942/2946.  It does **not** use:

- the unproved continuation wall of THM-2997;
- nonvanishing of one selected `36x36` Macaulay chart at physical widths;
- construction or interpolation of `P_9`, `P_10`, or any later global `P_j`;
- positivity or ULC of the wall-stripped norm core.

The result is the pure resultant's coefficientwise formal corner.  It is not
a physical integer width satisfying `2^M=3^M=0`, not the raw selected minor
(which carries extraneous flag factors), and not a statement about the
wall-stripped core.  A separately primitive-normalized resultant can alter
(3) by an `n`-independent rational unit; all positive-degree log jets and
(5)--(21) are unchanged.  Fixing the standard resultant convention removes
even that harmless unit.

The hypothesis that `(a,b,c)` is fixed is load-bearing.  If a lower offset
moves with `M`, the surviving `R_low` need not be width-independent and the
quotient acquires an additional transport term.

## 8. Exact evidence

The companion independently checks:

- all `344` inclusion faces in the three forms (`31` survive, `313` terminal
  faces die at the corner);
- the character ledger `(46,20,26)` and quotient at every `M=3..40`;
- the general formulas (23)--(26) for every `k=3..10`, including the exact
  `5/2`, `26/20`, and `154/172` controls;
- the leading determinant `2` and resultant constant `2^24`;
- strict generalized-Vandermonde positivity on all `286` triples
  `0<=a<b<c<=12`;
- the exact logarithmic recurrence and Bernoulli slots through order `32`;
- every nonconstant coefficient in all eight frozen THM-3030 C-corner tables,
  with their eight independent integration constants retained as controls;
- the `P_10` and `P_12` predictions in (21).

Reproduce with normal and optimized Python:

```text
python 04-computation/gmc_formal_corner_resultant_width_quotient_thm3040.py
python -O 04-computation/gmc_formal_corner_resultant_width_quotient_thm3040.py
```

Both runs equal the stored ten-line transcript byte-for-byte.

**QED, subject only to the pending independent immutable-file audit.**
