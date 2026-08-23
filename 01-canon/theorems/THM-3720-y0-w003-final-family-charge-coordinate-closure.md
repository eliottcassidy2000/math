---
id: THM-3720
title: "W003 final-family charge-coordinate closure"
status: >
  PROVED + VERIFIED-EXACT.  The final residual W003 placement family D is
  Darboux-empty in the y=0 collision ring at every scale n>=3.  The triple
  row puts one hidden coefficient in a single charge coordinate U=HK^delta;
  the next two zero rows form a nonsingular differential system that forces
  the other hidden coefficient into C[U].  The scalar row then fails by
  divisibility or degree.  Combined with THM-3695, THM-3714, THM-3715, and
  THM-3717, this closes the complete W003 ray, not the full 3x4 cell, general
  quartic C3 data, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion checks the final
  placement census, triple primitive, charge polynomial, both differential
  rows, determinant, C(U) solution, scalar factor, and both gcd(n+2,7)
  outcomes in all seven residue classes.  An independent derivation checked
  every sign, the crucial delta*U*V coefficient, constant removal, determinant
  boundary, polynomial-pullback lemma, and final divisibility/degree split.
  Normal and optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-3714-y0-w003-scalar-fibre-02-two-orientation-euler-nonentry
  - THM-3715-y0-w003-odd-lowest-scalar-two-orientation-half-charge-nonentry
  - THM-3717-y0-w003-scalar-triple-square-law-euler-nonentry
related:
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
script: 04-computation/jacobian_y0_w003_final_charge_coordinate_thm3720.py
output: 05-knowledge/results/jacobian_y0_w003_final_charge_coordinate_thm3720.out
script_sha256: c70393f1020fb5b74a80e60beaffb228ff778a1019e537fd02e0e68b8c7c295c
output_sha256: 2ac646a64fff91c5801e118d0ca1034e2be71383ac828ceb0856755b776b34b2
hash_basis: LF-normalized bytes
---

# THM-3720 -- the final W003 family is empty

**PROVED + VERIFIED-EXACT.**  Work over `C` in the THM-3696 collision ring.
Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`; primes mean `d/db`.

## 1. The exact final placement

After the gates combined in THM-3715 and the family-C closure THM-3717, the
only W003 placement left is

```text
wt(P)=(1-3n,1-2n,1),
wt(Q)=(-n-2,-2,n-2,2n-2),                   n>=3,       (2)
```

with scalar fibre `13+21`.  Its fibres are

```text
00; 01+10; 02+11; 03+12+20; 13+21; 22; 23.             (3)
```

Every displayed coefficient is active; a zero coefficient falls under the
already closed smaller-support gates.

Set

```text
R=2n-1,
delta=gcd(3n-1,n+2)=gcd(n+2,7),
alpha=(3n-1)/delta,             beta=(n+2)/delta.       (4)
```

The singleton rows give nonzero constants `a,c,d,e,t` and nonconstant
`H,K in C[b]` with

```text
f0=aH^alpha,          g0=cH^beta,
f2=dK,                g2=eK^(n-2),       g3=tK^(2n-2).
                                                               (5)
```

THM-3696 membership gives `h|H` and `b|K`.  Thus

```text
E_delta=H'K+delta HK'                                  (6)
```

is nonzero, with leading multiplier `deg H+delta deg K`.
Write `M=f1` at weight `-R` and `L=g1` at weight `-2`.

## 2. The triple creates one charge coordinate

The zero triple `03+12+20`, multiplied by `K^(n+1)`, is exactly the derivative
of

```text
e(n-2)MK^R +(2n-2)atH^alpha K^(3n-1)
             -cdH^beta K^(n+2).                        (7)
```

At a root of `K`, every term in `(7)` vanishes.  Hence its integration
constant is zero.  Define

```text
U=HK^delta,
A=-2(n-1)at/[e(n-2)],             B=cd/[e(n-2)],
P(U)=AU^alpha+BU^beta.                                  (8)
```

Then `(7)` is equivalent in `C(b)` to

```text
M=K^(-R)P(U).                                           (9)
```

The right side is a Laurent presentation of the original polynomial `M`;
no polynomiality is discarded.  Its two monomials have the same charge:

```text
q+R=delta p.                                            (10)
```

## 3. The two lower zero rows are a nonsingular system

Put

```text
Z=K^2L in C[b],                    V=Z'/U' in C(b).      (11)
```

The polynomial `U` is nonconstant, so `U'` is nonzero.  Substitute `(5),(9)`
in the zero double `02+11`; after cancelling its common power of `K` and
dividing by `U'`, it becomes

```text
R P V-2P_U Z+ae(n-2)alpha U^(alpha-1)=0.               (12)
```

The lowest zero double `01+10` similarly becomes

```text
a alpha U^(alpha-1)(delta U V-2Z)
 +c beta U^(beta-1)(RP-delta U P_U)=0.                 (13)
```

Equations `(12),(13)` are a linear system over `C(U)` for `(V,Z)`.  Its
determinant is

```text
Delta_sys
 =2a alpha U^(alpha-1)(delta U P_U-RP).                (14)
```

It does not vanish.  Indeed

```text
delta U P_U-RP = An U^alpha+B(3-n)U^beta,              (15)
alpha-beta=(2n-3)/delta>0.                              (16)
```

At `n=3`, the first term in `(15)` is nonzero and the second vanishes.  At
`n>3`, the two active terms have distinct exponents.  Thus `(14)` is nonzero
for every `n>=3`, and Cramer's rule gives

```text
V in C(U),                       Z in C(U).             (17)
```

This is the decisive rigidity step: the zero rows eliminate all transverse
branch freedom of the last hidden coefficient.

## 4. Polynomial pullback lemma

We also know `Z in C[b]`.  For any nonconstant `U in C[b]`,

```text
C[b] intersect C(U)=C[U]             inside C(b).       (18)
```

Here is the short proof.  Write a rational function in lowest terms as
`q(T)=A(T)/B(T)`.  If `B` is nonconstant, choose a complex root `lambda` of
`B`; coprimality gives `A(lambda)!=0`.  The nonconstant polynomial
`U(b)-lambda` has a complex root, so `q(U(b))` has a pole and cannot lie in
`C[b]`.  Therefore `B` is constant.  This proves `(18)`.

Applying `(18)` to `(17)`, there is a polynomial `Q` with

```text
Z=Q(U),                         L=K^(-2)Q(U).           (19)
```

## 5. The scalar row fails by divisibility or degree

The scalar double `13+21`, using `(9),(19)`, now has the exact form

```text
W_(-R,R-1)(K^(-R)P(U),tK^(R-1))
 +W_(1,-2)(dK,K^(-2)Q(U))

 =K^(-1)U'[(R-1)tP_U-dQ_U]
 =K^(delta-2)E_delta Phi(U),                            (20)
```

where `Phi(U)=(R-1)tP_U-dQ_U` is a polynomial.  Since `delta` divides `7`,
there are only two cases.

If `delta=7`, the right side of `(20)` is

```text
K^5 E_7 Phi(U),                                         (21)
```

which is either zero or divisible by the nonunit `K`; it cannot equal one.

If `delta=1`, equality to one in `(20)` would give

```text
U'Phi(U)=K.                                             (22)
```

If `Phi=0`, this is already impossible.  Otherwise, with `m=deg Phi>=0`,

```text
deg(U'Phi(U))=(m+1)deg U-1
              >=deg U-1
               =deg H+deg K-1
               >deg K,                                 (23)
```

because `h|H` gives `deg H>=2`.  Equation `(22)` is impossible.  Thus the
scalar row never equals one.

## 6. W003 is complete

This proves that the final residual family

```text
P=(1-3n,1-2n,1), Q=(-n-2,-2,n-2,2n-2), n>=3,
scalar fibre 13+21                                      (24)
```

is empty.  THM-3715's exhaustive residual theorem left only families C and D;
THM-3717 closes C and this theorem closes D.  Together with THM-3695 and
THM-3714, every post-parity W003 placement at every scale is Darboux-empty in
the y=0 collision ring.

This is a complete ray theorem, not a closure of W004--W006, every `3 x 4`
word, the unrestricted quartic C3/cofactor lane, or `JC(2)`.

Run

```bash
python3 -B 04-computation/jacobian_y0_w003_final_charge_coordinate_thm3720.py
python3 -B -O 04-computation/jacobian_y0_w003_final_charge_coordinate_thm3720.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
