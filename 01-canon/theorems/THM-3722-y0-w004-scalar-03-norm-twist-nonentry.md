---
id: THM-3722
title: "W004 scalar-03 norm-twist nonentry"
status: >
  PROVED + VERIFIED-EXACT.  Every W004 placement in the named all-scale family
  with weights P=(1-3n,1-2n,1), Q=(-n-2,-2,n-2,3n-2) and scalar fibre
  03+21 is Darboux-empty in the y=0 collision ring.  The tail creates a
  quadratic norm twist; coprime-power, negative-Laurent, differential-system,
  and degree gates close every branch, with n=1,2 handled exactly.  This
  closes the complete named family, not all of W004, the full 3x4 cell,
  general quartic C3 data, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion checks the W004 fibre
  word, upper transport, quadratic primitive, shifted norm, zero-shift
  systems in every residue modulo 7, exceptional n=5 determinant/scalar, n=2
  compatibility, and n=1 Euler boundary.  An independent derivation checked
  every sign, constant, UFD/coprime-power step, Laurent exponent, determinant,
  pullback, degree floor, and exceptional-scale identity.  Normal and
  optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3720-y0-w003-final-family-charge-coordinate-closure
script: 04-computation/jacobian_y0_w004_scalar03_norm_twist_thm3722.py
output: 05-knowledge/results/jacobian_y0_w004_scalar03_norm_twist_thm3722.out
script_sha256: 4ba3fde1275fbbd307814a645cc605d557ac6a27453d8786bf32298ebf651ef5
output_sha256: 8beb8aae4e7e9e945fbb6a7a9df5b931ac5227004c1fbb9392e6283d232454a0
hash_basis: LF-normalized bytes
---

# THM-3722 -- a complete W004 norm-twist family is empty

**PROVED + VERIFIED-EXACT.**  Work over `C` in the THM-3696 collision ring.
Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`; primes mean `d/db`.

## 1. Exact family and singleton bases

The W004 ray has gaps `(1,2,1,1,2)` and fibres

```text
00; 01+10; 02+11; 12+20; 03+21; 13+22; 23.            (2)
```

We close the following actual-support placement at every `n>=1`:

```text
wt(P)=(1-3n,1-2n,1),
wt(Q)=(-n-2,-2,n-2,3n-2),
scalar fibre 03+21.                                    (3)
```

For `n>=3`, set

```text
R=2n-1,
delta=gcd(3n-1,n+2)=gcd(n+2,7),
alpha=(3n-1)/delta,              beta=(n+2)/delta.      (4)
```

The two singleton rows give nonzero constants `a,c,d,t` and nonconstant
`H,K in C[b]` with

```text
f0=aH^alpha,           g0=cH^beta,
f2=dK,                 g3=tK^(3n-2).                   (5)
```

THM-3696 membership gives `h|H` and `b|K`.  Put

```text
U=HK^delta,                     E_delta=H'K+delta HK'.  (6)
```

Both `U` and `K` are nonconstant, and `E_delta` is nonzero by its leading
coefficient.  Write `M=f1`, `L=g1`, and `N=g2`.

## 2. Upper transport and the quadratic norm

Put `P=K^R M`.  The zero double `13+22` integrates exactly to

```text
N=K^(n-2)(lambda+rho P),
rho=t(3n-2)/d.                                         (7)
```

The next zero double `12+20`, multiplied by `K^(n+1)`, is the derivative of

```text
(n-2)lambda P +rho(3n-3)P^2/2-cdU^beta.                (8)
```

At a root of `K`, both `P=K^RM` and `U=HK^delta` vanish, so the integration
constant is zero:

```text
(n-2)lambda P +rho(3n-3)P^2/2=cdU^beta.                (9)
```

This is the first genuine quadratic charge cover in the W00x analysis.  The
W003 final family collapsed to the rational base `C(U)` in THM-3720; here the
zero rows initially retain a norm-twist sheet.

## 3. Nonzero shift: every beta>=2 twist is impossible

Suppose `lambda!=0`.  Complete the square in `(9)`.  With nonzero constants
`k,C`, it becomes

```text
S^2=k^2+CU^beta.                                       (10)
```

If `beta>=2`, then

```text
(S-k)(S+k)=CU^beta.                                    (11)
```

The two factors on the left are coprime because their difference is the
nonzero scalar `2k`.  Unique factorization over `C[b]` therefore gives, after
absorbing scalar beta-th roots,

```text
S-k=V^beta,                    S+k=W^beta.              (12)
```

Thus `W^beta-V^beta=2k`.  Factoring the left side over `C`, every factor
`W-zeta V` must be a unit.  Two distinct beta-th roots of unity then force
`V`, and hence `W,U`, to be constant, a contradiction.

The only possible nonzero-shift exception has `beta=1`.  Since `delta` divides
`7` and `n>=3`, this occurs exactly at

```text
n=5,                 delta=7, alpha=2, beta=1, R=9.    (13)
```

Section 5 closes it.

## 4. Zero shift: polynomiality sees a forbidden Laurent charge

Now let `lambda=0`.  Equation `(9)` says `P^2=CU^beta`.  Put

```text
g=gcd(2,beta).
```

Prime valuations in the UFD `C[b]` give a nonconstant polynomial `W` and
nonzero scalars such that

```text
U=u0 W^(2/g),                   P=p0 W^(beta/g).        (14)
```

Set

```text
Z=K^2L,                         V=Z'/U'.                (15)
```

Logarithmic differentiation of `(14)` gives
`P'/P=(beta/2)U'/U`.  The two remaining zero doubles `02+11` and `01+10`
become

```text
R V-beta Z/U
 +a rho alpha(3n-2)U^(alpha-1)/2=0,                   (16)

a alpha U^(alpha-1)(delta U V-2Z)
 +c beta(3n-4)U^(beta-1)P/2=0.                        (17)
```

Their determinant is `-a alpha(3n-4)U^(alpha-1)`, so they give the exact
algebraic solution

```text
Z= -[a rho alpha delta(3n-2)/(2(3n-4))]U^alpha
   +[c beta R/(2a alpha)]U^(beta-alpha)P.               (18)
```

The first term is polynomial.  Under `(14)`, the second nonzero term is a
scalar times

```text
W^((3beta-2alpha)/g),
3beta-2alpha=(8-3n)/delta<0.                            (19)
```

It is a negative power of the nonconstant polynomial `W`, contradicting
`Z in C[b]`.  Hence every zero-shift tail is empty.

## 5. The exceptional n=5 rational norm

At `(13)`, equation `(9)` makes `U` a quadratic polynomial in `P`:

```text
U=[6rho P^2+3lambda P]/(cd),            rho=13t/d.      (20)
```

Let `theta=dU/dP` and now put `V=Z'/P'`.  The rows `02+11` and `01+10` are a
linear system over `C(P)` for `(V,Z)`.  Its determinant is

```text
2a alpha U^(alpha-1)(delta U-RP theta).                 (21)
```

The last factor is

```text
delta U-RP theta=-6P(11rho P+lambda)/(cd),              (22)
```

so `(21)` is nonzero.  Thus `Z in C(P)`.  Since `Z in C[b]` and `P` is a
nonconstant polynomial, the polynomial-pullback lemma from THM-3720 gives

```text
Z=Q(P) in C[P].                                         (23)
```

The scalar row then has the exact form

```text
K^(-1)P'[26at U theta-dQ_P].                            (24)
```

If it equalled one, then `P'Psi(P)=K` for a polynomial `Psi`.  But
`P=K^9M` and `M` is a nonzero polynomial, so

```text
deg(P'Psi(P))>=deg P-1>=9deg K-1>deg K.                (25)
```

This contradiction closes `n=5` when `lambda!=0`; its zero-shift case was
already closed by Section 4.

## 6. The n=2 differential-compatibility boundary

At `n=2`, put `P=K^3M` and `rho=4t/d`.  The upper row gives

```text
N=lambda+rho P,                                        (26)
```

and the next row gives

```text
P=sigma U^2,                    sigma^2=2cd/(3rho).     (27)
```

For `Z=K^2L` and `V=Z'/U'`, the two lower zero rows are

```text
3UV-4Z+10a rho U^5=0,
5aUV-10aZ+4c sigma U=0.                                (28)
```

Solving them algebraically yields

```text
Z=6c sigma U/(5a)-5a rho U^5,
V=8c sigma/(5a)-10a rho U^4.                           (29)
```

But differentiating the first expression in `(29)` and using `V=Z'/U'`
would require

```text
0=dZ/dU-V=-2c sigma/(5a)-15a rho U^4.                  (30)
```

The right side is a nonconstant polynomial with nonzero constant and leading
coefficients.  This closes `n=2` before the scalar row.

## 7. The n=1 proportional-row boundary

At `n=1`, the upper zero double integrates to

```text
K(tM-dN)=constant.                                     (31)
```

Evaluation at a root of `K` kills the constant, so `N=tM/d`.  In the next
zero double, the hidden bracket has equal weights and proportional
coefficients, hence vanishes.  What remains is

```text
-3cdH^2 E_1!=0,                                        (32)
```

a contradiction.

## 8. Scope

Sections 3--7 close every scale of the named W004 placement `(3)`.  This is
one complete all-scale family and the first explicit algebraization/no-go for
the quadratic charge twist.  Other W004 placements, W005--W006, arbitrary
`3 x 4` supports, unrestricted quartic C3/cofactor data, and `JC(2)` remain
open.

Run

```bash
python3 -B 04-computation/jacobian_y0_w004_scalar03_norm_twist_thm3722.py
python3 -B -O 04-computation/jacobian_y0_w004_scalar03_norm_twist_thm3722.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
