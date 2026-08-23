---
id: THM-3727
title: "W004 dual scalar-13 ternary Kummer nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every W004
  placement in the named all-scale family with weights
  P=(-n-2,-2,2n-2), Q=(1-4n,1-3n,1-2n,1) and scalar fibre 13+22 is
  Darboux-empty in the y=0 collision ring.  Its endpoint gcd divides 9, and
  the end rows create a negative-charge Kummer sheet that is nonpolynomial
  except at n=2; the remaining scale is differentially incompatible.  This
  does not close all of W004, the full 3x4 cell, general quartic C3 data, or
  JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion checks the W004 fibre
  word, singleton arithmetic in every residue modulo 9, both end-row
  transports, the negative-charge split, the complete n=2 middle-row Cramer
  system and compatibility polynomial, and the n=1 boundary.  Normal and
  optimized runs byte-match the frozen transcript.  An independent
  derivation checked every sign and constant, the rational-function UFD
  step for delta=1,3,9, negative-charge coprimality, every n=2 Cramer term,
  and the n=1 reduction.
depends_on:
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3724-y0-w004-scalar-13-kummer-twist-nonentry
script: 04-computation/jacobian_y0_w004_dual_scalar13_ternary_kummer_thm3727.py
output: 05-knowledge/results/jacobian_y0_w004_dual_scalar13_ternary_kummer_thm3727.out
script_sha256: 041fdbd779d656f3c9db4cd92d63c8ea7d65ac680651678811c4462e9903e0c5
output_sha256: 690d9e6c9ce80d6200bf6399a0132acd4aff659dbc2d33f80e66dd0e7d14432e
hash_basis: LF-normalized bytes
---

# THM-3727 -- the dual W004 ternary Kummer family is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over `C`
in the THM-3696 collision ring.  Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`; primes mean `d/db`.

## 1. Exact family and endpoint bases

The W004 supports `(0,n,3n)` and `(0,n,2n,4n)` have fibres

```text
00; 01+10; 02+11; 12+20; 03+21; 13+22; 23.            (2)
```

We close the opposite scalar-`13+22` orientation from THM-3724:

```text
wt(P)=(-n-2,-2,2n-2),
wt(Q)=(1-4n,1-3n,1-2n,1),
scalar fibre 13+22,                           n>=1.    (3)
```

First suppose `n>=2`.  Set

```text
delta=gcd(n+2,4n-1)=gcd(n+2,9),
alpha=(n+2)/delta,              beta=(4n-1)/delta,
T=2n-2,                         C=3n-1,
S=2n-3,                         m=2alpha-beta=(5-2n)/delta.
                                                               (4)
```

The singleton rows `00` and `23` give nonzero constants `a,c,d,t` and
nonconstant `H,K in C[b]` with

```text
f0=aH^alpha,          g0=cH^beta,
f2=dK^T,              g3=tK.                           (5)
```

THM-3696 membership gives `h|H` and `b|K`.  Put

```text
U=HK^delta,           Pi=K^2f1,
L=g1.                                                   (6)
```

## 2. The end rows and the negative Kummer charge

The upper row `03+21` is exactly

```text
(K^C L)'=[at/(dT)](U^alpha)'.                          (7)
```

At a root of `K`, the integrated left side and `U^alpha` both vanish.  Thus
the integration constant is zero.  With `rho=at/(dT)`,

```text
L=rho K^(-C)U^alpha=rho H^alpha K^(-S).                (8)
```

Polynomiality of the right side is already a nontrivial divisibility
condition, but we do not need to estimate it.  Substitution into the lowest
row `01+10` gives

```text
delta U Pi'-2Pi U'=D U^m U',
D=-a rho alpha S/(c beta).                             (9)
```

Since `delta*m-2=-S`, put

```text
A=a rho alpha/(c beta),                 F=Pi-AU^m.     (10)
```

Then `(9)` becomes

```text
delta U F'-2F U'=0,
hence                         F^delta=kappa U^2        (11)
```

in the rational function field `C(b)`, for some `kappa in C`.

## 3. Every n>=3 sheet is nonpolynomial

The number `delta` divides `9`, so it is odd and `gcd(delta,2)=1`.  If
`n>=3`, then

```text
delta*m=5-2n<0.                                       (12)
```

If `kappa=0`, equation `(11)` gives `Pi=AU^m`, a nonzero negative power of
the nonconstant polynomial `U`.  This contradicts `Pi in C[b]`.

Suppose `kappa!=0`.  Prime valuations in the UFD `C[b]` applied to `(11)`
show, after scalar rescaling over `C`, that

```text
U=v^delta,                    F=Bv^2,          B!=0,   (13)
```

for a nonconstant polynomial `v`.  Indeed, coprimality of `delta` and `2`
first forces every prime valuation of `U` to be divisible by `delta`, and
then forces `ord_p(F)>=0` at every finite prime and the displayed valuation
of `F`; hence the initially rational `F` is the polynomial `Bv^2`.
Equations `(10)` and `(12)` now read

```text
Pi=Av^(5-2n)+Bv^2.                                    (14)
```

Writing `s=2n-5>0`, the numerator of `(14)` over `v^s` is
`A+Bv^(s+2)`, which is coprime to `v` because `A!=0`.  Thus `(14)` is not a
polynomial.  This closes every `n>=3` before the middle or scalar rows.

## 4. The exceptional n=2 sheet is differentially incompatible

At `n=2`, one has

```text
delta=1,       alpha=4,       beta=7,       m=1,
Pi=AU+BU^2                                             (15)
```

for an arbitrary `B in C`, including zero.  Put `N=g2`, `Z=K^3N`, and
`Y=Z'/U'`.  The two middle rows `12+20` and `02+11` become

```text
2Pi Y-3Pi_U Z=14cdU^6,
4a(UY-3Z)+rho(-5U Pi_U+8Pi)=0.                        (16)
```

Their determinant in `(Y,Z)` is

```text
-12aAU !=0.                                            (17)
```

Cramer's rule gives

```text
Z= U[3A^2rho+ABrho U-2B^2rho U^2+28acdU^5]/(6Aa),

Y= [3A^2rho+4ABrho U-4B^2rho U^2+56acdU^5]/(4Aa).
                                                               (18)
```

But the necessary identity `dZ/dU=Y` would require

```text
0=dZ/dU-Y
 =[-3A^2rho-8ABrho U+168acdU^5]/(12Aa).                (19)
```

The constant and fifth-power coefficients in `(19)` are nonzero.  Since
`U` is nonconstant, `(19)` is impossible for every value of `B`.

## 5. The n=1 boundary

When `n=1`, the singleton `23` has weights `(0,1)`.  Its zero bracket is

```text
W_(0,1)(f2,g3)=f2'g3=0.                                (20)
```

The domain and actual support force `f2'=0`, so the weight-zero piece `f2`
is a literal scalar.  Delete it without changing the bracket.  The remaining
pair has at most two graded pieces in the first output and four in the
second, and is impossible by THM-3583.

## 6. The dyadic/ternary endpoint pairing

THM-3724 and this theorem close the two opposite W004 scalar-`13+22`
placements by the same end-row operation, but their Euclidean sidecars are

```text
THM-3724:  gcd(n-3,8)  in {1,2,4,8},
THM-3727:  gcd(n+2,9)  in {1,3,9}.                     (21)
```

Thus one orientation produces the dyadic Kummer tower and the other the
ternary tower.  This is an exact structural comparison of these two named
families, not a general modular-group theorem.  The asymmetry is useful:
the ternary orientation reverses the charge far enough that polynomiality
alone closes every generic scale, while the dyadic orientation needs the
middle-row compatibility of THM-3724.

## 7. Scope and reproduction

Sections 2--5 close every scale of the named W004 placement `(3)`.  Other
W004 placements, W005--W006, arbitrary `3 x 4` supports, unrestricted
quartic C3/cofactor data, and `JC(2)` remain open.

Run

```bash
python3 -B 04-computation/jacobian_y0_w004_dual_scalar13_ternary_kummer_thm3727.py
python3 -B -O 04-computation/jacobian_y0_w004_dual_scalar13_ternary_kummer_thm3727.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
