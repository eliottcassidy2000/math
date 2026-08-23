---
id: THM-3728
title: "W004 scalar-12 anchor-20 charge-normalization nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The W004
  all-scale placement with scalar fibre 12+20 and arm anchor 20 is
  Darboux-empty in the y=0 collision ring.  Even scales die by the inherited
  parity gate.  At odd scales at least three, the zero rows normalize the
  last Kummer remainder into C(V), where its required odd valuation is
  impossible; the zero Kummer sheet has a noncancellable top coefficient.
  The exceptional scale n=1 collapses to a constant-product derivative.
  Together with the audited tail census, this closes the entire stabilized
  W004 tail.  It does not close exceptional small-scale placements, W005--W006,
  general quartic C3 data, or JC(2).
source: root + lrc14-connection-miner + jc-quartic-c3-construction / 2026-08-22
audit: >
  PASS.  An independent hostile derivation checked the endpoint transports,
  complete Euler/Kummer remainder, charge-normalized second-row factor,
  rational-field inference, both Kummer sheets, the n=1 root and derivative
  collapse, the even-scale parity gate including n=2, and the stabilized-tail
  census dependency.  Normal and optimized runs byte-match the frozen
  transcript; the recorded hashes match.
depends_on:
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3722-y0-w004-scalar-03-norm-twist-nonentry
  - THM-3724-y0-w004-scalar-13-kummer-twist-nonentry
  - THM-3727-y0-w004-dual-scalar-13-ternary-kummer-nonentry
  - THM-3733-y0-w004-scalar-12-persistent-arm-pair-nonentry
  - THM-3735-y0-w004-scalar-02-persistent-arm-pair-nonentry
  - THM-3739-y0-w004-scalar-01-tail-pair-nonentry
script: 04-computation/jacobian_y0_w004_scalar12_anchor20_charge_normalization_thm3728.py
output: 05-knowledge/results/jacobian_y0_w004_scalar12_anchor20_charge_normalization_thm3728.out
script_sha256: 6d23c9f025acede93f7d5aab76fa9af66c71af832804199e6c4ed6cce73be282
output_sha256: cf19f17f4a86d9cbde1571fdb3e92e7168d500894149348156ac38b374a005d3
hash_basis: LF-normalized bytes
---

# THM-3728 -- the final stabilized W004 tail family is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
`C` in the THM-3696 collision ring.  Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`.  The W004 ray has fibre word

```text
00; 01+10; 02+11; 12+20; 03+21; 13+22; 23,            (2)
```

and the placement considered here is

```text
wt(P)=(1-3n,1-2n,1),
wt(Q)=(-2,n-2,2n-2,4n-2),                 n>=1.        (3)
```

The scalar fibre is `12+20`; address `20=(1,-2)` is the persistent arm.
At every even `n`, the negative singleton component `00` has gcd
`gcd(3n-1,2)=1`.  Its weight-`-2` arm coefficient therefore has common-base
exponent two, and THM-3613 rejects the placement.  Sections 1--3 treat every
odd scale.

## 1. Odd n at least three: endpoint transports

Fix odd `n>=3` and put

```text
alpha=(3n-1)/2,              m=(n-1)/2.               (4)
```

The singleton endpoint rows have gcds `2` and `1`, so for nonzero constants
`a,c,d,t` and nonconstant `H,K in C[b]`,

```text
f0=aH^alpha,       g0=cH,       f2=dK,
g3=tK^(4n-2).                                           (5)
```

THM-3696 gives `h|H` and `b|K`.  In particular

```text
V=HK^2,                         U=V^alpha              (6)
```

are nonconstant.  Write `R=f1,L=g1,M=g2` and set

```text
mu=(4n-2)t/d.                                          (7)
```

The rows `03+21` and `13+22` integrate completely to

```text
L=lambda K^(n-2)+a mu H^alpha K^(4n-3),
M=kappa K^(2n-2)+mu R K^(4n-3),                       (8)
```

with arbitrary `lambda,kappa in C`.  Subtracting the displayed particular
solutions leaves respectively the derivatives of `L/K^(n-2)` and
`M/K^(2n-2)`, so no homogeneous branch is missing.

## 2. The first lower row and its Kummer remainder

Put

```text
E2=H'K+2HK',                    P0=H^(alpha-1)K^(n-2).
                                                               (9)
```

The row `01+10=0` is exactly

```text
2HR'-(2n-1)RH'
 =(a alpha/c)H^(alpha-1)K^(n-3)
   [lambda(n-2)+a mu(4n-3)U]E2.                      (10)
```

Its complete solution in `C(b)` is

```text
R=P0[(a alpha/c)(lambda+a mu U)+f],                   (11)
f^2 V^(n-2)=chi,                         chi in C.     (12)
```

Indeed, the two terms in `(11)` are direct Euler particulars.  If `S=P0f`
is the remainder, its homogeneous equation is

```text
2HS'-(2n-1)SH'=0,
S^2=chi H^(2n-1).                                     (13)
```

Dividing by `P0^2` gives `(12)` because
`P0^2=H^(3n-3)K^(2n-4)`.

## 3. Charge normalization closes both Kummer sheets

Define

```text
r=R/P0,
D=lambda(n-2)+a mu(4n-3)U,
G=lambda(n-2)(2n-1)+a mu(4n-3)(5n-2)U.               (14)
```

Use `(10)` to eliminate `R'` from `02+11=0`.  The complete row factors as

```text
02+11 = [H^(alpha-2)K^(2n-5)E2/2]
 [G r+(a alpha/c)D^2+4a kappa alpha(n-1)V].           (15)
```

The prefactor is nonzero: `V'=KE2`, and `V` is nonconstant.  Also `G` is a
nonzero affine polynomial in `U`.  Consequently `(15)` puts `r in C(V)`.
The explicit particular in `(11)` already belongs to `C(V)`, hence

```text
f in C(V).                                             (16)
```

If `chi!=0`, take the `V=0` valuation in the rational function field `C(V)`.
Equation `(12)` gives

```text
2 ord_0(f)+(n-2)=0,                                   (17)
```

which is impossible because `n-2` is odd.

If `chi=0`, then `f=0`.  Substitute

```text
r=(a alpha/c)(lambda+a mu U)                          (18)
```

into the second factor of `(15)`.  As a polynomial in `V`, its `U^2`
coefficient is

```text
(a alpha/c)(a mu)^2(4n-3)(9n-5),                     (19)
```

which is nonzero in characteristic zero.  The remaining term involving
`kappa` has degree only one in `V`, whereas `U^2=V^(2alpha)` and
`alpha>=4`.  Thus the zero Kummer sheet is impossible as well.  The scalar
row is never reached: the zero rows already contradict each other.

## 4. The exceptional scale n=1

At `n=1`, the endpoint rows give

```text
f0=aH,          g0=cH,          f2=dK,          g3=tK^2.
                                                               (20)
```

Put `mu=2t/d`.  The two upper rows integrate to

```text
KL=a mu HK^2+lambda,              M=mu KR+kappa.      (21)
```

Because `K` is nonconstant, it has a root; evaluation of the first identity
there gives `lambda=0`.  Hence `L=a mu HK`.  The row `02+11=0` is now

```text
0=(a mu/K^2)(HK^3R)'.                                 (22)
```

Thus the nonzero polynomial `HK^3R` is constant.  This is impossible because
`H` and `K` are nonconstant and `R` is the retained nonzero coefficient.
This closes `n=1`.

## 5. Consequence and exact scope

THM-3739's audited tail census showed that `(3)` was the only stabilized
W004 tail placement not already captured by THM-3722, THM-3724, THM-3727,
THM-3733, THM-3735, and THM-3739.  Therefore the entire stabilized W004 tail
is empty.  Exceptional small-scale schemes not belonging to the stabilized
families, W005--W006, arbitrary `3 x 4` supports, unrestricted quartic
C3/cofactor data, and `JC(2)` remain open.

Run

```bash
python3 -B 04-computation/jacobian_y0_w004_scalar12_anchor20_charge_normalization_thm3728.py
python3 -B -O 04-computation/jacobian_y0_w004_scalar12_anchor20_charge_normalization_thm3728.py
```

Both commands must agree byte for byte with the frozen transcript.
