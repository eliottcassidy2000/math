---
id: THM-3735
title: "W004 scalar-02 persistent-arm pair nonentry"
status: >
  PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.  Both named
  W004 all-scale placements with scalar fibre 02+11 and a persistent arm at
  address 11 are Darboux-empty in the y=0 collision ring.  One orientation
  creates a harmless high Kummer sheet but forces second-order arm death;
  the opposite orientation kills every generic scale by arm order and its
  lone n=2 survivor by incompatible local valuations.  This is not proved
  canon until audit promotion, and it does not close all of W004, the full
  3x4 cell, general quartic C3 data, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion checks both absolute
  placements, endpoint gcd tables modulo 3, 8, and 7, both generic upper and
  lowest-row signs across complete residue windows, strengthened arm orders,
  the n=2 primitive and exact local-order gap, and both n=1 deletions.  Normal
  and optimized runs byte-match the frozen transcript.  Independent hostile
  audit remains open.
depends_on:
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3733-y0-w004-scalar-12-persistent-arm-pair-nonentry
script: 04-computation/jacobian_y0_w004_scalar02_persistent_arm_pair_thm3735.py
output: 05-knowledge/results/jacobian_y0_w004_scalar02_persistent_arm_pair_thm3735.out
script_sha256: 100465cfdcf3bfcb80788496255d3def51ef13575b15f0d713c21ddc18e952d7
output_sha256: c832ccec20cdebc55c196b7c0c29ba5cb39cdbf1f2d0f7c63a3932d3c38e4a26
hash_basis: LF-normalized bytes
---

# THM-3735 -- both persistent-arm W004 scalar-02 families are empty

**PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.**  Work over
`C` in the THM-3696 collision ring.  Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`; primes mean `d/db`.  The W004
fibre word is

```text
00; 01+10; 02+11; 12+20; 03+21; 13+22; 23.            (2)
```

We close the two placements

```text
X: wt(P)=(-n-2,-2,2n-2),
   wt(Q)=(1-n,1,n+1,3n+1),
   scalar arm 11=(-2,1);

Y: wt(P)=(1-n,1,2n+1),
   wt(Q)=(-n-2,-2,n-2,3n-2),
   scalar arm 11=(1,-2).                       n>=1.   (3)
```

## 1. Family X: a high Kummer sheet

Take `n>=2` in Family X and set

```text
delta=gcd(n+2,n-1)=gcd(n-1,3),
alpha=(n+2)/delta,              beta=(n-1)/delta,

epsilon=gcd(2n-2,3n+1)=gcd(n+3,8),
r=(2n-2)/epsilon,               s=(3n+1)/epsilon,
k=(n+3)/epsilon,
m=(n+5)/delta,                  ell=3/delta.            (4)
```

The singleton rows give nonzero constants `a,c,d,t` and nonconstant
`H,K in C[b]` with

```text
f0=aH^alpha,          g0=cH^beta,
f2=dK^r,              g3=tK^s.                         (5)
```

THM-3696 membership gives `h|H` and `b|K`.  Put `L=g1`, `M=f1`, and

```text
E=epsilon H'K+delta HK',
rho=ats/(dr)=at(3n+1)/[d(2n-2)].                       (6)
```

The upper row `03+21` has the complete solution

```text
L=rho H^alpha K^k+J,                 J^epsilon=lambda K.
                                                               (7)
```

Indeed the first term is a direct particular solution.  Subtracting it
leaves `K'J-epsilon KJ'=0`, equivalently `(J^epsilon/K)'=0` in `C(b)`.
No vanishing of this high Kummer sheet is assumed.

The lowest row `01+10` is

```text
delta H M'-2MH'=
 D0 H^mK^(k-1)E+D1 H^ell(H'J+delta HJ'),              (8)

D0=a rho alpha k/(c beta),             D1=a alpha/(c beta).
```

The operator on the left sends `H^mK^k` to `kH^mK^(k-1)E`, because
`delta*m-2=n+3=epsilon*k`.  It sends `H^ell J` to
`H^ell(H'J+delta HJ')`, because `delta*ell-2=1`.  Therefore

```text
M=A0H^mK^k+B0H^ell J+F,
A0=a rho alpha/(c beta),        B0=a alpha/(c beta),
F^delta=kappa H^2.                                    (9)
```

## 2. Family X dies at both arms

Here `delta` is `1` or `3`.  For `delta=1`, one has `m=n+5>=7` and
`ell=3`; also `(9)` gives `F=kappa H^2`.  Thus every term of `M` is
divisible by `h^2`.

For `delta=3`, regularity of the weight `-(n+2)=-3alpha` coefficient gives

```text
h^ceil(3alpha/2) | H^alpha,
hence                                      h^2|H.     (10)
```

Now `m=(n+5)/3>=3` and `ell=1`.  If `kappa=0`, every term in `(9)` is again
`h^2`-divisible.  If `kappa!=0`, unique factorization gives
`H=v^3`, `F=Bv^2`; equation `(10)` gives `h|v`, so the same conclusion
holds.  Hence, uniformly,

```text
h^2|M.                                                 (11)
```

Address `11=(-2,1)` therefore has zero value at both Danielewski arms.
The other scalar address is `02=(-n-2,n+1)`, whose complementary magnitude
is at least four and hence also vanishes there by THM-3696.  The scalar row
cannot equal one.

## 3. Family Y at every generic scale

Now take `n>=3` in Family Y.  Set

```text
delta=gcd(n-1,n+2)=gcd(n-1,3),
alpha=(n-1)/delta,              beta=(n+2)/delta,

epsilon=gcd(2n+1,3n-2)=gcd(n-3,7),
r=(2n+1)/epsilon,               s=(3n-2)/epsilon,
k=(n-3)/epsilon,
rho=ats/(dr)=at(3n-2)/[d(2n+1)].                       (12)
```

The endpoint coefficients again have the form `(5)`.  The upper row now has

```text
L=rho H^alpha K^k+J,                 J^epsilon K^2=lambda.
                                                               (13)
```

The homogeneous equation is `2K'J+epsilon KJ'=0`.  Since the nonconstant
polynomial `K` has a root, evaluation of `(13)` there forces `lambda=0`,
and the domain then gives `J=0`.

If `delta=1`, then `alpha=n-1>=2`, so `h^2|H^alpha`.  If `delta=3`, the
same regularity argument as `(10)` gives `h^2|H`.  Thus

```text
h^2|L.                                                 (14)
```

The persistent address `11=(1,-2)` dies at both arms.  At `n=3`, the other
address `02` is also an arm address `(-2,1)`, but its negative coefficient
is `f0=aH^2`, so it dies as well.  For `n>=4`, address `02` has complementary
magnitude at least three and vanishes by THM-3696.  This closes every
`n>=3`.

## 4. Family Y at n=2: an exact local-order gap

At `n=2`, the endpoint coefficients are

```text
f0=aH,       g0=cH^4,       f2=dK^5,       g3=tK^4.   (15)
```

The upper row integrates exactly to

```text
K^2L=rho HK+lambda,                  rho=4at/(5d).     (16)
```

Evaluation at a root of `K` gives `lambda=0`.  Hence, for a polynomial `G`,

```text
H=KG,                         L=rho G.                 (17)
```

Because `L` has weight `-2`, membership gives `h|G`.  At either arm
`zeta in {+-1}`, address `02=(-1,0)` contributes zero.  If the scalar row
were one, address `11=(1,-2)` would force

```text
ord_zeta(G)=1,                 ord_zeta(M)=0.          (18)
```

After `(17)`, the lowest row `01+10`, up to multiplication by `-1`, is

```text
a rho G(2K'G+KG')
 +4cK^3G^3[KG M'+M(K'G+KG')]=0.                       (19)
```

Put `q=ord_zeta(K)>=0`.  Under `(18)`, the first summand in `(19)` has
exact order `q+1`: the leading coefficient in its parenthesis is
`2q+1!=0`.  The second has exact order `4q+3`: the leading coefficient in
its bracket is `q+1!=0`.  These orders can never agree.  Thus `(19)` is
impossible, closing `n=2`.

## 5. The n=1 boundaries

For Family X at `n=1`, singleton `23` has weights `(0,4)` and makes `f2` a
scalar.  For Family Y, singleton `00` has weights `(0,-3)` and makes `f0` a
scalar.  Deleting the scalar leaves at most two graded pieces in the first
output and four in the second, impossible by THM-3583.

## 6. Scope and reproduction

Sections 1--5 close both named W004 scalar-`02+11` families in `(3)` at
every scale.  Other W004 placements, W005--W006, arbitrary `3 x 4` supports,
unrestricted quartic C3/cofactor data, and `JC(2)` remain open.

Run

```bash
python3 -B 04-computation/jacobian_y0_w004_scalar02_persistent_arm_pair_thm3735.py
python3 -B -O 04-computation/jacobian_y0_w004_scalar02_persistent_arm_pair_thm3735.py
```

Both commands must agree byte for byte with the frozen transcript.  Audit
promotion is required before this proof candidate enters proved canon.
