---
id: THM-3733
title: "W004 scalar-12 persistent-arm pair nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Both named W004
  all-scale placements with scalar fibre 12+20 and a persistent arm at
  address 12 are Darboux-empty in the y=0 collision ring.  In the (-2,1)
  orientation, endpoint transport forces second-order arm vanishing; in the
  (1,-2) orientation it forces an uncancellable negative bivariate charge.
  This does not close the distinct anchor-20 family, all of W004, the full
  3x4 cell, general quartic C3 data, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion checks both absolute
  placements, all endpoint gcd residues modulo 5, 4, and 3, both primitive
  and lowest-row signs across complete residue windows, the strengthened
  delta=5 regularity order, the arm-death and negative-Laurent alternatives,
  and both n=1 deletions.  An independent derivation reconstructed both
  primitives and lowest operators, the delta=1/5 valuation split, the exact
  arm-law use, the rational-function normalization and reduced numerator,
  and both 2x4 boundaries.  Normal and optimized runs byte-match the frozen
  transcript.
depends_on:
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3728-y0-w004-scalar-12-charge-normalization-nonentry
script: 04-computation/jacobian_y0_w004_scalar12_persistent_arm_pair_thm3733.py
output: 05-knowledge/results/jacobian_y0_w004_scalar12_persistent_arm_pair_thm3733.out
script_sha256: bf3ad842e40730bbdffb38178eacfa9404a683d107af014fed2dbf14c557a83f
output_sha256: e0abe2c70417427ab91b18df049659169260bc4dc21813a0817b3abd199a046a
hash_basis: LF-normalized bytes
---

# THM-3733 -- both persistent-arm W004 scalar-12 families are empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over `C`
in the THM-3696 collision ring.  Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`; primes mean `d/db`.  The W004
supports `(0,n,3n)` and `(0,n,2n,4n)` have fibres

```text
00; 01+10; 02+11; 12+20; 03+21; 13+22; 23.            (2)
```

We close the following two actual-support placements:

```text
A: wt(P)=(-n-2,-2,2n-2),
   wt(Q)=(1-2n,1-n,1,2n+1),
   scalar arm 12=(-2,1);

B: wt(P)=(1-n,1,2n+1),
   wt(Q)=(-2n-2,-n-2,-2,2n-2),
   scalar arm 12=(1,-2).                       n>=1.   (3)
```

## 1. Family A endpoint transport

First take `n>=2` in Family A.  Set

```text
delta=gcd(n+2,2n-1)=gcd(n+2,5),
alpha=(n+2)/delta,              beta=(2n-1)/delta,

epsilon=gcd(2n-2,2n+1)=gcd(n-1,3),
r=(2n-2)/epsilon,               s=(2n+1)/epsilon,
k=3/epsilon,                    m=5/delta.              (4)
```

The singleton rows `00` and `23` give nonzero constants `a,c,d,t` and
nonconstant `H,K in C[b]` with

```text
f0=aH^alpha,          g0=cH^beta,
f2=dK^r,              g3=tK^s.                         (5)
```

THM-3696 membership gives `h|H` and `b|K`.  Put `L=g1`, `M=f1`, and

```text
E=epsilon H'K+delta HK'.                               (6)
```

For

```text
gamma=(n-1)/epsilon,             q=(n+2)/epsilon,
rho=ats/(dr)=at(2n+1)/[d(2n-2)],                       (7)
```

the upper row `03+21` is exactly

```text
(K^gamma L)'=rho(H^alpha K^q)'.                        (8)
```

At a root of `K`, the integration constant vanishes, and `q-gamma=k`.
Therefore

```text
L=rho H^alpha K^k.                                    (9)
```

The lowest row `01+10` now becomes

```text
delta H M'-2MH'=D H^m K^(k-1)E,
D=a rho alpha k/(c beta).                             (10)
```

Since `delta*m-2=3=epsilon*k`, the displayed differential operator sends
`H^mK^k` to `kH^mK^(k-1)E`.  Hence, with

```text
A0=a rho alpha/(c beta),              F=M-A0H^mK^k,   (11)
```

equation `(10)` is equivalent to

```text
delta H F'-2FH'=0,
hence                         F^delta=kappa H^2        (12)
```

for some `kappa in C`.

## 2. Family A dies at both Danielewski arms

Here `delta` is either `1` or `5`.  If `delta=1`, then `m=5` and `(12)`
simply says `F=kappa H^2`; since `h|H`, both terms in `(11)` are divisible
by `h^2`.

Suppose `delta=5`.  Then `n+2=5alpha`, and regularity of the weight
`-(n+2)` coefficient gives

```text
h^ceil(5alpha/2) | H^alpha.                            (13)
```

Because `ceil(5alpha/2)>2alpha`, equation `(13)` forces `h^3|H`.  If
`kappa=0`, this already makes `(11)` divisible by `h^2`.  If `kappa!=0`,
unique factorization in `(12)` gives, after scalar rescaling,

```text
H=v^5,                         F=Bv^2.                 (14)
```

Now `h|v`, so again `h^2|F` and `h^2|M`.  Thus in every case

```text
h^2 | M.                                               (15)
```

At `b=+-1`, THM-3696 says that only complementary weights `(-2,1)` (or
the reversed orientation) can contribute to a scalar bracket.  Family A's
address `12=(-2,1)` has negative coefficient `M`; `(15)` makes its arm value
zero.  The other scalar address is

```text
20=(2n-2,1-2n),                                       (16)
```

which has complementary magnitude at least three and therefore also
vanishes at both arms.  The scalar row consequently evaluates to zero at
`b=+-1`, contradicting the required value one.

## 3. Family B has a negative bivariate charge

Now take `n>=2` in Family B and set

```text
delta=gcd(n-1,2n+2)=gcd(n-1,4),
alpha=(n-1)/delta,              beta=(2n+2)/delta,

epsilon=gcd(2n+1,2n-2)=gcd(n-1,3),
r=(2n+1)/epsilon,               s=(2n-2)/epsilon,
k=3/epsilon,                    m=-4/delta.             (17)
```

Again the endpoints have the form `(5)`, with these new exponents.  Put
`L=g1`, `M=f1`, and retain `E` from `(6)`.  This time

```text
gamma=(n+2)/epsilon,             q=(n-1)/epsilon,
rho=ats/(dr)=at(2n-2)/[d(2n+1)].                       (18)
```

The upper and lowest rows are

```text
(K^gamma L)'=rho(H^alpha K^q)',
L=rho H^alpha K^(-k),                                 (19)

delta H M'+MH'=-D H^mK^(-k-1)E,
D=a rho alpha k/(c beta).                             (20)
```

Since `delta*m+1=-3=-epsilon*k`, the operator in `(20)` sends
`H^mK^(-k)` to `-kH^mK^(-k-1)E`.  Therefore

```text
A0=a rho alpha/(c beta),              F=M-A0H^mK^(-k),
delta H F'+FH'=0,
hence                         F^delta H=kappa.         (21)
```

If `kappa=0`, then `M=A0/[H^(4/delta)K^k]`, which is not a polynomial.
If `kappa!=0`, prime valuations in `(21)` force, after scalar rescaling,

```text
H=v^delta,                       F=B/v,          B!=0. (22)
```

Consequently

```text
M=[A0+Bv^3K^k]/[v^4K^k].                              (23)
```

Every irreducible factor of the nonconstant denominator `vK` sees the
numerator in `(23)` reduce to the nonzero scalar `A0`.  The fraction is
therefore reduced and nonpolynomial.  This contradiction closes Family B
before its scalar row.

## 4. The n=1 boundaries

In Family A at `n=1`, singleton `23` has weights `(0,3)` and forces the
weight-zero piece `f2` to be a scalar.  In Family B, singleton `00` has
weights `(0,-4)` and forces `f0` to be a scalar.  Deleting the scalar in
either case leaves at most two graded pieces in the first output and four in
the second, impossible by THM-3583.

## 5. Scope and reproduction

Sections 1--4 close both named W004 scalar-`12+20` families in `(3)` at
every scale.  The distinct anchor-`20` family reserved at THM-3728 is not a
dependency and is not claimed here.  Other W004 placements, W005--W006,
arbitrary `3 x 4` supports, unrestricted quartic C3/cofactor data, and
`JC(2)` remain open.

Run

```bash
python3 -B 04-computation/jacobian_y0_w004_scalar12_persistent_arm_pair_thm3733.py
python3 -B -O 04-computation/jacobian_y0_w004_scalar12_persistent_arm_pair_thm3733.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
