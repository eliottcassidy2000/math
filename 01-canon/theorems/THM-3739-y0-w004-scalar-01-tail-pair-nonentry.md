---
id: THM-3739
title: "W004 scalar-01 tail-pair nonentry"
status: >
  PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.  Both W004
  all-scale tail placements with scalar fibre 01+10 are Darboux-empty in the
  y=0 collision ring.  Even scales die by the inherited singleton-parity
  gate.  At every odd scale, two zero rows turn the remaining coefficient
  into a Mobius first integral and then force the cube of a nonconstant
  affine polynomial to be a monomial of incompatible degree.  This is not
  proved canon until audit promotion and does not close all exceptional W004
  placements, W005--W006, general quartic C3 data, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion reconstructs the W004
  word, both absolute placements, the n=1 and even-scale inherited gates,
  both endpoint gcd factorizations, all four transport signs, both pairs of
  reduced differential rows in forty exact odd-scale controls, both first
  integrals, and the final degree gaps.  Normal and optimized runs byte-match
  the frozen transcript.  Independent hostile audit remains open.
depends_on:
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3728-y0-w004-scalar-12-charge-normalization-nonentry
  - THM-3733-y0-w004-scalar-12-persistent-arm-pair-nonentry
  - THM-3735-y0-w004-scalar-02-persistent-arm-pair-nonentry
script: 04-computation/jacobian_y0_w004_scalar01_tail_pair_thm3739.py
output: 05-knowledge/results/jacobian_y0_w004_scalar01_tail_pair_thm3739.out
script_sha256: df7c9c42482104aad1e3066d5855a4a06ea6004e185653d745bf7ce568402d28
output_sha256: f46a9c93c312f1e018663c14bb249d66454851560f2e89bea58b5d7ca7bbd433
hash_basis: LF-normalized bytes
---

# THM-3739 -- the lowest-scalar W004 tail pair is empty

**PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.**  Work over
`C` in the THM-3696 collision ring.  Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`; primes mean `d/db`.  The W004 ray
has fibre word

```text
00; 01+10; 02+11; 12+20; 03+21; 13+22; 23.            (2)
```

We close the two opposite placements on its lowest double:

```text
A: wt(P)=(-2,n-2,3n-2),
   wt(Q)=(1-n,1,n+1,3n+1),       arm 01=(-2,1);

B: wt(P)=(1-n,1,2n+1),
   wt(Q)=(-2,n-2,2n-2,4n-2),    arm 10=(1,-2).         (3)
```

At `n=1`, singleton `00` has weights `(-2,0)` in A and `(0,-2)` in B,
contradicting the THM-3606 singleton gate.  At every even `n`, the two
negative endpoint weights have gcd `gcd(2,n-1)=1`; the negative arm
coefficient therefore has common-base exponent two.  THM-3613 rejects both
placements, including `n=2`.  It remains to treat odd `n>=3`.

Put

```text
m=(n-1)/2.                                             (4)
```

THM-3696 supplies `h|H` and `b|K`, so `H,K`, and `HK^2` are nonconstant.

## 1. Orientation A: two rows force an impossible cube

The endpoint singleton rows have gcds

```text
gcd(2,n-1)=2,             gcd(3n-2,3n+1)=1.           (5)
```

Hence, for nonzero constants `a,c,d,t`,

```text
f0=aH,       g0=cH^m,       f2=dK^(3n-2),
g3=tK^(3n+1).                                           (6)
```

Write `R=f1,L=g1,M=g2` and

```text
mu=t(3n+1)/[d(3n-2)].                                  (7)
```

The rows `03+21` and `13+22` integrate completely to

```text
L=lambda K+a mu H K^3,
M=kappa K^(n+1)+mu R K^3.                              (8)
```

For example, after subtracting the displayed particular solution, the first
row is `K'L-KL'=0`, and the second is the derivative of
`M/K^(n+1)`.  Thus `lambda,kappa in C` are the exhaustive constants, not an
ansatz restriction.

In the fraction field put

```text
U=HK^2,                 phi=R/K^(n-2),
A=lambda+3a mu U,
Psi=(n+1)kappa+3mu phi,
C_A=cd(3n-2)m.                                           (9)
```

Direct substitution into the next two zero rows gives, with no omitted
factor,

```text
02+11 = K^(n-1)[A phi'+a Psi U']=0,                    (10)

12+20 = K^(2n-1)[Psi phi'-C_A U^(m-1)U']=0.           (11)
```

The first equation has an elementary Mobius first integral.  Since
`A'=3a mu U'`,

```text
(A phi)'=-a(n+1)kappa U',
A phi=eta-a(n+1)kappa U,
A Psi=N_A:=(n+1)kappa lambda+3mu eta.                  (12)
```

Use `(10)` in `(11)` and cancel the nonzero polynomial `U'`.  Equations
`(10)--(12)` imply the polynomial identity

```text
C_A A^3 U^(m-1)=-a N_A^2.                             (13)
```

The left side has degree `m+2>=3` as a polynomial in the nonconstant `U`,
because `A` has nonzero linear coefficient `3a mu`.  The right side is a
constant.  This is impossible, including when `N_A=0`.  Notice that the
merged scale `n=3` is included: then `(13)` says that the cube of a
nonconstant affine polynomial is constant.

## 2. Orientation B: the dual degree mismatch

Now the endpoint gcds are

```text
gcd(n-1,2)=2,             gcd(2n+1,4n-2)=1,           (14)
```

and we may write

```text
f0=aH^m,       g0=cH,       f2=dK^(2n+1),
g3=tK^(4n-2).                                           (15)
```

Again put `R=f1,L=g1,M=g2`, and now set

```text
mu=t(4n-2)/[d(2n+1)].                                  (16)
```

The same two upper rows integrate completely to

```text
L=lambda K^(n-2)+a mu H^m K^(2n-3),
M=kappa K^(2n-2)+mu R K^(2n-3).                       (17)
```

Define in `C(b)`

```text
V=HK^2,                 U=H^mK^(n-1)=V^m,
phi=R/K,
A=(n-2)lambda+a mu(2n-3)U,
Psi=2(n-1)kappa+mu(2n-3)phi,
C_B=cd(2n+1).                                           (18)
```

This time the two lower zero rows reduce exactly to

```text
02+11 = K^(n-1)[A phi'+a Psi U']=0,                   (19)

12+20 = K^(2n-1)
  [Psi phi'-C_B U'/(mV^(m-1))]=0.                     (20)
```

Since `A'=a mu(2n-3)U'`, equation `(19)` integrates to

```text
A phi=eta-2a(n-1)kappa U,
A Psi=N_B:=2(n-1)(n-2)kappa lambda+mu(2n-3)eta.       (21)
```

Combining `(19)--(21)` and cancelling `U'` yields

```text
C_B A^3=-am N_B^2 V^(m-1).                            (22)
```

If `N_B=0`, this contradicts nonzero `C_B,A`.  Otherwise the left side has
degree `3m` in the nonconstant `V`, while the right side has degree `m-1`.
Since `m>=1`, equality is impossible.  This independently includes `n=3`,
where A and B are the same actual ordered support and `(22)` compares a
nonconstant cube with a constant.

## 3. Exact frontier after the closure

The THM-3613 stabilized W004 ledger has twelve placements.  Its parity gate
removes the two even-parity versions of `(3)` and three other tail schemes.
Among the surviving schemes, THM-3722, THM-3724, THM-3727, THM-3733, and
THM-3735 close every family except `(3)` and the separately reserved
THM-3728 anchor-`20` family.  Thus, after this theorem is audited, THM-3728 is
the only uncaptured stabilized-tail placement.  This is a tail census, not a
claim that every exceptional scale or all of W004 is already closed.

Other exceptional W004 placements, W005--W006, arbitrary `3 x 4` supports,
unrestricted quartic C3/cofactor data, and `JC(2)` remain open.

## 4. Reproduction

Run

```bash
python3 -B 04-computation/jacobian_y0_w004_scalar01_tail_pair_thm3739.py
python3 -B -O 04-computation/jacobian_y0_w004_scalar01_tail_pair_thm3739.py
```

Both commands must agree byte for byte with the frozen transcript.  Audit
promotion is required before this proof candidate enters proved canon.
