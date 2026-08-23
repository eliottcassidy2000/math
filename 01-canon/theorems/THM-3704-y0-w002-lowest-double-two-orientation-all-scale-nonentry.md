---
id: THM-3704
title: "W002 lowest-double two-orientation all-scale nonentry in the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the W002
  support ray with gaps
  (n,n;n,n,2n), the complete placement family whose scalar fibre is the
  lowest double 01=10 is empty in the y=0 collision ring, in both arm
  orientations and at every positive scale where the placement survives the
  inherited singleton gate.  The singleton endpoints give two UFD bases;
  the top two collision rows integrate exactly; and the scalar row then has
  the common nonunit Euler-Wronskian factor H'K+dHK', where
  d=gcd(2,n-1).  This closes one named two-orientation placement family, not
  the whole word W002, the whole 3x4 cell, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  PASS.  The independent hostile audit rechecked both absolute placements,
  every small-scale singleton exit, both negative common-power exponent
  systems, the connected positive gcd-one components, and independently
  derived all four rational-derivative integrations, including the active
  n=2 branch.  It also checked module divisibility, Euler degree, and the
  exact two-orientation scope.  Normal and optimized companion runs
  byte-match the stored transcript.
depends_on:
  - THM-3603-three-by-four-additive-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
script: 04-computation/jacobian_y0_w002_lowest_double_all_scale_thm3704.py
output: 05-knowledge/results/jacobian_y0_w002_lowest_double_all_scale_thm3704.out
script_sha256: 8d66b3915b0bc76d0a359d958cf2883a1a8906f8bdf3460b65e73f1c0492e00f
output_sha256: 8b9eea84575f237edd791dde54de99b47117930e5be557f30494f12a7a7e6179
hash_basis: LF-normalized bytes
---

# THM-3704 -- the lowest scalar double on W002 has an Euler factor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is a
collision-ring-specific all-scale closure inside the first `3 x 4` frontier.
It uses only three of the four multi-address equations: the unused row is not
silently assumed.

All rings are over `C`.  Use the THM-3696 presentation

```text
h=1-b^2,
R=C[e,ce,bc] subset D=C[b,c,e]/(c^2e-h),
wt(b,c,e)=(0,1,-2).                                    (1)
```

For homogeneous pieces put

```text
W_(r,s)(F,G)=sF'G-rFG',
{c^rF,c^sG}=c^(r+s+1)W_(r,s)(F,G).                    (2)
```

## 1. The W002 ray and its exact boundary

The THM-3603 word `W002` is the ray

```text
(X,Y;U,V,W)=(n,n;n,n,2n),                 n>=1,        (3)
```

with fibres

```text
00; 01+10; 02+11+20; 12+21; 03+22; 13; 23.            (4)
```

We close the complete family in which the lowest double `01+10` is the
scalar fibre and either of its two addresses is given an arm orientation.
There are two absolute weight placements:

```text
orientation A, (P_0,Q_1)=(-2,1):
  wt(P)=(-2,n-2,2n-2),
  wt(Q)=(1-n,1,n+1,3n+1);                              (5)

orientation B, (P_1,Q_0)=(1,-2):
  wt(P)=(1-n,1,n+1),
  wt(Q)=(-2,n-2,2n-2,4n-2).                            (6)
```

At `n=1`, the singleton `00` has weights `(-2,0)` in A and `(0,-2)`
in B, so both placements fail the inherited zero/nonzero singleton gate.  At
`n=2`, A fails on singleton `13`, whose weights are `(0,7)`.  Orientation B
does survive at `n=2` and is treated explicitly in Section 4.  Thus the live
ranges are

```text
A: n>=3,                         B: n>=2.              (7)
```

## 2. Orientation A

Fix `n>=3` and put

```text
delta=gcd(2,n-1),       epsilon=2/delta,
m=(n-1)/delta.                                            (8)
```

The singleton row `00` and the common-power lemma give, for nonzero constants
`a,b_0`,

```text
f_0=a H^epsilon,                    g_0=b_0 H^m.        (9)
```

The positive singleton rows `13` and `23` form one component.  Its three
weights have gcd

```text
gcd(n-2,2n-2,3n+1)=gcd(n-2,2,7)=1.                  (10)
```

Hence, for nonzero constants `c_0,d_0,t_0`,

```text
f_1=c_0 K^(n-2),
f_2=d_0 K^(2n-2),
g_3=t_0 K^(3n+1).                                     (11)
```

Write `g_1=L,g_2=M` and set `F=H^epsilon`.  The top row `03+22=0` is

```text
W_(-2,3n+1)(aF,t_0K^(3n+1))
 +W_(2n-2,n+1)(d_0K^(2n-2),M)=0.                     (12)
```

It integrates exactly to

```text
M=kappa K^(n+1)+gamma F K^(n+3),
gamma=a t_0(3n+1)/(d_0(2n-2)).                        (13)
```

Substitute this into `12+21=0`.  Its complete polynomial solution is

```text
L=lambda K-beta F K^3,
beta=c_0(n-2)gamma/(d_0(2n-2)).                       (14)
```

Indeed, the two ODEs reduce respectively to derivatives of `FK^2` and
`FK^3`; the constants `kappa,lambda` are the complete integration constants.

Now use the scalar row `01+10=1`.  Define

```text
E_delta(H,K)=H'K+delta HK'.                            (15)
```

Because `epsilon delta=2` and `m delta=n-1`, direct substitution of
`(9),(11),(14)` gives the exact factorization

```text
1=E_delta(H,K) [
    a epsilon H^(epsilon-1)(lambda-3beta H^epsilon K^2)
   -b_0 c_0(n-2)m H^(m-1)K^(n-3)
  ].                                                   (16)
```

This is already impossible; the unused middle triple row `02+11+20=0`
cannot repair a factor of the scalar equation.

## 3. Orientation B for `n>=3`

Keep `(8)`.  The negative singleton and the positive singleton component now
give

```text
f_0=a H^m,                         g_0=b_0 H^epsilon,
f_1=c_0K,                          f_2=d_0K^(n+1),
g_3=t_0K^(4n-2).                                      (17)
```

With `F=H^m`, the same top-down integration of `03+22=0` and `12+21=0`
gives, with the signs shown,

```text
M=kappa K^(2n-2)+gamma F K^(3n-3),
gamma=a t_0(4n-2)/(d_0(n+1)),

L=lambda K^(n-2)-beta F K^(2n-3),
beta=c_0 gamma/(d_0(n+1)).                            (18)
```

The scalar row factors as

```text
1=E_delta(H,K) [
    a m H^(m-1){lambda(n-2)K^(n-3)
                 -beta(2n-3)H^mK^(2n-4)}
   -b_0c_0 epsilon H^(epsilon-1)
  ].                                                   (19)
```

Again the scalar unit has acquired the same Euler-Wronskian factor.

## 4. The active `n=2` boundary in orientation B

At `n=2`, one has `delta=1,m=1,epsilon=2`.  Equations `(17)--(18)` remain
polynomial and specialize to

```text
f_0=aH,                g_0=b_0H^2,
f_1=c_0K,              f_2=d_0K^3,             g_3=t_0K^6,

M=kappa K^2+gamma HK^3,
L=lambda-beta HK.                                      (20)
```

The two collision rows still vanish exactly.  The scalar equation is

```text
1=-H(H'K+HK')(a beta+2b_0c_0),                        (21)
```

which is impossible for the same reason as the general factorization.

## 5. Why the Euler factor is never a unit

In A, the actual weight `-2` coefficient is `aH^epsilon`; in B it is
`b_0H^epsilon`.  THM-3696 gives `h|H^epsilon`.  Since `h=1-b^2` is
squarefree,

```text
h|H.                                                    (22)
```

In A, `K^(n-2)` is a positive-weight coefficient with `n-2>=1`; in B,
`K` itself has weight one.  The positive module line of THM-3696 therefore
gives

```text
b|K.                                                    (23)
```

Thus `deg H>=2` and `deg K>=1`.  The leading coefficient of `(15)` is

```text
(deg H+delta deg K)lc(H)lc(K),                         (24)
```

which is nonzero.  Hence

```text
deg E_delta=deg H+deg K-1>=2.                          (25)
```

It cannot divide the unit in `(16),(19)`, or `(21)`.

## 6. Exact scope

The theorem closes both arm orientations of the scalar placement on the
lowest W002 double `01+10`, at every scale where that retained placement
exists.  Output exchange is automatic.  Other scalar fibres and placements
on W002 remain open unless closed elsewhere; this is not a closure of the
whole word.  The proof is all-degree in the coefficient polynomials and
all-scale in `n`; the companion is an identity audit, not a bounded search.

## Reproduction

```bash
python3 -B 04-computation/jacobian_y0_w002_lowest_double_all_scale_thm3704.py
python3 -B -O 04-computation/jacobian_y0_w002_lowest_double_all_scale_thm3704.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
