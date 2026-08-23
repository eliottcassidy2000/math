---
id: THM-3715
title: "Odd-scale W003 lowest-scalar two-orientation half-charge nonentry"
status: >
  PROVED + VERIFIED-EXACT.  Every post-parity W003 placement whose scalar
  fibre is the lowest double 01+10 is Darboux-empty in the y=0 collision ring.
  Such placements occur only at odd n>=3: the two orientations merge into one
  actual support at n=3 and are distinct thereafter.  The upper rows impose a
  half-charge law on every hidden monomial, forcing the entire scalar row to
  contain H'K+2HK'.  This closes the complete named scalar-fibre pair, not all
  of W003, the full 3x4 cell, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  SELF-AUDITED.  The exact companion reconstructs the post-parity
  actual-support census, preserves the n=3 two-arm merger, verifies both
  all-odd-scale upper-row integrations, and checks every monomial charge and
  the full scalar factorizations.  Normal and optimized runs byte-match the
  frozen transcript.  Independent hostile audit remains open.
depends_on:
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3714-y0-w003-scalar-fibre-02-two-orientation-euler-nonentry
script: 04-computation/jacobian_y0_w003_odd_lowest_scalar_thm3715.py
output: 05-knowledge/results/jacobian_y0_w003_odd_lowest_scalar_thm3715.out
script_sha256: dba9fdc35256d78fefa5d81db9fa1f6f71e26e1e3ce25e098bae1762330f4999
output_sha256: 33c6df7b6bb4561ba80b0c416e2b1d000f82b3b5d8815806336e4d4875eafcc3
hash_basis: LF-normalized bytes
---

# THM-3715 -- the odd W003 lowest scalar fibre is empty

**PROVED + VERIFIED-EXACT.**  Work over `C` in the THM-3696 collision ring.
Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

## 1. Exact odd placement scope

The W003 fibres are

```text
00; 01+10; 02+11; 03+12+20; 13+21; 22; 23.            (2)
```

On scalar fibre `01+10`, the inherited gates and THM-3613 parity obstruction
leave no actual support at `n=1` or at even `n`.  For odd `n>=3`, the two
orientations are

| family | `wt(P)` | `wt(Q)` |
|---|---|---|
| O1 | `(1-n,1,2n+1)` | `(-2,n-2,2n-2,3n-2)` |
| O2 | `(-2,n-2,3n-2)` | `(1-n,1,n+1,2n+1)` |

At `n=3` they are the same actual support pair

```text
P=(-2,1,7),                       Q=(-2,1,4,7),         (3)
```

with both arm addresses retained.  At every odd `n>=5` they are distinct.

Write

```text
n=2m+1,                               m>=1.            (4)
```

## 2. Orientation O1

The singleton components give

```text
f_0=aH^m,                         g_0=b_0H,
f_2=d_0K^(2n+1),                 g_2=e_0K^(2n-2),
g_3=t_0K^(3n-2).                                      (5)
```

The positive gcd is one: a common divisor of the first two weights can only
divide three, while `3n-2` removes that possibility.  Let `L=f_1` and
`M=g_1`.  The zero triple `03+12+20` and upper double `13+21` integrate to

```text
L=kappa K-alpha H^mK^n+beta HK^3,
M=lambda K^(n-2)+rho LK^(n-3),                        (6)

alpha=(3n-2)a t_0/[2(n-1)e_0],
beta=(2n+1)b_0d_0/[2(n-1)e_0],
rho=(3n-2)t_0/[(2n+1)d_0].                            (7)
```

Every monomial `H^pK^q` in `L` satisfies

```text
q=2p+1,                                                (8)
```

and every monomial in `M` satisfies

```text
q=n-2+2p=2m-1+2p.                                     (9)
```

Define

```text
E_2(H,K)=H'K+2HK'.                                    (10)
```

For a term of `L`, `(8)` gives

```text
W_(1,-2)(H^pK^q,H)=-(2p+1)H^pK^(q-1)E_2,             (11)
```

up to its scalar coefficient.  For a term of `M`, `(9)` gives

```text
W_(1-n,n-2)(H^m,H^pK^q)
       =m q H^(m+p-1)K^(q-1)E_2,                     (12)
```

again up to its scalar coefficient.  Therefore the full scalar row is

```text
W_(1-n,n-2)(f_0,M)+W_(1,-2)(L,g_0)=E_2 Phi_1         (13)
```

with `Phi_1` polynomial.

## 3. Orientation O2

The singleton components now give

```text
f_0=aH,                           g_0=b_0H^m,
f_2=d_0K^(3n-2),                 g_2=e_0K^(n+1),
g_3=t_0K^(2n+1).                                      (14)
```

Here the positive gcd is one already from `gcd(n+1,2n+1)=1`.  Put

```text
U=kappa-alpha HK^2+beta H^mK^(n-1),                  (15)

alpha=(2n+1)a t_0/[(n+1)e_0],
beta=(3n-2)b_0d_0/[(n+1)e_0].                         (16)
```

The same two upper rows integrate completely to

```text
f_1=K^(n-2)U,
g_1=K(lambda+rho U),
rho=(2n+1)t_0/[(3n-2)d_0].                            (17)
```

Every monomial of `f_1` satisfies `(9)`, and every monomial of `g_1`
satisfies `(8)`.  The two termwise identities dual to `(11),(12)` give

```text
W_(-2,1)(f_0,g_1)+W_(n-2,1-n)(f_1,g_0)=E_2 Phi_2.    (18)
```

Thus O2 also has a nonunit factor in its scalar row.

## 4. The factor is a nonunit and the merger is retained

Negative coefficient membership and squarefreeness give `h|H`.  The
positive common-power component gives `b|K`.  Hence

```text
deg E_2=deg H+deg K-1>=2,                             (19)
```

with nonzero leading multiplier `deg H+2deg K`.  Equations `(13),(18)`
cannot equal one.

At `n=3`, `(3)` has two arm addresses rather than two coefficient systems.
The scalar fibre consists exactly of the two bracket terms in `(13)` or
`(18)`, and both termwise charge identities hold.  Thus the full merged
scalar row, not merely one chosen anchor, is divisible by `E_2`.

## 5. Sharpened W003 frontier and scope

This theorem closes the complete post-parity lowest-scalar pair.  Together
with THM-3695 and THM-3714, only the following two W003 placement families
remain:

```text
C: scalar 03+12+20,
   P=(-n-2,-2,2n-2), Q=(1-2n,1-n,1,n+1),       n>=2;

D: scalar 13+21,
   P=(1-3n,1-2n,1), Q=(-n-2,-2,n-2,2n-2),      n>=3.  (20)
```

This exact residual statement does not assert that C or D is nonempty as a
Darboux system.  W003, W004--W006, arbitrary `3 x 4` supports, general
quartic C3 data, and `JC(2)` remain open.

Run

```bash
python3 -B 04-computation/jacobian_y0_w003_odd_lowest_scalar_thm3715.py
python3 -B -O 04-computation/jacobian_y0_w003_odd_lowest_scalar_thm3715.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
