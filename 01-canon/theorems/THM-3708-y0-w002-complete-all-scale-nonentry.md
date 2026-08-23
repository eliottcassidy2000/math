---
id: THM-3708
title: "Complete all-scale W002 nonentry in the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT.  The full W002 support ray
  (n,n;n,n,2n) is Darboux-empty in the y=0 collision ring at every positive
  integral scale.  An exact actual-support census has 0,3,4,6 placements at
  n=1,2,3,n>=4.  THM-3704 closes the two lowest-double orientations.  Two
  arm-square gates, two generic Euler factorizations, two exceptional K|H
  peels, and one nonscalar-to-scalar bracket transfer close every residual
  placement.  This closes the whole named word W002, not the full 3x4 cell or
  JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion independently enumerates
  actual supports rather than anchor labels, checks the n=2,3 mergers and all
  stable families, and verifies every generic residue branch, exceptional
  identity, sign, Euler factor, and triple-to-scalar transfer.  Normal and
  optimized runs byte-match the frozen transcript.  Independent derivations
  checked the census, every C/D integration, both E transfer branches, and the
  nonunit conclusions.  The audit caught and repaired the false intermediate
  inference logged as MISTAKE-443; no theorem conclusion changed.
depends_on:
  - THM-3603-three-by-four-additive-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-3704-y0-w002-lowest-double-two-orientation-all-scale-nonentry
related:
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
script: 04-computation/jacobian_y0_w002_complete_all_scale_thm3708.py
output: 05-knowledge/results/jacobian_y0_w002_complete_all_scale_thm3708.out
script_sha256: a138674db44eaebc2b273a329ab313397a691231dc099248b894a9594abcd31d
output_sha256: b230ba6e8d00676d81a09f1387a85ed0dd292fcb62f2450f3aab6394be05df1e
hash_basis: LF-normalized bytes
---

# THM-3708 -- the entire W002 ray is empty

**PROVED + VERIFIED-EXACT.**  This theorem completes one of the seven
sumset-seven rays in the first `3 x 4` support frontier.  The proof retains
actual supports: two eligible arm addresses on one support are not counted as
two different coefficient systems.

Work over `C` in the THM-3696 collision ring

```text
h=1-b^2,
R=C[e,ce,bc] subset D=C[b,c,e]/(c^2e-h),
wt(b,c,e)=(0,1,-2).                                    (1)
```

For homogeneous coefficient polynomials use

```text
W_(r,s)(F,G)=sF'G-rFG'.                                (2)
```

## 1. Exact placement census

The W002 ray has gaps and fibres

```text
(X,Y;U,V,W)=(n,n;n,n,2n),                    n>=1,

00; 01+10; 02+11+20; 12+21; 03+22; 13; 23.             (3)
```

Anchor every nonsingleton fibre at `(-2,1)` or `(1,-2)` and impose the
same-sign-or-both-zero gate on the three singleton cells `00,13,23`.  After
deduplicating equal actual supports, the exhaustive list is:

| family | scalar fibre | `wt(P)` | `wt(Q)` | live scales |
|---|---|---|---|---|
| A | `01+10` | `(-2,n-2,2n-2)` | `(1-n,1,n+1,3n+1)` | `n>=3` |
| B | `01+10` | `(1-n,1,n+1)` | `(-2,n-2,2n-2,4n-2)` | `n>=2` |
| C | `02+11+20` | `(-2,n-2,2n-2)` | `(1-2n,1-n,1,2n+1)` | `n>=3` |
| D | `02+11+20` | `(1-n,1,n+1)` | `(-n-2,-2,n-2,3n-2)` | `n>=2` |
| E | `12+21` | `(1-n,1,n+1)` | `(-2n-2,-n-2,-2,2n-2)` | `n>=2` |
| F | `03+22` | `(-2,n-2,2n-2)` | `(1-4n,1-3n,1-2n,1)` | `n>=3` |

Thus the numbers of actual placements are

```text
n=1: 0,              n=2: 3,              n=3: 4,
n>=4: 6.                                                  (4)
```

At `n=3`, A=B and C=D as actual support pairs, although each merged pair has
two eligible arm addresses.  At `n=2`, only B,D,E occur.  The companion
reconstructs the fibre and anchor equations and verifies this census exactly;
the affine signs stabilize after the displayed thresholds.

THM-3704 closes A and B at every live scale.  It remains to close C--F.

## 2. The two arm-square families

### 2.1 Family F

On singleton `00`, the negative weights are `(-2,1-4n)`.  Their gcd is one,
so common-power rigidity gives

```text
f_0=aH^2,                         g_0=b_0H^(4n-1).      (5)
```

The address `03=(-2,1)` is the unique arm-capable address in the scalar
fibre.  Membership gives `h|f_0`, hence squarefreeness gives `h|H`; therefore
`f_0` has arm order at least two.  Its scalar bracket vanishes at `b=+-1`,
while every other scalar address vanishes there by THM-3696.  Family F is
empty.

### 2.2 Family C away from its merger

Here singleton `00` has weights `(-2,1-2n)`, again of gcd one, and

```text
f_0=aH^2,                         g_0=b_0H^(2n-1).      (6)
```

For `n>=4`, `02=(-2,1)` is the unique arm-capable scalar address, so the same
order-two argument closes C.

At `n=3`, the second scalar address `11=(1,-2)` is also arm-capable.  The
positive singleton component gives

```text
f_1=c_0K,                 f_2=d_0K^4,                 g_3=t_0K^7. (7)
```

Writing `g_2=M,g_1=L`, the rows `03+22=0` and `12+21=0` integrate to

```text
M=kappa K+gamma H^2K^3,       gamma=7at_0/(4d_0),
L=-beta H^2,                  beta=c_0gamma/(4d_0).    (8)
```

The scalar triple then factors exactly as

```text
1=(H'K+HK') [
  2aH(kappa+3gamma H^2K^2)+2c_0beta H-20b_0d_0H^4K^3
  ].                                                         (9)
```

Section 6 proves that the displayed Euler factor is a nonunit.  This closes
the exceptional merger and hence all of C.

## 3. Family D: a generic Euler factor and an exceptional divisibility peel

First take `n>=3` and put

```text
delta=gcd(n-1,3),       alpha=(n-1)/delta,
omega=(n+2)/delta.                                      (10)
```

The singleton components give, with nonzero constants,

```text
f_0=aH^alpha,                   g_0=b_0H^omega,
f_1=c_0K,                       f_2=d_0K^(n+1),
g_3=t_0K^(3n-2).                                      (11)
```

For `g_2=M,g_1=L`, the top two collision rows integrate completely to

```text
M=kappa K^(n-2)+gamma H^alpha K^(2n-3),
gamma=a t_0(3n-2)/(d_0(n+1)),

L=-beta H^alpha K^(n-3),
beta=c_0gamma/(d_0(n+1)).                              (12)
```

There is no omitted polynomial homogeneous solution in the second line: it
would be a scalar multiple of `K^-2`, while `K` is nonconstant.  Define

```text
E_delta(H,K)=H'K+delta HK'.                            (13)
```

The scalar row is

```text
1=E_delta(H,K) [
 a alpha H^(alpha-1){
    kappa(n-2)K^(n-3)
   +gamma(2n-3)H^alpha K^(2n-4)
   +(c_0beta/a)K^(n-3)}
 -b_0d_0(n+1)omega H^(omega-1)K^n
 ].                                                    (14)
```

At `n=2`, the generic display would have a negative `K` exponent, and the
polynomiality condition contains real information.  The top row gives

```text
M=kappa+gamma HK,               gamma=4at_0/(3d_0).   (15)
```

The next row integrates to

```text
K^2L=rho-beta HK,               beta=c_0gamma/(3d_0). (16)
```

Since `K` is nonconstant, evaluate `(16)` at a root of `K`: `rho=0`.  Hence
`K|H`; write `H=KJ`.  Then `L=-beta J`, and the scalar row becomes

```text
1=(J'K+2JK') [
 a gamma K^2J+c_0beta-12b_0d_0K^6J^3
 ].                                                    (17)
```

This closes D at every scale.

## 4. Family E: transfer the unresolved bracket into the scalar row

Take `n>=3` and put

```text
delta=gcd(n-1,4),       alpha=(n-1)/delta,
epsilon=2(n+1)/delta,   q=(n+3)/delta.                 (18)
```

The integer `q` exists because `delta|4`.  Singleton rigidity gives

```text
f_0=aH^alpha,                   g_0=b_0H^epsilon,
f_1=c_0K,                       f_2=d_0K^(n+1),
g_3=t_0K^(2n-2).                                      (19)
```

Let `g_2=M,g_1=L`.  The upper double integrates to

```text
M=gamma H^alpha K^(n-3),
gamma=a t_0(2n-2)/(d_0(n+1)).                         (20)
```

The lower double has the complete form

```text
L=L_0+eta H^qK,                   eta=b_0c_0epsilon/(a alpha),
W_(1-n,-n-2)(H^alpha,L_0)=0.                           (21)
```

No classification of the possible common-power root in `L_0` is needed.
Put

```text
U=W_(1,-n-2)(K,L_0).                                  (22)
```

Using

```text
(H^alpha)'K+(n-1)H^alpha K'
       =alpha H^(alpha-1)E_delta(H,K),                 (23)
```

direct expansion of the remaining triple row and scalar row gives

```text
c_0U+E_delta Psi=0,
d_0(n+1)K^nU+E_delta Phi=1,                            (24)
```

for explicit polynomials `Psi,Phi`.  The point is that the *same* unresolved
bracket `U` occurs in both equations: the identity

```text
W_(n+1,-n-2)(K^(n+1),L_0)
       =(n+1)K^n W_(1,-n-2)(K,L_0)                    (25)
```

is exact.  Eliminate `U` from `(24)`:

```text
c_0=E_delta[c_0Phi-d_0(n+1)K^nPsi].                   (26)
```

This is impossible because `c_0` is a nonzero scalar and `E_delta` is a
nonunit.

At `n=2`, top-row integration first gives

```text
K^2M=rho+gamma HK,             gamma=2at_0/(3d_0).    (27)
```

As before, a root of `K` forces `rho=0` and `K|H`.  Put `H=KJ`; then
`M=gamma J`.  The lower row has

```text
L=L_0+eta K^6J^5,              eta=6b_0c_0/a,
W_(-1,-4)(KJ,L_0)=0.                                  (28)
```

With `U=W_(1,-4)(K,L_0)`, the triple and scalar rows are

```text
c_0U+(J'K+2JK')Psi=0,
3d_0K^2U+(J'K+2JK')Phi=1.                             (29)
```

The same elimination closes the active `n=2` boundary.

## 5. Exact placement exhaustion

Sections 2--4 close C,D,E,F.  THM-3704 closes A,B.  The census in Section 1
is exhaustive, including the `n=2` boundary and the `n=3` mergers.  Therefore

```text
no retained W002 3 x 4 Darboux pair exists in R at any n>=1.              (30)
```

## 6. Every displayed Euler factor is a nonunit

In every generic family, a negative coefficient is a positive power of `H`.
THM-3696 gives `h|H^r`, hence squarefreeness gives `h|H`.  A positive
weight-one coefficient is a nonzero multiple of `K`, so THM-3696 gives
`b|K`.  Consequently

```text
deg(H'K+delta HK')=deg H+deg K-1>=2,                  (31)
```

with nonzero leading multiplier `deg H+delta deg K`.

In the exceptional cases `H=KJ`.  Here `h|H` gives `deg H>=2`, while the
positive weight-one coefficient gives `b|K`, so `K` is nonconstant.  Hence

```text
deg(J'K+2JK')=deg J+deg K-1=deg H-1>=1,              (32)
```

with nonzero leading multiplier `deg J+2deg K`.  This argument deliberately
does not divide `h|H` through `K`: the positive module allows `K` itself to
contain `h`.  None of the factors in `(9),(14),(17),(26),(29)` is a unit.

## 7. Frontier and reproduction

This proves complete nonentry for the named oriented word W002 and its output
transpose inside the y=0 collision ring.  It does not close the other
sumset-seven words, arbitrary `3 x 4` supports, general quartic C3 data, or
`JC(2)`.

Run

```bash
python3 -B 04-computation/jacobian_y0_w002_complete_all_scale_thm3708.py
python3 -B -O 04-computation/jacobian_y0_w002_complete_all_scale_thm3708.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
