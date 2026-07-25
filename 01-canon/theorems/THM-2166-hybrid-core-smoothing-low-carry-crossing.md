---
id: THM-2166
title: "Hybrid whole-core smoothing and low-carry crossing"
status: >
  PROVED + VERIFIED-EXACT. For a zero-measure defect-six split with seven
  retained speeds in {1,...,13}, a common scalar frequency has nonzero
  magnitude at most 708. Its six far coefficients have height 298 and every
  nonzero one is a 7-unit. The core side has an arithmetic representation
  on at most two retained speeds with coefficient height 57. Height 57 is
  certificate-minimal for representing the full possible carry interval.
  Combining any nonzero far coefficient with THM-2169's relation on that
  deletion gives a fixed-base-7 rank-two carrier with a 7-unit anchored
  minor and at most 28442762 joint carry pairs.
  This is not Fourier-vector sparsity: the natural core product has nonzero
  support-three coefficients. Two independent exact implementations verify
  the 1,716-core geometry, Jackson ledger, and sparse carry bank.
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2145
  - THM-2163
  - THM-2169
related:
  - THM-2054
  - THM-2162
  - THM-2167
  - THM-2171
script: 04-computation/lrc14_hybrid_core_low_carry_referee_codex_20260724.py
output: 05-knowledge/results/lrc14_hybrid_core_low_carry_referee_codex_20260724.out
script_sha256: 025c9d2a516587c12bba120fce06bb6579869bc0bf37577b2d507f336703c4a8
output_sha256: 044c687196e37851f81ccbee6695cc01b780d7fef022e348326649726e9c1e30
independent_script: 04-computation/lrc14_hybrid_whole_core_smoothing_thm2166.py
independent_output: 05-knowledge/results/lrc14_hybrid_whole_core_smoothing_thm2166.out
independent_script_sha256: 246d8775ac7e8332b45bfd227d77c958aad697480964328fc07d4c598c75be34
independent_output_sha256: e101b90e7fc045e1e05c44d5765d3df67b265d97e5c56edf8b7a89e453bf896d
hash_basis: working-tree bytes (LF)
---

# THM-2166 -- hybrid whole-core smoothing and low-carry crossing

At radius `1/14`, put

```text
G={x in R/Z: ||x||>=1/14},
G_A={t: at in G for every a in A}.                    (1)
```

Let `E` be a seven-subset of `{1,...,13}` and let `F` be a disjoint
six-set of distinct positive integer speeds. If

```text
mu(G_E intersect G_F)=0,                              (2)
```

then there are an integer `nu` and coefficient vectors `(a_f)_(f in F)`,
`(b_e)_(e in E)` such that

```text
0<|nu|<=708,
sum_(f in F) a_f f=nu,
|a_f|<=298,             a_f!=0 implies 7 does not divide a_f,
sum_(e in E) b_e e=-nu,
#{e:b_e!=0}<=2,         |b_e|<=57.                   (3)
```

Thus `sum a_f f+sum b_e e=0` is genuinely crossing. We call `nu` its
**cut carry**.

## 1. The hybrid smoothing inequality

Use THM-2145's normalized Jackson kernel and safe-interval smoother

```text
J_N=F_N^2/integral F_N^2,
q_N=J_N*1_G,
eta_N=2 integral ||x||J_N(x)dx.                      (4)
```

For the far block take `N=150` and set

```text
P_F(t)=product_(f in F)q_150(ft).                     (5)
```

The factor height is `2N-2=298`. The product telescope, Haar invariance,
and the six-speed floor used in THM-2145 give

```text
||P_F-1_(G_F)||_1<=6eta_150,
integral P_F>=61/273-6eta_150.                        (6)
```

For the core, first assemble its exact one-dimensional safe set. Let
`beta_E=mu(G_E)` and let `K_E` be the number of its positive-length circular
components; isolated weak-safe points are deliberately omitted. Put

```text
R_E=J_355*1_(G_E).                                    (7)
```

Then

```text
0<=R_E<=1,        integral R_E=beta_E,
Fourier(R_E,k)=0 for |k|>708.                         (8)
```

A translation by `x` changes a union of `K_E` circular intervals by at most
`2K_E||x||` in symmetric-difference measure. Hence

```text
||R_E-1_(G_E)||_1<=K_E eta_355.                       (9)
```

This is the BV-component distinction of THM-2162: an isolated point matters
for a weak witness but cancels from every `L^1` variation.

For `[0,1]`-valued `P,R` and indicators `f,g`,

```text
PR-fg=(P-f)R+f(R-g),
|PR-fg|<=|P-f|+|R-g|.                                 (10)
```

Apply (10) with `(P,R,f,g)=(P_F,R_E,1_(G_F),1_(G_E))`. By (2),

```text
integral P_F R_E<=6eta_150+K_E eta_355.               (11)
```

Therefore a nonzero common Fourier frequency follows if

```text
(61/273-6eta_150)beta_E
   >6eta_150+K_E eta_355.                             (12)
```

## 2. Exact all-core certificate

The danger-comb arrangement for speeds `1,...,13` has `178` rational
boundary points and `177` open cells. The primary exact sweep computes
`(beta_E,K_E)` on every one of the `binomial(13,7)=1716` cores.

The Jackson odd-mode formula, using `pi<355/113`, gives

```text
eta_150<439/156250,
eta_355<371/312500.                                   (13)
```

Substituting these caps into the left side of (12) minus the right, its
unique minimum is

```text
41050267/1222741406250>0                              (14)
```

at

```text
E_*=(1,5,7,8,9,11,13),
beta_(E_*)=45107/229320,       K_(E_*)=20.            (15)
```

The same core has four isolated weak-safe points, correctly charged zero.
The rational-`pi` ledger is negative with `J_354`, so `355` is the first of
these two adjacent Jackson choices to close this certificate. No optimality
among other kernels is claimed.

The independent companion repeats every load-bearing finite step by a
different path:

- iterative exact interval intersection versus global arrangement-cell runs;
- direct integer convolution of Fejer coefficients versus the cubic Jackson
  formula; and
- ordinary integer sets versus the primary bit masks for sparse carries.

The two core paths agree on every mass and positive-component count. The
independent component distribution is

```text
K_E : 12  14  16  17  18  20  22  24  26  28
count:  4  31 171   1 452 584 262 152  54   5.        (16)
```

Normal and optimized executions reproduce both stored transcripts.

## 3. Frequency extraction

Equations (6)--(15) imply

```text
(integral P_F)(integral R_E)>integral P_F R_E.         (17)
```

Finite Fourier orthogonality says the zero-mode contribution on the right is
the left side of (17). Thus some nonzero integer `nu` satisfies

```text
Fourier(P_F,nu) Fourier(R_E,-nu)!=0,
0<|nu|<=708.                                          (18)
```

Expanding the nonzero coefficient of `P_F` supplies integers `a_f` with

```text
sum_f a_f f=nu.                                       (19)
```

By THM-2145, the nonzero Fourier support of `q_150` consists exactly of
the `k` with

```text
0<|k|<=298,             7 does not divide k.          (20)
```

Zero modes are allowed in the product, so (19)--(20) give the far part of
(3).

## 4. Sparse arithmetic encoding of the core

For height `B`, define

```text
C_B(E)=union_({e,f} subset E)
       {ce+df: |c|,|d|<=B}.                           (21)
```

An exhaustive exact computation over all `1,716` cores and all `1,417`
possible carries proves

```text
[-708,708] subset C_57(E)          for every E.        (22)
```

It also proves sharpness for this unrestricted carry interval. At height
`56` exactly three cores fail:

```text
(1,2,3,4,5,6,7)       misses +/-699,+/-705,+/-706,
(1,2,3,4,5,6,8)       misses +/-701,
(1,2,4,5,6,8,10)      misses +/-701.                  (23)
```

Equations (18), (22) supply the core coefficients in (3).

There is also a conceptual, nonsharp proof of support two. Every seven-subset
`E` either contains `1`, or, after partitioning `{2,...,13}` into

```text
(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),
```

pigeonhole forces a consecutive pair `(r,r+1)`. One may represent `-nu` by
`-nu*1`, or by `nu*r-nu*(r+1)`. This has height `708`; (22) is the exact
height refinement.

## 5. Hostile type check: support two is not Fourier sparsity

The whole-core scalar smoother forgets labelled coefficient-vector origins.
On the natural seven-torus,

```text
Phi(x_1,...,x_7)=product_(j=1)^7 1_G(x_j)
```

has support-three coefficient

```text
Fourier(Phi,(1,1,1,0,0,0,0))
=Fourier(1_G,1)^3 Fourier(1_G,0)^4
=(-sin(pi/7)/pi)^3(6/7)^4!=0.                        (24)
```

Thus higher-support core Fourier coefficients do not vanish. The valid
mechanism is

```text
scalar cutoff |nu|<=708
  + a finite additive property of seven-subsets of {1,...,13}
  -> a support-two arithmetic re-encoding.            (25)
```

This distinction is load-bearing for any generalization.

## 6. Defect-six and radix consequences

If exactly one far coefficient in (19) is nonzero, then

```text
|a_f|f=|nu|<=708,
```

so that far speed is at most `708`. Otherwise at least two far speeds have
a height-`298`, `7`-unit near-cancellation of nonzero magnitude at most
`708`.

The full coefficient vector `m=(a,b)` satisfies

```text
||m||_1<=6*298+2*57=1902.                             (26)
```

THM-2163 gives the generic radix bound `|kappa_j|<1902`. For a sharper
dyadic tail, divide only the far speeds:

```text
F=2^j Z_j+R_j,       D_j=Z_j mod 2,
lambda_j=(a.R_j-nu)/2^j=-a.Z_j.                       (27)
```

Then

```text
lambda_0=-nu,
lambda_(j+1)=(lambda_j+a.D_j)/2.                      (28)
```

For `j>=4` all core speeds have terminated, so `lambda_j` is the full
relation carry and only six far owners remain. Moreover

```text
|lambda_j|
 <=((2^j-1)||a||_1+|nu|)/2^j
 <=1788-1080/2^j<1788,                               (29)
```

hence `|lambda_j|<=1787`.

This is a much smaller state than THM-2145's raw cut-carry range, but it does
not bound radix depth. The owner mask remains essential, the multi-far branch
still permits arbitrarily large speeds, and defect six and LRC(14) remain
open.

## 7. A fixed-base-7 rank-two lift

The `7`-unit support in (20) does more than control one relation. Since
`nu!=0`, choose a far label `f` with

```text
a_f!=0,                       7 does not divide a_f. (30)
```

Apply THM-2169 to the deletion of this label. After primitive normalization
it gives a second relation

```text
0!=u in Lambda(E union F),
u_f=0,                  ||u||_infinity<=1247.        (31)
```

Primitivity makes some `u_j` a `7`-unit. Necessarily `j!=f`, and the
anchored minor is

```text
det ((m_f,m_j),(u_f,u_j))=a_f u_j!=0 mod 7.          (32)
```

Thus the two relations have rank two over `F_7`: the defect-six branch has
a **fixed** base seven, rather than THM-2167's row-adaptive prime. If `m`
is itself primitive-normalized, (32) persists because its coefficient gcd
divides the `7`-unit `a_f`.

The state bounds are explicit. Equation (26) and positivity give at most

```text
1902-1=1901
```

carries for `m`. Since `u_f=0` and the other twelve coefficients form a
primitive vector of height `1247`,

```text
||u||_1<=11*1247+1246=14963,
```

so `u` has at most `14962` carries. Therefore

```text
joint carry pairs <=1901*14962=28442762,              (33)
unrestricted digit-fibre size =7^11=1977326743.       (34)
```

After sorting there are at most fourteen owner suffixes, giving
`398198668` carry-owner states over the whole path. With THM-2171's
quotient-tie sidecar, the effective ordered algebraic path cap is

```text
26*28442762=739511812.                                (35)
```

From digit level two onward every core speed has terminated because
`E subset {1,...,13}` and `7^2=49`; only six far owners remain. The tail
owner/tie pair changes at most six owner labels and six quotient cuts, so
it assumes at most thirteen values. The corresponding ordered tail cap is

```text
13*28442762=369755906.                                (36)
```

While the two pivot owners survive, the owner-restricted digit fibre has
exactly `7^(|O|-2)` elements.

There is also a useful pump boundary. A repeated augmented THM-2171 state
at levels `2<=j<k` preserves the literal core `E`: all of its digits lie
below the deleted block. It preserves both relations and hence the affine
far cut `sum_f a_f f=nu`; active far speeds remain at least `49`, so the
core/far split stays disjoint, and the tie mask preserves strict order.
Primitivity is automatic without normalization: no integer `d>=2` has
seven multiples in `{1,...,13}`, so the unchanged seven-set `E` has gcd
one. The only target loss is therefore lonely-phase geometry. A parity
sidecar repairs literal residues modulo `14`, but no fixed residue depth
controls the unbounded phase target of THM-2174. Live rank can also fall
when the anchored owners terminate.

This is a structural specialization, not a smaller automaton than the
adaptive THM-2167 carrier and not a depth bound.

QED.
