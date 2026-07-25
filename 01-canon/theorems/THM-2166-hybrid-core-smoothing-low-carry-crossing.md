---
id: THM-2166
title: "Hybrid whole-core smoothing and low-carry crossing"
status: >
  PROVED + VERIFIED-EXACT. In every zero-measure defect-six 6+7 split, with
  the retained seven-core contained in {1,...,13}, there is a nonzero scalar
  crossing frequency of size at most 708. The six far coefficients have
  height at most 298 and every nonzero one is prime to 7. The core side can
  be arithmetically re-encoded on at most two retained speeds, with coefficient
  height at most 57, which is certificate-minimal for representing the full
  interval of possible carries. This support-two conclusion is not Fourier-vector
  sparsity: the natural core product has nonzero support-three coefficients.
  The exact 1,716-core mass/component sweep and the Jackson ledger have two
  independent implementations.
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2145
  - THM-2162
  - THM-2163
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
G={x in R/Z: ||x||_(R/Z)>=1/14},
G_A={t: at in G for every a in A}.                    (1)
```

Let

```text
E subset {1,...,13},       |E|=7,                     (2)
F={six distinct positive integer speeds},             (3)
E intersect F=empty.
```

Suppose

```text
mu(G_(E union F))=0.                                  (4)
```

Then there are an integer `nu`, far coefficients `(a_f)_(f in F)`, and
core coefficients `(b_e)_(e in E)` such that

```text
0<|nu|<=708,                                          (5)
sum_(f in F) a_f f=nu,                                (6)
|a_f|<=298, and a_f!=0 implies 7 does not divide a_f, (7)
sum_(e in E) b_e e=-nu,                               (8)
#{e:b_e!=0}<=2,             |b_e|<=57.                (9)
```

In particular

```text
sum_(f in F) a_f f+sum_(e in E)b_e e=0               (10)
```

is genuinely crossing: both labelled restrictions have nonzero value. The
number `nu` is the **cut carry**.

## 1. Smooth the far block factorwise and the core as one set

Use the normalized Jackson kernels and interval smoothers of THM-2145:

```text
J_N=F_N^2/integral F_N^2,       q_N=J_N*1_G,          (11)
eta_N=2 integral ||x|| J_N(x) dx.                     (12)
```

For the far block use `N=150`, so `q=q_150` has degree `298`,

```text
0<=q<=1,       ||q-1_G||_1<=eta_150,
eta_150<439/156250.                                   (13)
```

Put

```text
Q_F(t)=product_(f in F) q(ft).                        (14)
```

The product telescope and Haar invariance give

```text
||Q_F-1_(G_F)||_1<=6 eta_150.                         (15)
```

By the six-speed floor used in THM-2145,

```text
integral Q_F>=mu(G_F)-6eta_150
             >=61/273-6eta_150.                      (16)
```

Now retain the actual one-dimensional core geometry. Let `K_E` be the
number of **positive-length circular components** of `G_E`; isolated weak-safe
points are not counted. Define

```text
R_E=J_355*1_(G_E).                                    (17)
```

Then

```text
0<=R_E<=1,             integral R_E=mu(G_E),
Fourier(R_E,nu)=0 for |nu|>708.                       (18)
```

If a circle set is a union of `K` positive-length intervals, translation by
`x` changes its indicator in `L^1` by at most `2K||x||`. Averaging against
`J_355` therefore gives

```text
||R_E-1_(G_E)||_1<=K_E eta_355.                       (19)
```

This is the whole-core counterpart of THM-2162's BV component split. An
isolated point has measure zero and its two endpoint terms cancel; charging
it as a component would use the wrong carrier.

## 2. The hybrid crossing inequality

For any `[0,1]`-valued `Q,R` and indicators `f,g`,

```text
QR-fg=(Q-f)R+f(R-g),
|QR-fg|<=|Q-f|+|R-g|.                                 (20)
```

Apply this with

```text
(Q,R,f,g)=(Q_F,R_E,1_(G_F),1_(G_E)).
```

Equation (4) says `integral fg=0`, so (15), (19), and (20) imply

```text
integral Q_F R_E<=6eta_150+K_E eta_355.               (21)
```

Consequently a strict common-frequency certificate follows whenever

```text
(61/273-6eta_150)mu(G_E)
   >6eta_150+K_E eta_355.                             (22)
```

This is where smoothing the core as a whole helps: its mean is exactly
`mu(G_E)`. Only its one-dimensional BV boundary is charged on the right.

## 3. Exact 1,716-core and Jackson ledger

The companion computes `(mu(G_E),K_E)` for every

```text
binomial(13,7)=1716
```

core in two independent ways:

1. iterative intersection of exact closed rational safe intervals, deleting
   isolated intervals only when counting `K_E`; and
2. a sweep over the `178` exact global boundary points and `177` open cells,
   where mass is the sum of safe-cell widths and `K_E` is the cyclic number
   of safe-cell runs.

The two ledgers agree on every core. The positive-component distribution is

```text
K_E : number of cores
 12 :   4
 14 :  31
 16 : 171
 17 :   1
 18 : 452
 20 : 584
 22 : 262
 24 : 152
 26 :  54
 28 :   5.                                            (23)
```

For `N=355`, the exact Jackson coefficient formula and a second direct
integer convolution of the Fejer coefficients agree at every mode. With
`pi<355/113`, the odd-mode ledger gives

```text
C_0=29826035,
sum_(1<=k<=707, k odd) C_k/k^2>36709039,
eta_355<11872/10000000=371/312500.                    (24)
```

The same independent coefficient check at `N=150` reproduces (13).

Insert the rational caps from (13), (24) into the left side of (22) minus
the right side, for all `1,716` cores. Its exact minimum is

```text
41050267/1222741406250>0,                             (25)
```

uniquely at

```text
E=(1,5,7,8,9,11,13),
mu(G_E)=45107/229320,          K_E=20.                (26)
```

That core also has four isolated weak-safe points, which correctly contribute
zero to (19). Equations (22)--(26) prove

```text
(integral Q_F)(integral R_E)>integral Q_F R_E.         (27)
```

The margin in (25) is small enough that both the positive-component convention
and the exact per-core pairing of mass with component count are load-bearing;
combining the global mass minimum with the global component maximum would
lose the certificate.

The same exact rational-`pi` ledger is negative with `J_354` on the core.
Thus `N=355` is the first of these two adjacent Jackson choices to certify
this inequality; no optimality among other positive kernels is asserted.

The independent companion named in the header repeats every load-bearing
finite claim by a different path: iterative rational interval intersection
versus the arrangement run count, direct integer convolution of the Fejer
coefficients versus the cubic Jackson formula, and ordinary integer sets
versus the primary bit masks for the height-`57` carry bank. Normal and
optimized executions match its stored transcript.

## 4. Extract the scalar frequency and the far coefficients

Both factors in (27) are finite trigonometric polynomials. Fourier
orthogonality gives

```text
integral Q_F R_E
 =sum_nu Fourier(Q_F,nu) Fourier(R_E,-nu).             (28)
```

The zero-mode term on the right is the left side of (27). Hence at least one
nonzero integer `nu` has

```text
Fourier(Q_F,nu) Fourier(R_E,-nu)!=0.                  (29)
```

Equation (18) gives `0<|nu|<=708`.

Expand the nonzero coefficient of `Q_F`. At least one convolution summand
supplies integers `(a_f)` with

```text
sum_f a_f f=nu.                                       (30)
```

THM-2145 proves that the nonzero Fourier support of `q_150` consists exactly
of the integers `k` with

```text
0<|k|<=298,             7 does not divide k.          (31)
```

Zero modes are also allowed in the product. Thus (30)--(31) prove (6)--(7).
Notice that the core coefficient in (29) is used only to bound the scalar
frequency and make it nonzero; no coefficient-vector origin has yet been
assigned to it.

## 5. The support-two core encoding

Every seven-subset `E` of `{1,...,13}` obeys the following elementary
dichotomy:

```text
1 in E,
or E contains one whole pair among
(2,3),(4,5),(6,7),(8,9),(10,11),(12,13).              (32)
```

Indeed, if `1` is absent, seven elements selected from the six displayed
pairs force a full pair by pigeonhole.

If `1 in E`, set

```text
b_1=-nu,             b_e=0 otherwise.                 (33)
```

If `(r,r+1) subset E`, set

```text
b_r=nu,              b_(r+1)=-nu,                    (34)
```

and all other core coefficients to zero. Then

```text
nu*r-nu*(r+1)=-nu.                                   (35)
```

This already proves a support-two encoding of height `708`, but the exact
finite carrier is much smaller. For an integer `B`, define

```text
C_B(E)=union_({e,f} subset E)
       {ae+bf: |a|,|b|<=B}.                           (36)
```

An exhaustive integer-set computation over all `1,716` cores proves

```text
[-708,708] subset C_57(E)       for every E.           (37)
```

The bound is sharp for representing the entire unrestricted carry interval:
at height `56` exactly three cores fail,

```text
(1,2,3,4,5,6,7)       misses +/-699,+/-705,+/-706,
(1,2,3,4,5,6,8)       misses +/-701,
(1,2,4,5,6,8,10)      misses +/-701.                  (38)
```

Thus (5) and (37) supply coefficients satisfying (8)--(9), completing the
theorem. Equations (32)--(35) remain the conceptual reason support two is
always possible; (37) is the exact height refinement.

## 6. Hostile type check: this is not core Fourier sparsity

The strongest tempting interpretation of (9) is false. On the natural
seven-torus, put

```text
Phi(x_1,...,x_7)=product_(j=1)^7 1_G(x_j).            (39)
```

Its coefficient at the support-three vector

```text
(1,1,1,0,0,0,0)
```

is

```text
Fourier(1_G,1)^3 Fourier(1_G,0)^4
=(-sin(pi/7)/pi)^3(6/7)^4!=0.                        (40)
```

Thus core Fourier coefficients do **not** vanish outside support two. The
whole-core operation first restricts the labelled product to the
one-dimensional speed line and then smooths that scalar function. Different
coefficient vectors have already collided at a common scalar frequency, so
their owner support cannot be recovered.

The valid statement is precisely the repaired one proved above:

```text
scalar Fourier cutoff |nu|<=708
  + the special additive geometry of seven-subsets of {1,...,13}
  -> a new support-two arithmetic encoding.           (41)
```

This distinction is essential if the theorem is generalized to other core
sizes or banks.

## 7. Consequences and limits

Compared with THM-2145, the cut carry falls from `20860` to `708`, at the
same time as the core restriction contracts to two coefficients of height
`57`. This is a post-extraction re-encoding, not preservation of the
factorwise coefficient owners.

If exactly one far coefficient in (6) is nonzero, then

```text
|a_f|f=|nu|<=708,
```

so that far speed is itself at most `708`. Otherwise at least two far speeds
have a height-`298`, `7`-unit near-cancellation of nonzero magnitude at most
`708`.

The full relation satisfies

```text
||a,b||_1<=6*298+2*57=1902.                           (42)
```

Therefore THM-2163's generic radix carries obey `|kappa_j|<1902`. There is
also a sharper dyadic tail. Write

```text
F=2^j Z_j+R_j,             D_j=Z_j mod 2,
lambda_j=(a.R_j-nu)/2^j=-a.Z_j.                       (43)
```

Then

```text
lambda_0=-nu,
lambda_(j+1)=(lambda_j+a.D_j)/2.                      (44)
```

For `j>=4` all core speeds have terminated, so `lambda_j` is the full
relation carry and only the six far owners remain. Moreover

```text
|lambda_j|
 <=((2^j-1)||a||_1+|nu|)/2^j
 <=1788-1080/2^j<1788,                               (45)
```

hence the integer carry satisfies `|lambda_j|<=1787`. The owner-mask warning
of THM-2163 remains unchanged: neither finite bound controls the number of
radix levels.

The theorem does not prove that the extracted relation is independent of
previous relations, prevent a one-far/one-core reduced-pair tautology, bound
all six far speeds, or close defect six. What it proves is the sharper finite
object that a next descent may use:

```text
one nonzero scalar cut carry in [-708,708],
six height-298 far digits with nonzero digits prime to 7,
and at most two height-57 explicitly encoded core owners. (46)
```

QED.
