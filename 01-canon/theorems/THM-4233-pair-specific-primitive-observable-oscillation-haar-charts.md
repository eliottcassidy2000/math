---
id: THM-4233
title: "Pair-specific primitive-observable oscillation Haar charts"
status: >
  PROVED ANALYTIC REDUCTION + FINITE-EXACT + INDEPENDENTLY AUDITED FIXED-PAIR
  CERTIFICATE. For the displayed thirty-label pool, every nine-body is
  Haar-safe after adjoining every common dilation of the coprime pair
  (3713,5149). A second chart proves the same conclusion for every dilation
  of (5k+2,7k+3), exactly when k>=748 for the sufficient primitive-oscillation
  gate. Its period-seven primitive formula is proved for every k>=1. The
  adjacent coprime family has limiting primitive oscillation 36/2401 and is a
  rigorous hostile control. Arbitrary finite pair entry and LRC(14) remain
  OPEN.
source: codex-lrc14-pair-specific-oscillation-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
related:
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
primary_script: 04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.cpp
primary_output: 05-knowledge/results/lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.out
independent_audit_script: 04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_independent_audit_thm4233.cpp
independent_audit_output: 05-knowledge/results/lrc14_pair_specific_primitive_observable_gcd_one_chart_independent_audit_thm4233.out
primary_script_sha256: 9c02b42562d3b911e9df2a553844c6e1acb1ab6d64dee176d9d76c861f54dcc9
primary_output_sha256: 03e0c222ee7ba9e2384de76c33217c05887990c64ec24c7635ca59bc05b9df07
independent_audit_script_sha256: b6691337c81199e7459fbabfeec10e49d0899ee40368d7b1d729deabd7a0ad20
independent_audit_output_sha256: e43dc469c08d89a5e5ca20e0c734bf56f44d35aa84d5d4cac8ef304464487253
exact_tail_script: 04-computation/lrc14_resonance_zero_exact_tail_thm4233.py
exact_tail_output: 05-knowledge/results/lrc14_resonance_zero_exact_tail_thm4233.out
exact_tail_independent_audit_script: 04-computation/lrc14_resonance_zero_exact_tail_thm4233_independent_audit.py
exact_tail_independent_audit_output: 05-knowledge/results/lrc14_resonance_zero_exact_tail_thm4233_independent_audit.out
exact_tail_script_sha256: 721dbbf9d11a0814f7f8edfec43bc38d073abc5d915e93ad34af18b9bbbbe6b5
exact_tail_output_sha256: 56b6aa79076cf4ebbb51c12883323631f61615f72f5687068b7f553ba711b319
exact_tail_independent_audit_script_sha256: 27a7fb1b283b864480754c9ef3d839dcd3463a0de26b99442b21962e21ab9531
exact_tail_independent_audit_output_sha256: 9e6eeaa21bb98d2895664d616f7b39b181fd1f943b306caa9cc6771728f7ab2a
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary integer common-grid sweep and an independent
  reduced-rational safe-teeth intersection/integration audit reproduce the
  fixed-pair mass, primitive extrema, oscillation, threshold margins,
  collision ticks, and hostile controls. Their independently serialized
  semantic ledgers are 4725493ac81fe903 and 635148be446ddc28. O0/O3 replays
  match; the independent path also passes ASan/UBSan. For the exact tail, a
  second common-grid sweep checks k=1..1000 and far controls, while an ordered
  safe-teeth intersection independently checks k=1..850 and the same crossing;
  their ledgers are 388cbaa2d0d3242d and 83c13873d6a23452. The infinite
  residue law itself is supplied by the analytic seven-macrostate proof below.
---

# THM-4233 -- pair-specific primitive-observable oscillation Haar charts

**PROVED ANALYTIC REDUCTION + FINITE-EXACT + INDEPENDENTLY AUDITED FIXED-PAIR
CERTIFICATE; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Retain the thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive label set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

> **Fixed primitive-ratio ray.** For every integer `g>=1` and every
> `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {3713g,5149g}))>=4/63.                  (3)
> ```

In particular, `g=1` is a genuinely coprime comparable pair. Its coverage is
subsumed by the [THM-4231](THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift.md)
cofinite region because its maximum label exceeds `770`. Its gcd is below the
[THM-4228](THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray.md)
gate, and the pair is outside the
[THM-4227](THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge.md)
scale-separated wedge, so the exact primitive-observable certificate remains
a method sidecar not supplied by those earlier mechanisms.

There is also an analytic infinite chart. For `k>=1`, define

```text
u_k=5k+2,                 v_k=7k+3.                      (4)
```

> **Canceled-resonance ray family.** For every `k>=748`, every integer
> `g>=1`, and every `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {g u_k,g v_k}))>=4/63.                  (5)
> ```

The conclusion `(5)` is numerically subsumed by THM-4231, since at `k>=748`
already `u_k,v_k>1290>770`. Its new content is the cyclotomic-zero mechanism that makes
the primitive-observable oscillation `O(1/k)` without a common divisor.

For every body in `(3)` or `(5)`,
[THM-4150](THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer.md)
gives, for every positive integer `c` and every two distinct positive odd
integers `a,b`, an `x in R/Z` such that

```text
min_(w in 2c(B union {gu,gv}) union {a,b})||wx||>=1/14, (6)
```

where `(u,v)` is respectively `(3713,5149)` or `(u_k,v_k)`. These are
infinite thirteen-speed LRC(14) families, not a proof of LRC(14).

The closest proved mechanism is THM-4228's periodic-observable discrepancy
lemma and its exact depth-eight deck `E_3467`. The canonical hostile is the
adjacent coprime family `(n,n+1)`, whose primitive oscillation does not tend
to zero. The corrected near miss is to interpret failure of a sufficient
oscillation gate as literal Haar danger. The least-used sidecar is the
centered primitive of the **joint** two-comb observable, rather than either
marginal comb. The live concept board is

```text
depth-eight repair deck | pair density | primitive oscillation
endpoint collisions | slow-fast averaging | cyclotomic zeros.             (7)
```

## 2. Pair-specific transfer through the THM-4228 deck

Let

```text
j(t)=1_(||t||>=1/14),
A_(u,v)=G_u intersect G_v,
beta_(u,v)=mu(A_(u,v)),                                  (8)
```

where changing the values of `j` at its two endpoints is immaterial. Define
the one-periodic centered primitive

```text
H_(u,v)(t)=integral_0^t(1_(A_(u,v))(s)-beta_(u,v))ds,
omega_(u,v)=max H_(u,v)-min H_(u,v).                    (9)
```

If `U` is a union of `c_U` positive-length circular intervals, the oriented-
lift proof of THM-4228, equation `(24)`, gives for every integer `g>=1`

```text
mu(U intersect G_(gu) intersect G_(gv))
 >=beta_(u,v)mu(U)-c_U omega_(u,v)/g.                   (10)
```

There is no independence assumption in `(10)`: it is the identity
`G_(gu) intersect G_(gv)=m_g^(-1)(A_(u,v))`, followed component by component
by the endpoint cocycle of `H_(u,v)`.

Put

```text
beta_0=66/91,
T=(1650/8281)/3467=1650/28710227.                       (11)
```

For `R in binom(P,8)`, write

```text
U_R=G_(P\R),       M_R=mu(U_R),       c_R=#components(U_R). (12)
```

THM-4228 constructs a finite deck `E_3467` such that

```text
tau(E_3467)>9,
beta_0 M_R-c_R T>=4/63              for every R in E_3467. (13)
```

Consequently, any distinct positive pair `(u,v)` satisfying

```text
beta_(u,v)>=beta_0,             omega_(u,v)<=T           (14)
```

supplies every dilation `g>=1`. Indeed, for any `B in binom(P,9)`, `(13)`
provides `R in E_3467` disjoint from `B`. Equations `(10)`, `(13)`, and `(14)`
give

```text
mu(G_((P\R) union {gu,gv}))>=4/63.                     (15)
```

Since `B union {gu,gv}` is contained in the label set in `(15)`, safe-set
monotonicity proves the corresponding conclusion for `B`. More generally,
the same proof works when `omega_(u,v)/g<=T`; `(14)` is the convenient gate
that includes the primitive dilation `g=1`.

## 3. Exact coprime certificate for `(3713,5149)`

Set

```text
u=3713,                 v=5149,                 gcd(u,v)=1,
N=14uv=267655318.                                         (16)
```

Both `u` and `v` exceed `max(P)=290`. On the integer grid `y=t/N`, the
circular endpoint multisets of `G_u` and `G_v` are

```text
v(14i+1), v(14i-1) mod N,        0<=i<u,
u(14j+1), u(14j-1) mod N,        0<=j<v.                  (17)
```

There are `17724=2(u+v)` raw circular entries. Because
`u+v=8862=14(633)`, precisely two endpoint pairs coincide, at ticks

```text
95591185,                 172064133.                      (18)
```

These are null collisions, not positive-length cells. After adjoining the
linear sentinels `0,N`, the sweep has `17726` raw points, `17724` unique
points, and `17723` open cells. It finds

```text
safe cells=7595,          unsafe cells=10128,
safe components=7595,
S=safe ticks=196644720.                                  (19)
```

Thus

```text
beta_(u,v)=S/N=98322360/133827659,
beta_(u,v)-66/91=16387038/1739759567>0.                 (20)
```

For completeness, the exact primitive computation can be stated without any
floating-point convention. Let

```text
0=e_0<e_1<...<e_m=N
```

be the unique sorted linear points from `(17)` and the two sentinels, put
`d_i=e_(i+1)-e_i`, and let `a_i in {0,1}` be the joint-safe value on that
open cell. Starting with `C_0=0`, set

```text
C_(i+1)=C_i+d_i(Na_i-S).                                (21)
```

Then `H_(u,v)(e_i/N)=C_i/N^2`; since `H_(u,v)` is affine on each cell, its
extrema occur among these endpoint values. The exhaustive recurrence `(21)`
gives

```text
min C_i=-2057535171948 at t=207775767=55959N/72086,
max C_i= 2057535171948 at t= 59879551=16127N/72086,
N^2=71639369253681124,                                  (22)

min H=-138535899/4823550313337,
max H= 138535899/4823550313337,
omega_(u,v)=277071798/4823550313337.                    (23)
```

The decisive strict margin is

```text
T-omega_(u,v)
 =82934716896/2826229070241355051>0.                    (24)
```

Equations `(20)`, `(23)`, and `(24)` verify `(14)`. Section 2 proves `(3)`
for every `g>=1`, and THM-4150 proves `(6)`. This is the fixed-pair theorem.

The exact artifacts retain the hostile and positive controls

```text
(1,13):       beta=66/91, omega=990/8281,
(1,2):        beta=11/14, omega=11/98,
(2584,4181):  omega=42967355/553336008694>T.             (25)
```

The last line is a certificate-gate failure, not an unsafe-set witness. The
primary semantic ledger is `4725493ac81fe903`; the independent ledger is
`635148be446ddc28`. The primary path is the integer-grid recurrence
`(17)--(23)`. The audit path independently intersects reduced rational safe
teeth and integrates the resulting rational primitive. Neither path uses
floating point, imports the other's geometry, or treats the pair-selection
process as a proof of minimality.

## 4. The adjacent coprime family is a hostile limit

Let

```text
A_n=G_n intersect G_(n+1),
beta_n=mu(A_n),
H_n(t)=integral_0^t(1_(A_n)(s)-beta_n)ds.               (26)
```

Then

```text
lim_(n->infinity) osc(H_n)=36/2401.                     (27)
```

Here is a direct slow-fast proof. Define the autocorrelation

```text
C(z)=integral_(R/Z) j(x)j(x+z)dx.                       (28)
```

The danger arc has length `1/7`, so, with `||z||` the circle distance,

```text
C(z)=5/7+(1/7-||z||)_+,
integral_(R/Z)C(z)dz=36/49.                             (29)
```

On a cell `[a/n,(a+1)/n]`, substitute `t=(a+x)/n`. Apart from null walls,

```text
j(nt)j((n+1)t)=j(x)j(x+t).                              (30)
```

Replacing `t=(a+x)/n` in the second factor by `a/n` costs at most `4/n` in
the `x`-integral: if the two values differ, `x+a/n` lies within `1/n` of one
of the two walls of `j`. After the outer factor `1/n`, each complete cell
therefore costs at most `4/n^2`. Summing complete
cells, treating one final partial cell identically, and using the uniform
Riemann convergence of the continuous piecewise-linear function `C` proves

```text
integral_0^t j(ns)j((n+1)s)ds
 -> integral_0^t C(s)ds                                 (31)
```

uniformly for `0<=t<=1`. At `t=1`, `(31)` also gives `beta_n->36/49`.
Therefore `H_n` converges uniformly to

```text
K(t)=integral_0^t(C(s)-36/49)ds.                        (32)
```

On the first side of the triangular peak,
`C(s)-36/49=6/49-s`; its positive area on `[0,6/49]` is
`18/2401`. By reflection symmetry around zero, `K` has maximum `18/2401`,
minimum `-18/2401`, and oscillation `36/2401`. Oscillation is continuous
under uniform convergence, proving `(27)`.

In particular,

```text
36/2401-T=21012378/1406801123>0.                        (33)
```

Thus primitive comparable pairs do not acquire the tiny THM-4228 oscillation
gate merely by moving to large adjacent labels. This identifies why a pair-
specific cancellation, not size alone, is needed below the common-gcd or
cofinal-quadrant routes.

## 5. A cyclotomic-zero family with vanishing oscillation

### 5.1. Exact zero fast mean

For the pair `(4)`, note first that

```text
7u_k-5v_k=-1,                                           (34)
```

so every pair is coprime. Put

```text
F(x,y)=j(5x+2y)j(7x+3y)-36/49.                          (35)
```

The Fourier coefficients of `j` are

```text
c_0=6/7,
c_m=-sin(pi m/7)/(pi m)       for m!=0.                 (36)
```

In the `x`-average of the product in `(35)`, a Fourier pair `(a,b)` can
survive only when `5a+7b=0`, hence `(a,b)=(7h,-5h)`. If `h!=0`, its first
coefficient is `c_(7h)=0`; if `h=0`, the contribution is `c_0^2=36/49`.
Fejer approximation justifies the calculation for the step function. Thus

```text
integral_0^1 F(x,y)dx=0                 for every y.    (37)
```

This is the load-bearing cyclotomic cancellation.

### 5.2. Exact density quasipolynomial

Let

```text
L_k=14u_kv_k,
k=7q+r,                  0<=r<=6,                       (38)
d_0,...,d_6=(16,0,-16,-18,-6,6,18).                    (39)
```

Then the joint density is exactly

```text
beta_(u_k,v_k)=36/49+d_r/(7L_k).                        (40)
```

To prove `(40)`, let `D_m=(R/Z)\G_m` up to null endpoints and write
`gamma_k=mu(D_(u_k) intersect D_(v_k))`. The exact two-comb overlap formula
from THM-4228 is

```text
L_k gamma_k=W_k,
W_k=w(0)+2sum_(t>=1)w(t),
w(t)=[min(2u_k,u_k+v_k-14t)]_+.                         (41)
```

For `k=7q+r`, the plateau ends at

```text
floor((v_k-u_k)/14)=q,                                  (42)
```

and the positive support ends at

```text
h=6q+b_r,
(b_0,...,b_6)=(0,1,2,2,3,4,5).                         (43)
```

Consequently,

```text
W_k=2u_k+2[q(2u_k)+sum_(t=q+1)^h(u_k+v_k-14t)].         (44)
```

Substitution of `k=7q+r` into the finite arithmetic sum gives, residue by
residue,

```text
7W_k-2u_kv_k=d_r.                                       (45)
```

Finally, inclusion-exclusion and `L_k=14u_kv_k` give

```text
beta_(u_k,v_k)-36/49
 =5/7+W_k/L_k-36/49
 =(7W_k-2u_kv_k)/(7L_k),                                (46)
```

which is `(40)`. Since `|d_r|<=18`, `L_k>=980` for `k>=1`, and

```text
36/49-66/91=6/637>9/3430>=18/(7L_k),                   (47)
```

every member of the family satisfies `beta_(u_k,v_k)>beta_0`.

### 5.3. Exact primitive law and the sharp gate transition

The cyclotomic cancellation admits an exact event proof, not merely the
coarse norm estimate below. Retain `L_k=14u_kv_k`, and put

```text
S_k=L_k beta_(u_k,v_k)=(72u_kv_k+d_r)/7.              (48a)
```

For a tick `x in [0,L_k]`, let `Q_k(x)` be the joint-safe tick length in
`[0,x]` and define the two opposite raw primitives

```text
C_k(x)=L_k Q_k(x)-S_k x=L_k^2 H_(u_k,v_k)(x/L_k),
D_k(x)=S_k x-L_k Q_k(x)=-C_k(x).                       (48b)
```

For `k=7q+r`, define the following seven quadratics:

```text
R_0=25284q^2+ 2997q+   67,
R_1=25284q^2+10011q+  990,
R_2=25284q^2+17368q+ 2924,
R_3=25284q^2+24823q+ 6045,
R_4=25284q^2+32033q+10136,
R_5=25284q^2+39243q+15231,
R_6=25284q^2+46453q+21330.                             (48c)
```

Then, for every `k>=1`, the exact primitive law is

```text
max C_k=-min C_k=A_k=2u_k R_r(q),
omega_(u_k,v_k)=R_r(q)/(49u_k v_k^2).                  (48d)
```

The minimizing tick for `C_k` is

```text
x_- =u_k(98q-13) mod L_k,             r=0,1,
x_- =u_k(98q+85),                     r=2,...,6,       (48e)
x_+ =L_k-x_-.
```

Here `x_+` is the maximizing tick. The modular convention includes `k=1`,
where `x_-=889` and `x_+=91` on the grid `L_1=980`.

#### 5.3.1. Seven finite macrostates

The identity `5v_k-7u_k=1` is the entire reason the exact event order is
finite. On the first half-circle the macro index ranges through

```text
0<=h<=q+floor((5q+r)/2).                               (48e')
```

Index a safe tooth of the `v_k`-comb by `j=7h+s`, `0<=s<=6`. In common
ticks it is

```text
I_(h,s)=[u_k(98h+14s+1),u_k(98h+14s+13)].             (48f)
```

On the first half-circle, the only `u_k`-danger tooth that can meet this
interval has index `i=5h+a_s`, where

```text
(a_0,...,a_6)=(0,1,2,2,3,4,5).                        (48g)
```

Indeed, `u_k/v_k=5/7-1/(7v_k)`. The center displacement, `u_k`-danger
center minus `I_(h,s)` center, is

```text
s=0: 14h-35k-14,       s=1: 14h-7k,
s=2: 14h+21k+14,       s=3: 14h-49k-14,
s=4: 14h-21k,          s=5: 14h+7k+14,
s=6: 14h+35k+28.                                      (48h)
```

The safe and danger half-widths are `6u_k` and `v_k`, so their difference
and sum are `23k+9` and `37k+15`. Comparing `(48h)` with those two numbers
gives the complete overlap word

```text
                 s=0  1  2  3  4  5  6
h<=q-1            L   C  C  0  C  C  R
h=q, r<=2         L   C  R  0  C  C  0
h=q, r>=3         L   C  C  0  C  C  0
h>=q+1            L   C  R  0  C  C  0,              (48i)
```

where `L,C,R,0` mean partial-left, contained, partial-right, and disjoint.
The nonzero overlap lengths are

```text
O_0=14h+2k+1,       O_C=14k+6,
O_(2,R)=16k+1-14h,  O_(6,R)=2k-13-14h.                (48j)
```

Thus a seven-tooth block contains `360k+156` joint-safe ticks before the
boundary and `360k+148` after it. This is a finite event-order derivation of
the purported quasipolynomial; no interpolation or eventuality assertion is
being used.

Below, `V_s` and `R_s` denote respectively the left and right endpoint of
`I_(h,s)`, while `U_s^+` and `U_s^-` denote the intervening `u_k`-comb safe
start and safe end in this block. This notation names endpoint events, not
new runners.

#### 5.3.2. The global positive deficit

The function `D_k` rises on unsafe cells and falls on safe cells. Its local
maxima are therefore safe starts. Let `X_h` be its value at the `s=6` safe
start. Directly subtracting the ten other possible starts in the pre-block
and using

```text
b_k=S_k/L_k=beta_(u_k,v_k)>=2511/3430                 (48k)
```

gives the following positive lower bounds after division by `L_k`:

| candidate | pre-block lower bound for `(X_h-D)/L_k` |
|---|---:|
| `U^+_0` | `(13738k+11167)/3430` |
| `V_1` | `(599k+318)/49` |
| `U^+_1` | `3(1432k+1031)/686` |
| `V_2` | `6(285k+163)/245` |
| `U^+_2` | `(29222k+19763)/3430` |
| `V_3` | `(425k+366)/245` |
| `V_4` | `2(1285k+612)/245` |
| `U^+_4` | `9(296k+1149)/3430` |
| `V_5` | `(1285k+612)/245` |
| `U^+_5` | `(10406k+14639)/3430` |

In the post-block the corresponding list is

```text
(3021k-5983)/3430,       (354k+73)/49,
(10763k-1685)/3430,      (485k-247)/245,
(425k+366)/245,          2(1285k+612)/245,
9(296k+1149)/3430,       (1285k+612)/245,
(10406k+14639)/3430.                                  (48l)
```

It is positive for `k>=2`; `k=1` is the displayed direct control. The common
`r>=3` hybrid starts have one of the pre-block types already listed. Thus the
only blockwise contender is `X_h`. From `(48j)` and `(48a)`,

```text
X_(h+1)-X_h=L_k(2+d_r/v_k)>0             before,
X_(h+1)-X_h=L_k(-4+d_r/v_k)<0            after.        (48m)
```

At the one hybrid boundary, `X_q-X_(q-1)` is respectively

```text
-98(21q-1)(35q+2),       -98(5q+1)(49q+10),
 14(35q+12)(49q+1),       28(35q+17)(49q+15),
196(7q+4)(35q+22),        28(35q+27)(49q+41),
 28(35q+32)(49q+54)                                  (48n)
```

for `r=0,...,6`. The next difference is negative in every residue; its seven
factorizations are

```text
-56(35q+2)(49q-1),       -392(5q+1)(49q+10),
-392(7q+3)(35q+12),       -28(35q+17)(98q+57),
 -28(35q+22)(98q+65),     -28(35q+27)(98q+73),
 -28(35q+32)(98q+81).                                (48o)
```

Consequently the global safe-start maximum is at `h=q-1` for `r=0,1` and
at `h=q` for `r=2,...,6`. Substitution gives `D_k(x_-)=A_k` with `(48c)`.

#### 5.3.3. The opposite extremum

It remains essential not to confuse the largest positive deficit with the
oscillation. Let `E_h` be `D_k` at the `s=0` safe end and put
`E_0=D_k(13u_k)`. In a pre-block, subtracting `E_h` from the ten other safe
ends gives the positive lower bounds

```text
9(952k+197)/3430,       (1285k+612)/245,
(826k-2525)/3430,       2(1285k+612)/245,
(425k+366)/245,         3(9128k+2299)/3430,
6(285k+163)/245,        23(854k+113)/3430,
(599k+318)/49,          (11900k-1699)/3430,           (48p)
```

for `k>=7`; `k<=6` is an exact direct check. In a post-block the nine exact
gaps from `E_h`, in endpoint order, are

```text
-14(1-b_k)h+(-16+26b_k)k-9+13b_k,
-2((23-35b_k)k+9-14b_k),
-14(1-b_k)h+(-90+124b_k)k-41+55b_k,
-14h+(-150+210b_k)k-65+84b_k,
(-28+14b_k)h+(-152+222b_k)k-74+97b_k,
-14h+(-196+280b_k)k-83+112b_k,
(-28+14b_k)h+(-226+320b_k)k-106+139b_k,
-14h+(-242+350b_k)k-101+140b_k,
-14h+(-302+420b_k)k-125+168b_k.                      (48q)
```

The pre-block baselines themselves satisfy

```text
E_(h+1)-E_h=L_k(2+d_r/v_k)>0,                          (48q')
```

so every pre-block baseline is at least `E_0`; the pre-block table includes
the common candidates of the `r>=3` hybrid. For the post-block comparison,
write `h=q+n`, `n>=1`. Then

the block baseline is

```text
(E_h-E_0)/L_k=2q+10n+(q+n)d_r/v_k,             r<=2,
(E_h-E_0)/L_k=2q+10n+5-2r+(q+n)d_r/v_k,        r>=3. (48r)
```

After adding `(48r)` to `(48q)`, three rows increase in `n` and six decrease.
It is therefore enough to test `n=1` and
`n=floor((5q+r)/2)`. Clearing positive denominators leaves coefficientwise
positive polynomials, except for the following seven elementary survivors:

```text
392q^2-81q-3,
512785q^3+41587q^2+524q-4,
344715q^3+7973q^2-904q-4,
8q-1,                 392q^2+143q-1,
72030q^3+111916q^2+43088q-297,
72030q^3+124656q^2+44094q-8847.                       (48s)
```

Each is positive at `q=1` and has positive forward difference thereafter.
The boundary `r<=2` uses the same calculation at `n=0`. For `r>=3`, the
only extra hybrid end has positive gap from `E_0` equal to

```text
28(7q+3)(35q+17)(119q+66),
28(35q+22)(833q^2+1056q+334),
28(35q+27)(833q^2+1293q+502),
28(35q+32)(833q^2+1530q+702)                          (48t)
```

in residues `3,4,5,6`. Hence `E_0` is the first-half minimum. Explicitly,

```text
E_0=u_k(13d_r-7028k^2-5686k-1146)/7<0.               (48u)
```

Finally, `A_k+E_0` factors, for `r=0,...,6`, as

```text
28q(35q+2)(49q+11),          196q(5q+1)(49q+10),
14(35q+12)(98q^2+67q+1),      28(q+1)(35q+17)(49q+15),
196(q+1)(7q+4)(35q+22),       28(q+1)(35q+27)(49q+41),
28(q+1)(35q+32)(49q+54).                              (48v)
```

Thus `-E_0<=A_k`, with equality only for `k=1`. Reflection of both combs
gives `D_k(L_k-x)=-D_k(x)`, so the second half supplies `-A_k` and no larger
absolute extremum. This completes the proof of `(48d)--(48e)`.

The exact formula is strictly decreasing with `k`. Cross-multiplying the
successive residue expressions in `(48d)` gives positive denominators and,
up to positive common factors, the following coefficient-positive
numerators (the last row compares residue six at `q` with residue zero at
`q+1`):

```text
r=0: 2(3277365q^3+1000825q^2+78768q+1454),
r=1:   5966485q^3+5339726q^2+1527989q+138652,
r=2: 284122335q^4+473196283q^3+290531045q^2+77845533q+7667748,
r=3: 304710910q^4+682023258q^3+568329489q^2+208868341q+28551678,
r=4: 304710910q^4+857824478q^3+902405315q^2+420381265q+73168566,
r=5: 304710910q^4+1033625698q^3+1311680461q^2+738025641q+155354760,
r=6: 312946340q^4+1239225729q^3+1836580907q^2+1207368844q+297073440.
                                                                    (48w)
```

At the unique crossing,

```text
omega_747=4575651/79563540224,
omega_747-T=12591073311/326326757250667264>0,

omega_748=144518186/2516324606159,
T-omega_748=336393488/8724097409553253>0.              (48x)
```

Therefore

```text
omega_(u_k,v_k)<=T                  iff k>=748.         (48y)
```

Together with `(47)`, Section 2 proves `(5)` for every `k>=748` and every
common dilation `g>=1`; THM-4150 gives the odd-tail consequence `(6)`. The
failure at `k=747` is failure of this sufficient transfer gate, not an
unsafe-set witness.

### 5.4. Coarse analytic mechanism

By `(37)`, define the one-periodic fast primitive

```text
P(x,y)=integral_0^x F(s,y)ds.                           (48)
```

For fixed `y`, the uncentered indicator in `(35)` has mean `36/49`.
THM-4228's sharp primitive estimate therefore gives

```text
osc_x P(.,y)<=(36/49)(13/49)=468/2401.                 (49)
```

Moreover, translating the first safe comb by `2 delta` changes it on measure
at most `4|delta|`, and translating the second by `3 delta` changes it on
measure at most `6|delta|`. Hence

```text
|P(x,y)-P(x,y')|<=10|y-y'|.                             (50)
```

The function `P` is Lipschitz and piecewise affine; hence its absolutely
continuous restriction to the line `(kt,t)` obeys the chain rule away from
finitely many wall crossings, and integration yields

```text
integral_0^t F(ks,s)ds
 =P(kt,t)/k-(1/k)integral_0^t partial_y P(ks,s)ds.      (51)
```

Because `P(0,y)=0`, `(49)` implies `|P(x,y)|<=468/2401`. The first term on
the right of `(51)` therefore has oscillation at most `936/(2401k)`, while
the primitive of `partial_y P` has oscillation at most `10`. Centering from
`36/49` to the exact density `(40)` adds at most
`|d_r|/(7L_k)`. Thus

```text
omega_(u_k,v_k)
 <=24946/(2401k)+18/(7L_k)
 <=24946/(2401k)+9/(1715k^2).                           (52)
```

The right side decreases with `k`. At `k=180785`, it is

```text
22549313113/392362010781125,                            (53)
```

and the exact threshold margin is

```text
T-22549313113/392362010781125
 =28971847951/229893926442909103375>0.                  (54)
```

Equations `(47)`, `(52)`, and `(54)` give an independent, deliberately coarse
analytic proof of the gate for every `k>=180785`. The finite macrostate
argument in Section 5.3 is what sharpens this to the exact transition `748`.

### 5.5. General unimodular resonance-zero template

The cancellation above is not isolated. Let `a,c` be positive integers and
`b,d` be integers such that

```text
ad-bc=+-1,                  and                  7|a or 7|c. (55)
```

For every integer `k>=1` for which

```text
u=ak+b,                     v=ck+d                         (56)
```

are distinct and positive, they are coprime, because every common divisor
divides `cu-av=cb-ad=+-1`. Define

```text
F_(a,b,c,d)(x,y)=j(ax+by)j(cx+dy)-36/49.                (57)
```

The `x`-resonances are `(ch,-ah)`. Condition `(55)` makes one of their two
nonzero Fourier coefficients a multiple-of-seven coefficient of `j`, hence

```text
integral_0^1 F_(a,b,c,d)(x,y)dx=0.                       (58)
```

The same primitive proof as `(48)--(52)` has fast oscillation at most
`468/2401` and slow Lipschitz constant at most `2(|b|+|d|)`. On the other
hand, Fourier orthogonality in the two danger combs forces the modes
`(vh,-uh)`; product-to-sum and
`sum_(h>=1)cos(2pi hz)/h^2=pi^2 B_2({z})` give

```text
beta_(u,v)-36/49
 =[B_2({(v-u)/14})-B_2({(v+u)/14})]/(uv),               (59)
```

where braces denote fractional part and `B_2(z)=z^2-z+1/6` on `[0,1)`.
Since
`-1/12<=B_2<=1/6`, equations `(58)--(59)` prove the uniform compiler bounds

```text
|beta_(u,v)-36/49|<=1/(4uv),
omega_(u,v)
 <=[936/2401+2(|b|+|d|)]/k+1/(4uv).                    (60)
```

Whenever the lower density and upper oscillation in `(60)` pass `(14)`, the
THM-4228 deck supplies every common dilation. The family `(4)` is the
instance `(a,b,c,d)=(5,2,7,3)`; its exact density `(40)` and primitive law
`(48d)` yield the sharp displayed cutoff.

## 6. Exact audit and reproduction

The fixed-pair computation has the explicit finite universe `(17)--(21)`:
all joint endpoints, all resulting cells, and every primitive prefix. Its
inherited filters are exactly the two safe-comb predicates; there is no
search truncation or floating comparison. The positive control `(1,2)`, the
sharp-density hostile `(1,13)`, and the failed-gate Fibonacci predecessor in
`(25)` test the density, primitive normalization, and inequality direction.

For the fixed pair, the primary implementation uses the common integer grid
and midpoint cell classification. The independent implementation builds
reduced-rational safe teeth, intersects intervals directly, and integrates
its own rational primitive. For the canceled family, a new common-grid path
checks the formula at every `k=1,...,1000` and at `5000,10000,50000`; its
independent path intersects the two ordered safe-tooth families without
constructing or classifying the endpoint-cell union, checking every
`k=1,...,850` and the same far controls. These finite sweeps audit, but do
not replace, the infinite macrostate proof in Section 5.3. Reproduce from the
repository root with

```bash
clang++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.cpp \
  -o /tmp/thm4233-primary
/tmp/thm4233-primary

clang++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_independent_audit_thm4233.cpp \
  -o /tmp/thm4233-independent
/tmp/thm4233-independent

python3 \
  04-computation/lrc14_resonance_zero_exact_tail_thm4233.py

python3 \
  04-computation/lrc14_resonance_zero_exact_tail_thm4233_independent_audit.py
```

All four output streams must byte-match their frozen files.
Repeating both compilations with `-O0` gives the same two output streams; the
independent executable also passes combined AddressSanitizer/UBSan replay.

## 7. Connection contract, scope, and open boundary

```text
source:       one primitive pair observable A_(u,v) and the THM-4228 deck
target:       a lawful fixed-pair or canceled-family depth-eight repair
map:          A_(u,v) -> (beta,centered primitive H,omega), then pull back
              by the common dilation and intersect each base component
preserved:    literal pair ratio, exact density, endpoint phase, base mass,
              component count, common dilation, and target threshold
destroyed:    endpoint addresses after reduction to (beta,omega)
sidecar:      exact primitive extremizers or an analytic fast primitive
decisive test: beta>=66/91 and omega/g<=1650/28710227.   (61)
```

The proof and failure boundaries are:

1. The exact computation certifies `(3713,5149)` and its dilations; it makes
   no claim that this is the first, smallest, or unique coprime pair passing
   the gate.
2. Equation `(27)` is an obstruction to one proof mechanism, not a claim
   that adjacent pairs are unsafe. Common dilation still reduces their
   discrepancy, and THM-4231 handles their sufficiently large coordinates.
3. The canceled family has the exact period-seven primitive law `(48d)` and
   passes this sufficient gate exactly for `k>=748`. For `k<=747`, gate
   failure is not danger; another repair or direct safe-set argument may work.
4. The deck `E_3467` is sufficient, not necessary. Missing `(14)` says only
   that this transfer does not fire.
5. Arbitrary finite pair entry, replacement of the fixed pool, and full
   LRC(14) remain open.

**QED.**
