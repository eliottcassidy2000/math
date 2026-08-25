---
id: THM-4092
title: "Parity-weighted antipodal k-comb density compiler"
status: >
  PROVED + PROVED RELATIVE TO THM-2061/2066/2072 + VERIFIED-EXACT +
  INDEPENDENTLY VERIFIED-EXACT. An even antipodal danger comb has density
  weight one and an odd comb weight two. A safe interval J of length L
  absorbs any finite outlier bank V of total weight W<7 whenever
  sum_(v in V)(7-omega(v))/v < 7L(7-W), and compiles an adaptive even clock.
  Exact AP5--AP8 intervals yield new eleven-core families with three through
  six outliers; in particular AP8 plus any three distinct speeds at least 281
  closes every dyadic two-odd-tail seam through all common dilations. This is
  not LRC(14), an arbitrary-core supplier, or a converse at weighted arity 7.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-883-fragmentation-lemma-726-rigorized
  - THM-1155-threespeed-exhaustive-and-ceiling
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
  - THM-4087-lrc-literal-open-two-comb-compiler
related:
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-4081-lrc-antipodal-height12-obstruction-and-six-speed-floor
  - MISTAKE-126
  - MISTAKE-274
script: 04-computation/lrc_parity_weighted_antipodal_kcomb_thm4092.py
output: 05-knowledge/results/lrc_parity_weighted_antipodal_kcomb_thm4092.out
independent_audit_script: 04-computation/lrc_parity_weighted_antipodal_kcomb_thm4092_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_parity_weighted_antipodal_kcomb_thm4092_independent_audit.out
script_sha256: cc9f1e1d210cb9cb0f8da356d1f391cbf5d196548d7c8dfa9f341f61913757cf
output_sha256: fd01f330831e7df484583323bc1e44bf7ed5cfd1059c59726ca801a03d1f68e5
independent_audit_script_sha256: 8aab44a830b07b614f8460be07df70b2659a3b6874a75fa920b3e12d666dc7d7
independent_audit_output_sha256: 2880e844bf7e19ce785f83511f457bb40cc6f05ae250f98981e29c705eed1dd8
hash_basis: raw LF bytes
---

# THM-4092 -- the parity-weighted antipodal comb compiler

Put `delta=1/14`. For a finite set `D` of positive integer speeds, retain the
closed two-phase safe set

```text
G_D^+-={theta in R/Z:
          ||d theta||>=delta and ||d(theta+1/2)||>=delta
          for every d in D}.                                      (1)
```

THM-4087 treated two new speeds by bounding connected components of their
literal danger union. The present theorem uses the least-used sidecar from
THM-1155 instead: interval fragmentation. The necessary new coordinate is
parity. It permits more outliers, but its conclusion is a density certificate,
not a component classification.

## 1. Parity is the exact antipodal density weight

For a speed `v`, let

```text
U_v={theta:min(||v theta||,||v(theta+1/2)||)<delta}.          (2)
```

Define

```text
omega(v)=1  if v is even,
omega(v)=2  if v is odd.                                    (3)
```

The literal open teeth of `(2)` have centres and half-width

```text
centre k/(omega(v)v),       half-width 1/(14v).              (4)
```

Indeed, when `v` is even the two phases in `(2)` agree modulo one. When `v`
is odd they differ by a half-turn, so their two root sets interlace to the
`2v`-grid. Thus the period, tooth length, and duty cycle are

```text
p_v=1/(omega(v)v),
a_v=1/(7v),
rho_v=a_v/p_v=omega(v)/7.                                   (5)
```

This is why an odd outlier costs exactly twice an even outlier. A naked order
on the speeds retains none of `(3)--(5)`.

## 2. Sharp one-period interval discrepancy

Let `I` be any lifted real interval of length `L`. Then

```text
boxed:
mu(I intersect U_v)
 <= omega(v)L/7 + (7-omega(v))/(49v).                        (6)
```

To prove `(6)`, write `L=mp_v+r`, with `0<=r<p_v`. Every interval of length
`p_v` contains exactly `rho_v p_v` danger mass. An interval of residual
length `r` meets at most

```text
min(r,rho_v p_v)                                             (7)
```

danger mass: if it meets two adjacent teeth it must also span the intervening
gap, and its total tooth portion is still at most `(7)`. The maximum of

```text
min(r,rho_v p_v)-rho_v r
```

over `0<=r<=p_v` is `rho_v(1-rho_v)p_v`. Substitution from `(5)` gives

```text
rho_v(1-rho_v)p_v=(7-omega(v))/(49v),                       (8)
```

which proves `(6)`. The constant is sharp for every `v`: take `I` to be one
whole tooth, of length `rho_v p_v`, centred at a root. Both exact referees
check this equality separately for every tested speed.

The coarser THM-1155 error `1/(7v)` is valid but misses the factors `6/7` and
`5/7` in `(8)`. Those factors materially improve the all-height thresholds.

## 3. The weighted absorption criterion

Let `J subset G_D^+-` be a closed circular interval of length `L>0`, and let
`V` be a finite set of distinct positive speeds outside `D`. Put

```text
W=sum_(v in V) omega(v).                                     (9)
```

Assume

```text
boxed:
W<7
and
sum_(v in V) (7-omega(v))/v < 7L(7-W).                     (10)
```

Then

```text
boxed: J minus union_(v in V) U_v is nonempty.              (11)
```

If the union covered `J`, subadditivity and `(6)` would give

```text
L <= (W/7)L + (1/49)sum_(v in V)(7-omega(v))/v,
```

or equivalently the reverse weak inequality to `(10)`. This contradiction
proves `(11)`. Every point in `(11)` lies in `G_(D union V)^+-`.

The exact reciprocal condition `(10)` is the primary theorem. It permits a
small outlier to be compensated by larger ones. A convenient uniform
corollary is the following. Let `k=|V|`, let `o` be the number of odd speeds,
so `W=k+o`, and suppose every `v in V` is at least `B`. Since

```text
sum_(v in V)(7-omega(v))/v <= (7k-W)/B,                    (12)
```

it suffices that

```text
boxed:
k+o<7
and
B > (7k-k-o)/(7L(7-k-o)).                                  (13)
```

The strict inequality is part of the certificate. No conclusion is asserted
at equality.

## 4. A survivor compiles an adaptive even owner clock

Suppose the endpoints of a real lift of `J` have even rational
presentations with denominators at most `Q`. The closed set in `(11)` is a
finite union of intervals and points. Choose an endpoint of one of its
components. It is either an endpoint of `J` or a tooth endpoint

```text
(14j +/- 1)/(14v)       for even v,
( 7j +/- 1)/(14v)       for odd v.                          (14)
```

Therefore the survivor has an even presentation

```text
H(V)=max({0} union V),
theta=r/N,        N<=max(Q,14H(V)).                          (15)
```

At both labels `r` and `r+N/2`, every speed in `D union V` is weak-safe. For
every odd tail class `z`, the two clock distances are complementary:

```text
|z(r+N/2)|_N=N/2-|zr|_N.                                   (16)
```

They cannot both be strictly below `N/7`. THM-2066 eligibility requires
exactly those two strict inequalities at every safe label, hence

```text
boxed: E_N(D union V)=R_N(D union V)=empty.                 (17)
```

By THM-2061/2066/2072, when `D union V` has eleven speeds, `(17)` closes every
dyadic seam with two distinct positive odd tails.

## 5. Four exact AP-core intervals

For odd `d`, the two-phase condition is exactly

```text
1/14 <= ||d theta|| <= 3/7;                                 (18)
```

for even `d` it is just `||d theta||>=1/14`. The following table gives four
closed intervals contained in the corresponding AP safe sets.

| core | interval `J_m` | length `L_m` | endpoint clock bound `Q_m` |
|---|---|---:|---:|
| `{1,...,8}` | `[11/49,13/56]` | `3/392` | `98` |
| `{1,...,7}` | `[4/35,13/98]` | `9/490` | `98` |
| `{1,...,6}` | `[4/35,1/7]` | `1/35` | `70` |
| `{1,...,5}` | `[4/35,1/7]` | `1/35` | `70` |

For completeness, the exact norm ranges on these intervals are:

```text
AP5/AP6:
d=1 [4/35,1/7]   d=2 [8/35,2/7]   d=3 [12/35,3/7]
d=4 [3/7,1/2]    d=5 [2/7,3/7]    d=6 [1/7,11/35]

AP7:
d=1 [4/35,13/98] d=2 [8/35,13/49] d=3 [12/35,39/98]
d=4 [16/35,1/2]  d=5 [33/98,3/7]  d=6 [10/49,11/35]
d=7 [1/14,1/5]

AP8:
d=1 [11/49,13/56] d=2 [22/49,13/28] d=3 [17/56,16/49]
d=4 [1/14,5/49]   d=5 [6/49,9/56]  d=6 [17/49,11/28]
d=7 [3/8,3/7]     d=8 [1/7,10/49].                          (19)
```

Each range follows by checking the interval endpoints together with the
integer and half-integer breakpoints of the tent function `||d theta||`.
Equations `(18)--(19)` prove all four inclusions without a sampling argument.

## 6. Explicit eleven-core families

Let `m in {5,6,7,8}`, put `k=11-m`, and adjoin distinct speeds

```text
v_1<...<v_k,       v_1>=B,                                  (20)
```

to `{1,...,m}`. The smallest integer `B` supplied by `(13)` is:

| core | outliers `k` | number odd `o` | sufficient `B` |
|---|---:|---:|---:|
| AP8 | 3 | 0 | 85 |
| AP8 | 3 | 1 | 106 |
| AP8 | 3 | 2 | 150 |
| AP8 | 3 | 3 | **281** |
| AP7 | 4 | 0 | 63 |
| AP7 | 4 | 1 | 90 |
| AP7 | 4 | 2 | 172 |
| AP6 | 5 | 0 | 76 |
| AP6 | 5 | 1 | 146 |
| AP5 | 6 | 0 | 181 |

The missing parity rows have `k+o>=7` and are not certified by this method.
For every displayed row, there is an even clock

```text
N<=14v_k,       E_N=R_N=empty.                              (21)
```

The all-parity headline is particularly simple:

```text
K={1,...,8,B,C,D},       281<=B<C<D                         (22)
```

always has a two-phase safe point and a clock `N<=14D`. Consequently, for
every `q>=1` and every two distinct positive odd integers `x,y`,

```text
2qK union {x,y}                                             (23)
```

is `1/14`-lonely.

To justify the common dilation in `(23)`, if `theta` is safe for the core,
use `theta/q`. At phase zero every product is unchanged. At the half phase,
an odd `q` preserves the old half phase and an even `q` collapses it to the
old zero phase. Both were safe. The clock becomes `qN<=14qv_k`. This is
forward certificate transport; it is not an even-dilation obstruction
equivalence.

## 7. Failure boundaries and hostile controls

The weighted ceiling in `(10)` is exact for the first-order method. At
`W=7`, the positive coefficient of `L` vanishes. The AP7 interval is in fact
literally covered by the four outliers

```text
{8,9,11,13},        weights 1+2+2+2=7.                      (24)
```

This does not prove that all weight-seven banks cover; it proves that the
ceiling is not a cosmetic algebraic artifact.

Each displayed parity row also has a low-height hostile which literally
covers the chosen core interval:

```text
AP8: o=0 (10,22,26), o=1 (9,10,26),
     o=2 (9,10,11), o=3 (9,11,13);
AP7: o=0 (8,10,12,26), o=1 (8,9,10,12),
     o=2 (8,9,10,11);
AP6: o=0 (8,10,14,22,26), o=1 (7,8,10,12,26);
AP5: o=0 (6,8,10,14,22,26).                                (25)
```

These are hostiles to this selected-interval certificate only. They are not
claimed to be antipodal obstructions, LRC counterexamples, or sharp global
height thresholds.

The theorem also does not supersede THM-4087: its connected-component proof
gives the much better all-parity AP9 two-outlier threshold `70`. Density is
the correct sidecar only once component chains make the two-comb argument
unavailable. THM-4081's zero-width hostile still shows why some positive core
interval is necessary.

## 8. Verification and connection ledger

Primary reproduction:

```bash
python3 -B 04-computation/lrc_parity_weighted_antipodal_kcomb_thm4092.py
python3 -B -O 04-computation/lrc_parity_weighted_antipodal_kcomb_thm4092.py
```

Independent reproduction:

```bash
python3 -B 04-computation/lrc_parity_weighted_antipodal_kcomb_thm4092_independent_audit.py
python3 -B -O 04-computation/lrc_parity_weighted_antipodal_kcomb_thm4092_independent_audit.py
```

The primary referee checks `34,000` exact hostile interval windows, `80`
sharp equalities, all ten threshold rows, `6,370` odd owner classes, dilation
controls, and all eleven literal covers `(24)--(25)`. The independent referee
uses original-theta cell classification rather than interval merging; it
checks `15,300` windows, `60` sharp equalities, `774` threshold-bank rows,
`7,518` owner classes, and the same hostile covers. Normal and optimized
outputs are byte-identical; neither script uses `assert` or floating-point
arithmetic.

The typed connection is:

```text
source:       one metric two-phase safe interval J
operation:    delete finitely many literal antipodal danger combs
preserved:    weak two-phase safety
lost data:    which tooth/component supplied the survivor
sidecar:      parity weight omega(v) plus reciprocal boundary tax
recovery:     choose a surviving arrangement endpoint
target:       an adaptive even owner clock with E_N=R_N=empty
cheap test:   exact wall sweep plus complementary odd-class replay.        (26)
```

Thus the advance is an infinite structured supplier for the dyadic seam, not
a proof of arbitrary-core supply, physical entry, or LRC(14).
