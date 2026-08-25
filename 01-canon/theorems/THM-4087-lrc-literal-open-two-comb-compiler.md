---
id: THM-4087
title: "LRC literal-open two-comb compiler"
status: >
  PROVED + PROVED RELATIVE TO THM-2061/2066/2072 + VERIFIED-EXACT +
  INDEPENDENTLY VERIFIED-EXACT. A closed antipodal-safe core interval of
  length L absorbs every two distinct new speeds B<C once B>=2/(7L).
  If its endpoints have bounded even rational presentations, a surviving
  arrangement endpoint gives an adaptive even owner clock
  N<=max(Q,14C) with E_N=R_N=empty. In particular every eleven-core
  {1,...,9,B,C}, 70<=B<C, closes every dyadic two-odd-tail seam, at a
  clock N<=14C; the certificate transports under every positive common
  dilation. This is an infinite two-outlier seam family, not LRC(14).
source: codex-frontier-synthesis-creative-20260825f / two-outlier clocks anchor
audit: >
  PASS. The primary Fraction-exact antipodal-comb engine checks 3,160
  literal-open component pairs, 499,500 all-parity proof branches, 6,105
  AP9 two-outlier rows, 1,756,510 odd owner classes, dilation controls,
  the H8+7m zero-width hostile, and the 48/123 local cover. The independent
  original-theta engine constructs both literal half-shift walls, checks
  1,225 wall/cell component pairs, 3,240 family rows, 2,004,002 odd owner
  classes, all twelve H8 endpoint walls, MISTAKE-274's open handoff, and
  the exact 48/123 hostile component. Normal and optimized outputs
  byte-match; both scripts have zero assert nodes and zero floating literals.
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
related:
  - THM-1094-exact-two-comb-component-theorem
  - THM-2440-sharp-two-comb-centred-window-radius
  - THM-4075-lrc14-divisor-complete-dyadic-owner-word-closure-through-30
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-4081-lrc-antipodal-height12-obstruction-and-six-speed-floor
  - MISTAKE-126
  - MISTAKE-273
  - MISTAKE-274
  - MISTAKE-382
  - MISTAKE-464
script: 04-computation/lrc_literal_open_two_comb_compiler_thm4087.py
output: 05-knowledge/results/lrc_literal_open_two_comb_compiler_thm4087.out
independent_audit_script: 04-computation/lrc_literal_open_two_comb_compiler_thm4087_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_literal_open_two_comb_compiler_thm4087_independent_audit.out
script_sha256: 1766d3dcc517be0128d65ac67a1a0099fe536b0fa61095412a23875250621374
output_sha256: ed506bcec5685b32ebfeb24ccb9f1bdd7615b57238c86d4c75f710384ad0bd92
independent_audit_script_sha256: 3fe60d9581fec2bbac7815bda45226cba858d238cc8c7eb0caebf16ccedaba4a
independent_audit_output_sha256: 49598d183205f1c1f34dbf298c900eef32ba839f2c35050f571c4bccf4c62036
hash_basis: raw LF bytes
---

# THM-4087 -- a literal-open two-comb compiler

Put `delta=1/14`. For a finite set `D` of positive integer speeds, define its
closed antipodal-safe set

```text
G_D^+-={theta in R/Z:
          ||d theta||>=delta and ||d(theta+1/2)||>=delta
          for every d in D}.                                (1)
```

The theorem uses the metric length of one closed interval in `(1)`. This is
strictly more faithful than retaining only one phase or the measure of the
whole safe set.

## 1. Literal antipodal danger combs

For a positive speed `v`, let

```text
U_v={theta:min(||v theta||,||v(theta+1/2)||)<delta}.          (2)
```

These are the literal **open** failures of weak antipodal safety. Their tooth
description depends on parity:

```text
v even: U_v = union_(k mod v)
                  (k/v-1/(14v), k/v+1/(14v)),

v odd:  U_v = union_(k mod 2v)
                  (k/(2v)-1/(14v), k/(2v)+1/(14v)).          (3)
```

Thus every tooth has length `1/(7v)`. Even speeds have root spacing `1/v`,
while odd speeds have root spacing `1/(2v)`. Endpoints in `(3)` are excluded;
two teeth that merely touch remain different open components.

### Literal two-comb component lemma

Let `B<C` be distinct positive integers. Every connected component of

```text
U_B union U_C                                                (4)
```

has circular length strictly less than

```text
boxed: 2/(7B).                                               (5)
```

To prove this, first compare one `C`-tooth with the gap between consecutive
`B`-teeth. The `C`-tooth has length `1/(7C)<1/(7B)`, whereas the smaller of
the two possible `B`-gaps is

```text
1/(2B)-1/(7B)=5/(14B).                                      (6)
```

Therefore no `C`-tooth meets two `B`-teeth. Since distinct `C`-teeth are
disjoint, a component of `(4)` contains at most one `B`-tooth.

If at most one `C`-tooth meets that `B`-tooth, the component length is at
most

```text
1/(7B)+1/(7C)<2/(7B).                                      (7)
```

Suppose instead that at least two `C`-teeth meet it. Their centres lie in the
`B`-tooth enlarged by one `C` half-width, an interval of length

```text
1/(7B)+1/(7C).                                              (8)
```

If `C` is odd, its centre spacing `1/(2C)` is at most `(8)`, forcing

```text
5B<=2C.
```

All attached `C`-teeth extend the `B`-tooth by at most one full `C`-tooth on
each side. Hence the component length is at most

```text
1/(7B)+2/(7C) <= 9/(35B) < 2/(7B).                         (9)
```

If `C` is even, its spacing `1/C` is at most `(8)`, forcing `6B<=C`; the
same extension bound is at most

```text
1/(7B)+2/(42B)=4/(21B)<2/(7B).                            (10)
```

Components containing no `B`-tooth are still shorter. This proves `(5)`.
The proof concerns literal open components, not closed coverage or
almost-everywhere coverage.

## 2. Safe-interval absorption

Let `J` be a closed circular interval of length `L>0` satisfying

```text
J subset G_D^+-.                                            (11)
```

Let `B<C` be two distinct positive speeds outside `D`. If

```text
boxed: B>=2/(7L),                                           (12)
```

then

```text
there exists theta in J with
theta in G_(D union {B,C})^+-.                              (13)
```

Indeed, if every point of `J` failed one of the two new speeds, then the
connected interval `J` would be contained in one connected component of the
open set `U_B union U_C`. That component would have length at least `L`.
But `(5)` and `(12)` give

```text
component length < 2/(7B) <= L,                            (14)
```

a contradiction. The strict inequality in `(14)` is why equality is allowed
in `(12)`. Any internal open-tooth handoff omitted by both combs is already a
valid weak-safe point, rather than a gap in the proof.

### Strict-point corollary

Suppose `R=max D` and a phase `theta_0` has strict margin `eta>0`:

```text
||d(theta_0+j/2)||>=delta+eta
    for d in D and j=0,1.                                  (15)
```

The one-Lipschitz norm makes the closed interval

```text
[theta_0-eta/R, theta_0+eta/R] subset G_D^+-.               (16)
```

It has length `2eta/R`. Therefore two distinct outliers are absorbed whenever

```text
boxed: min(B,C)>=R/(7eta).                                 (17)
```

This is the two-outlier continuation of THM-4079's one-outlier perturbation:
the required height doubles because the proof must escape a literal union of
two combs. Section 4 below shows why a positive interval, not mere nonempty
weak feasibility, is load-bearing.

## 3. Rational endpoints compile an adaptive owner clock

Assume the interval in `(11)` has a chosen real lift with endpoints

```text
a=r_0/N_0,       a+L=r_1/N_1,                              (18)
```

where `N_0,N_1` are positive even integers at most `Q`. Then the phase in
`(13)` may be chosen at an arrangement endpoint. More precisely, the closed
set

```text
J minus (U_B union U_C)                                    (19)
```

is nonempty. An endpoint of one of its components is either an endpoint in
`(18)` or an endpoint of a tooth in `(3)`. In the latter case it has the form

```text
(14k +/- 1)/(14v)       if v is even,
( 7k +/- 1)/(14v)       if v is odd,                       (20)
```

for `v in {B,C}`. Thus there is an even clock and label

```text
theta=r/N,       N<=max(Q,14C),                            (21)
```

such that both `r` and `r+N/2 mod N` lie in THM-2066's labelled safe packet

```text
A_N(D union {B,C}).                                        (22)
```

For every odd tail class `z mod 2N`, the two normalized clock distances in
`(22)` are complementary:

```text
|z(r+N/2)|_N = N/2-|zr|_N.                                (23)
```

They cannot both be strictly below `N/7`. Eligibility requires exactly those
two strict inequalities at every safe label. Consequently

```text
boxed: E_N(D union {B,C})=R_N(D union {B,C})=empty.         (24)
```

This is an adaptive clock manufactured from the surviving arrangement
endpoint; no nearest-integer approximation or residual positive margin is
needed.

By THM-2061/2066/2072, `(13)` or `(24)` implies that no two distinct positive
odd tails `x,y` make

```text
2(D union {B,C}) union {x,y}                               (25)
```

a strict LRC(14) counterexample, whenever the core in `(25)` has eleven
distinct speeds.

## 4. The exact AP9 two-outlier family

Take

```text
D={1,2,...,9}.                                              (26)
```

The connected component of `G_D^+-` containing `1/12` is exactly

```text
J_9=[4/49,3/35],             |J_9|=1/245.                  (27)
```

For odd speeds, antipodal safety is the band
`1/14<=||d theta||<=3/7`; for even speeds it is
`||d theta||>=1/14`. On `(27)` the exact ranges are

| `d` | `min ||d theta||` | `max ||d theta||` |
|---:|---:|---:|
| 1 | `4/49` | `3/35` |
| 2 | `8/49` | `6/35` |
| 3 | `12/49` | `9/35` |
| 4 | `16/49` | `12/35` |
| 5 | `20/49` | `3/7` |
| 6 | `17/35` | `1/2` |
| 7 | `2/5` | `3/7` |
| 8 | `11/35` | `17/49` |
| 9 | `8/35` | `13/49` |

Every required inequality follows. At the left endpoint speed `7` has upper
band equality, and immediately to the left it fails; at the right endpoint
speed `5` has upper equality, and immediately to the right it fails. This
proves the exact component assertion, including both closed endpoints.

Since

```text
2/(7|J_9|)=70,                                              (28)
```

Sections 2--3 prove the explicit all-height family:

```text
C_(B,C)={1,...,9,B,C},        70<=B<C,

there is theta in G_(C_(B,C))^+-,
and an even N<=14C with E_N=R_N=empty.                      (29)
```

The endpoints of `(27)` have the even presentations

```text
4/49=8/98,             3/35=6/70,                          (30)
```

so `(21)` gives the displayed clock bound. In particular, for every two
distinct positive odd integers `x,y`, the thirteen-speed set

```text
2{1,...,9,B,C} union {x,y}                                 (31)
```

is `1/14`-lonely. This replaces THM-4079's one large outlier over an AP10
body by two completely arbitrary-parity outliers over an AP9 body.

## 5. Dilation correctness

The certificate transports under **every** positive common dilation, without
a primitivity assumption. If `q>=1` and `theta` proves `(29)`, put

```text
theta_q=theta/q.                                           (32)
```

For every old speed `d`, the phase at label zero is unchanged:

```text
qd theta_q=d theta.                                        (33)
```

At the half label,

```text
qd(theta_q+1/2)=d theta+qd/2.                              (34)
```

If `q` is odd, `(34)` agrees modulo one with `d(theta+1/2)`; if `q` is even,
it agrees with `d theta`. Both were weak-safe. Hence `theta_q` is an
antipodal-safe phase for `q C_(B,C)`.

If `theta=r/N` is the clock phase in `(29)`, then

```text
theta_q=r/(qN),          qN<=14qC,                         (35)
```

and the same complementary-label proof gives `E_(qN)=R_(qN)=empty` for the
dilated core. Odd dilation preserves the entire parity anatomy; even dilation
is only this forward safety transport, not an obstruction equivalence. This
matches THM-4081's sharp odd/even boundary.

## 6. Failure boundaries and correction lineage

### Positive width is necessary

THM-4081's endpoint hostile

```text
H_8={1,3,5,8,11,13,23,36}                                 (36)
```

has transformed safe set exactly

```text
{1/7,2/7,3/7,4/7,5/7,6/7}.                               (37)
```

For every positive integer `m`, adjoining the arbitrarily large speed `7m`
kills all six points. If `7m` is odd, its transformed multiplier is `7m`; if
it is even, the multiplier is `7m/2`. Either is a multiple of seven, so its
norm at every point in `(37)` is zero. Therefore

```text
H_8 union {7m} is an antipodal obstruction for every m>=1. (38)
```

This is an obstruction to the antipodal certificate, not an LRC
counterexample. It proves that nonempty weak feasibility, even with
arbitrarily large outlier height, cannot replace the interval premise.

### A low two-outlier pair can cover the AP9 component

The threshold `70` is sufficient, not claimed least. The exact low hostile

```text
(B,C)=(48,123)                                              (39)
```

does show that the interval-width proof has a genuine lower boundary. Around
`theta_0=1/12`, speed `48` supplies the central tooth of radius `1/672`.
Speed `123` supplies symmetric teeth centred at offsets `+/-1/(4*123)` and
of radius `1/(14*123)`. Their inner endpoints overlap the central tooth, and
the resulting literal-open component is

```text
(139/1722,74/861)
 =(1/12-3/1148,1/12+3/1148).                              (40)
```

It strictly contains all of `J_9`. The full enlarged core still has safe
components elsewhere, so `(39)` is a counterexample only to blindly deleting
the height hypothesis from this particular interval bridge.

### Open seams are not filled

MISTAKE-274 corrected a closed/almost-everywhere two-comb radius that had been
assigned to a literal open component. The minimal seam is `{1,13}` at
`theta=1/14`: both standard danger inequalities are equalities, so neither
open comb contains the point. The proof of `(5)` uses strict overlap to join
teeth and treats equality as a safe separator. The independent audit replays
that seam and the genuinely overlapping `{1,14}` control. It also tests every
wall, so MISTAKE-382/464's endpoint-only failures cannot be hidden by
midpoint sampling. Formulae `(32)--(35)` explicitly pay MISTAKE-126's
dilation requirement.

## 7. Audits, connection ledger, and scope

The primary companion constructs `(3)` directly as open rational teeth. It
does not merge touching endpoints. It checks all `3,160` pairs `1<=B<C<=80`,
obtaining worst bound ratio `103/104` at `(77,78)`, and verifies the all-scale
proof inequalities over `499,500` ordered pairs through `1000`. It then tests
all `6,105` rows `70<=B<C<=180`, all `134,310` core/outlier phase
inequalities, and all `1,756,510` odd owner classes at the manufactured
clocks. It separately checks odd/even dilations, six hundred `H_8+7m` gates,
and `(40)`.

The independent companion never uses the parity-compressed teeth in `(3)`.
It constructs both literal original-theta walls

```text
theta=(k +/- 1/14)/d-h/2 mod 1,       h in {0,1},           (41)
```

tests every wall and intervening open cell, and joins adjacent bad cells only
when their common wall is itself strictly bad. It checks `1,225` component
pairs, `3,240` family rows, `71,280` literal family phase inequalities, and
`2,004,002` odd owner classes. It reconstructs all twelve isolated literal
`H_8` phases, zero safe open cells, the MISTAKE-274 handoff, and the five
literal cells in `(40)`.

Replay from the repository root:

```bash
python3 -B 04-computation/lrc_literal_open_two_comb_compiler_thm4087.py
python3 -B -O 04-computation/lrc_literal_open_two_comb_compiler_thm4087.py
python3 -B 04-computation/lrc_literal_open_two_comb_compiler_thm4087_independent_audit.py
python3 -B -O 04-computation/lrc_literal_open_two_comb_compiler_thm4087_independent_audit.py
```

The source object is a metric antipodal-safe interval. The map removes two
literal open danger combs and selects a surviving arrangement endpoint. It
preserves pointwise weak safety, both half labels, and the owner eligibility
predicate. Passing only to interval length forgets location and endpoint
owners; the rational-endpoint sidecar restores an adaptive labelled clock.
The cheapest decisive hostile is the zero-width set `(36)--(38)`.

This theorem extends the all-height structural surface beyond THM-4075's
bounded cores and THM-4079's one-outlier family. It closes exactly the
two-outlier AP9 dyadic seams `(31)`, including nonprimitive dilations. It does
not show that an arbitrary THM-2061 residual contains such an interval or AP9
body, does not handle three outliers (whose components can chain through more
than one lower-speed tooth), does not decide `kappa=7,8`, does not provide the
missing physical-entry map, and does not prove LRC(14). **QED.**
