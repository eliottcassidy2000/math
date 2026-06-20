---
id: HYP-2662
title: LRC14 far Delta Galois phase - trace and QR/NQR projections leave an intra-quadratic residual
status: OPEN; exact decomposition, simple QR/Galois bound refuted
source: codex-2026-06-19-S38
depends_on:
  - HYP-2653
  - HYP-2655
  - HYP-2657
  - HYP-2660
related:
  - HYP-2661
  - HYP-2644
  - HYP-2646
  - HYP-2648
  - HYP-2658
  - OPEN-Q-108
---

# HYP-2662 - Far Delta Galois Phase

## Claim

The far-element discrepancy

```text
w*Delta_w = sum_cells [G0(w*b - s/7) - G0(w*a - s/7)]
```

has an exact endpoint-level apex-prime decomposition:

```text
G0(y - s/7)
  = F7_trace_average(y)
    + chi_7(s) * QR_NQR_quadratic_channel(y)
    + intra_quadratic_residual(y,s).
```

The first term is the full `F_7^*` orbit average over nonzero sector phases.
The second is the quadratic-character projection separating QR sectors
`{1,2,4}` from NQR sectors `{3,5,6}`.  The third has zero average inside each
quadratic class.

This decomposition is exact and proof-facing, but the simple hope

```text
large resonant Delta is controlled by the trace or QR/NQR projection
```

is false on the tested HYP-2655 multiscale rows.  The dominant channel is the
intra-quadratic residual, i.e. unequal endpoint/breakpoint weights *within* a
quadratic class.  Thus the Galois phase remains relevant, but the proof must
bound the residual endpoint-weight imbalance jointly with the plateau margin.

## Evidence

Script:

```text
04-computation/lrc14_far_delta_galois_phase_codex_s38.py
```

Stored output:

```text
05-knowledge/results/lrc14_far_delta_galois_phase_codex_s38.out
```

The script verifies the decomposition identities through denominator `29` and
then scans named resonant cores with exact `Fraction` arithmetic.

Important exact rows:

```text
odd_struct_HYP2655, E'=(0,1,3,5,7,9,10,11), w=90
raw wDelta = 1370/539
trace      = 167/3234
quadratic  = 1735/3234
residual   = 1053/539
```

```text
multiscale_30_60, E'=(0,1,2,30,31,32,60,61), w=62
raw wDelta = 2804017/717360
trace      = 201751/2152080
quadratic  = 9073/10980
residual   = 114857/38430
```

```text
cluster_40_80, E'=(0,1,2,40,41,42,80,81), w=82
raw wDelta = 2697113/555660
trace      = 38267/246960
quadratic  = 964849/952560
residual   = 125399/34020
```

The scale-family trend is monotone in the tested rows and residual-dominant:

```text
M=10 raw=12557/10290, residual=216/245
M=20 raw=53197/20090, residual=90287/47355
M=30 raw=2804017/717360, residual=114857/38430
M=40 raw=2697113/555660, residual=125399/34020
M=50 raw=19373093/3216850, residual=6305213/1378650
```

## Interpretation

HYP-2657's Galois trace result is still load-bearing: it supplies an exact
phase projection and explains why the correction is real after full
multiplicative averaging.  But HYP-2662 shows that the far-element breakpoint
weights are not close enough to `F_7^*`-uniform, or even QR/NQR-uniform, for
that projection alone to bound resonance.

The incoming KPS HYP-2661 shell-1 carry-conservation law is the matching
plateau-side signal: the bounded near-collar rows below `426/35035` must keep
the full dyadic-1 tower `{1,2,4,8}`, while damaging that tower pays the AP
second threshold.  Read together, HYP-2661 constrains the tight mouth/plateau
templates and HYP-2662 isolates the far-resonant residual channel that remains
to be bounded against the shrinking plateau.

The right replacement target is:

```text
trace channel
+ QR/NQR channel
+ intra-quadratic endpoint-weight residual
+ plateau margin
< cap_k.
```

Equivalently, prove that the same multiscale structure that increases the
intra-quadratic residual also drives the plateau
`p0(E') + p1(E')/7` far below the AP cap.  This is the HYP-2655 joint
plateau/Delta recursion with a sharper Galois-phase ledger.

## Tournament Analysis

Vertices are phase/proof channels:

```text
F7_trace_average
QR_NQR_quadratic_channel
intra_quadratic_residual
raw_wDelta
sigma_bound
plateau_margin
```

Pairwise observable: which channel carries the exact resonant value of
`|w*Delta_w|`.

Switch/gauge: multiply missed-sector labels by `F_7^*` before scalarizing.

Tie Hamiltonian path:

```text
intra_quadratic_residual
> plateau_margin
> QR_NQR_quadratic_channel
> F7_trace_average
> raw_wDelta
> sigma_bound
```

Fingerprint: transitive residual-dominance proof-channel order; directed
`3`-cycles `0`.

Assumption challenge: the Galois phase is not a single scalar residue class or
runner-level tournament.  The useful vertices are endpoint-weight channels
after the `F_7^*` action; collapsing to QR/NQR too early destroys the dominant
residual that actually carries multiscale resonance.
