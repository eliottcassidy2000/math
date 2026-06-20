---
id: HYP-2664
title: LRC14 residual/plateau packet tax - far residual is paid by p1 boundary mass
status: OPEN; exact finite evidence and new proof target
source: codex-2026-06-19-S39
depends_on:
  - HYP-2662
  - HYP-2655
  - HYP-2653
  - HYP-2661
  - HYP-2663
related:
  - HYP-2648
  - HYP-2658
  - HYP-2660
  - OPEN-Q-108
---

# HYP-2664 - Residual/Plateau Packet Tax

## Claim

The HYP-2662 intra-quadratic residual should be treated as a boundary tax on
the single-missed-sector chain, not as a free endpoint-count error.

For `E=E' union {w}`,

```text
p0(E) = Phi(E') + Delta_w,
Phi(E') = p0(E') + p1(E')/7.
```

HYP-2662 decomposes `w*Delta_w` into trace, QR/NQR, and intra-quadratic
residual channels.  HYP-2664 proposes the next sharp inequality:

```text
positive far discrepancy Delta_w is controlled by c*p1(E')
```

with a proof-friendly target `c=1/3`, and a sharper residual-only target near
`c=1/4`.  The key point is that `p1(E')` is already present in the plateau, so
this turns the resonant endpoint problem into a joint plateau/Delta inequality.

## Evidence

Script:

```text
04-computation/lrc14_residual_plateau_packet_codex_s39.py
```

Stored output:

```text
05-knowledge/results/lrc14_residual_plateau_packet_codex_s39.out
```

The script aggregates endpoint terms by the actual phase packet
`y=frac(w*x)`, quotients each packet by QR/NQR class means, and measures
intra-class imbalance before scalarizing.

On the known HYP-2662 resonant families:

```text
max |intra_quadratic_residual|/(w*p1) = 29166/144607 ~= 0.201691
max |raw wDelta|/(w*p1)               = 2804017/11198173 ~= 0.250400
min CAP9 - Phi - |residual|/w         = 77424533/613602990 ~= 0.126180
```

So the residual channel stays well below `p1/4` in the tested rows.  The raw
far discrepancy slightly exceeds `p1/4` in one multiscale row, so the robust
proof target should be `raw positive Delta <= p1/3`, or a split proof with
residual `<= p1/4` plus a smaller trace/QR tax.

Bounded AP-window bank:

```text
E' = {0} + 7-subsets of [1,13], 1716 rows.
cap9 = 1979/4004.
minimum allowed c in p0 + (1/7+c)*p1 <= cap9:
  388929/718718 ~= 0.541143 at AP8.
```

In particular:

```text
c=1/4: 0 violations, max = 17417/41160 at AP8, slack = 418499/5885880.
c=13/51: 0 violations, max = 44539/104958 at AP8, slack = 2098409/30017988.
c=1/3: 0 violations, max = 13691/30870 at AP8, slack = 448069/8828820.
```

Thus a theorem of the form

```text
Delta_w <= p1(E')/3
```

would clear the entire bounded AP-window 9-row bank with exact slack, while the
wide/multiscale rows already have much smaller `Phi(E')` and large cap margin.

## Interpretation

This is the outside-the-box bridge between three recent ideas:

1. HYP-2662: the Galois phase does not collapse to trace/QR; the residual
   lives inside QR/NQR classes.
2. HYP-2663: the near-collar proof object is a root packet, not a replacement
   count.
3. HYP-2661: the bounded mouth side is protected by the dyadic-1 tower clamp.

The common pattern is packet closure before scalarization.  Far endpoints
should be grouped by `y=frac(w*x)`, then class-balanced inside QR/NQR, and only
then converted into a scalar discrepancy.  The packet defect is a boundary of
the one-missed-sector region, so its natural price is `p1(E')`.

## Proof Target

Prove one of the following:

```text
(A) Delta_w^+ <= p1(E')/3.
(B) residual^+/w <= p1(E')/4 and trace+QR positive tax <= p1(E')/12.
(C) A sharper packet coarea inequality using phase_packet_class_l1 and p1.
```

Then combine with the exact bounded-bank certificate:

```text
p0(E') + (1/7+1/3)*p1(E') <= cap9
```

on `E'={0}+7`-subsets of `[1,13]`.  Large-tail and nonlocal rows should route
through HYP-2655/HYP-2658 recursion, while near-collar small-plateau rows route
through HYP-2661/HYP-2663 packet rigidity.

## Tournament Analysis

Vertices are proof carriers:

```text
p1_boundary_mass
phase_packet_class_l1
residual_channel
raw_wDelta
raw_endpoint_count
```

Pairwise observable: which carrier predicts `Delta_w` without losing the
plateau slack.

Switch/gauge: aggregate endpoints by `y=frac(w*x)`, then quotient by QR/NQR
class means.

Hamiltonian path:

```text
p1_boundary_mass
> phase_packet_class_l1
> residual_channel
> raw_wDelta
> raw_endpoint_count
```

Challenged assumption: the far residual is not controlled by endpoint count or
by the Galois trace alone.  It is a packet-boundary imbalance, and the missing
mass `p1(E')` is the proof-facing currency.

## Honest Status

LRC(14) is not proved.  HYP-2664 supplies a much sharper target than the false
uniform `C` route: prove a `p1`-tax inequality for far discrepancy, then use
the bounded AP-window certificate and the HYP-2661/HYP-2663 packet rigidity.
