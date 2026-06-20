---
id: HYP-2664
title: LRC14 residual/plateau packet tax - far residual is paid by p1 boundary mass
status: OPEN; p1 currency supported, but constants require shell-gated packets
source: codex-2026-06-19-S39
depends_on:
  - HYP-2662
  - HYP-2655
  - HYP-2653
  - HYP-2661
  - HYP-2663
related:
  - HYP-2665
  - HYP-2668
  - HYP-2667
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

with an initial proof-friendly target `c=1/3`, and a sharper residual-only
target near `c=1/4`.  The key point is that `p1(E')` is already present in the plateau, so
this turns the resonant endpoint problem into a joint plateau/Delta inequality.

**S40/S41/S42 correction:** HYP-2665 refutes the raw `c=1/3` constant,
HYP-2667 refutes the raw `c=3/8` constant on the full bounded `B=13` bank,
and HYP-2668 refutes global `c=2/5` on the full bounded `B=14` bank.  The
`p1` currency remains useful, but only after the right packet gate: the current
target is `c=2/5` on the shell-1-full quotient, with shell-damaged rows routed
through HYP-2661/HYP-2666 first.

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
far discrepancy slightly exceeds `p1/4` in one multiscale row, later S40/S41
bounded rows exceed both `1/3` and `3/8`, and S42 has one shell-damaged B=14
row above `2/5`.  The robust proof target should now be shell-gated: prove raw
positive Delta `<= 2*p1/5` on the shell-1-full quotient, while shell-damaged
rows pay via tower/mouth rigidity.

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
HYP-2665 shows this theorem is too strong as a raw far-discrepancy claim:
`Delta_w^+/p1` reaches `4691/13076 ~= 0.358749` in a targeted bounded row.
HYP-2667 further shows `3/8` is too strong: the full B=13 bank reaches
`997/2562 ~= 0.389149` at `E'=(0,1,2,4,6,7,8,10), w=12`, while `2/5`
survives the full bank with minimum tax gap `139/2450`.
HYP-2668 then shows global `2/5` is too strong at B=14:
`7071/17584 ~= 0.402127` for `E'=(0,1,4,6,8,10,12,14), w=16`, but that row
misses shell-1 bit `2`; all shell-1-full B=14 rows remain below `2/5`.

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
(A') If E' is shell-1-full, Delta_w^+ <= 2*p1(E')/5.
(A'') If E' damages shell 1, route through HYP-2661/HYP-2666 tower/mouth tax first.
(B) residual^+/w <= p1(E')/4 and trace+QR positive tax <= p1(E')/12.
(C) A sharper packet coarea inequality using phase_packet_class_l1 and p1.
```

Then combine with the exact bounded-bank certificate:

```text
p0(E') + (1/7+2/5)*p1(E') <= cap9
```

on the shell-1-full bounded quotient; HYP-2668 is the current finite evidence
for this gate.  Large-tail and nonlocal rows should route
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
uniform `C` route, but HYP-2665/HYP-2667/HYP-2668 correct the constants:
prove a shell-gated `2p1/5` raw tax, or a packet-refined `p1` tax with
dyadic/shell gates, then use the bounded AP-window certificate and the
HYP-2661/HYP-2663 packet rigidity.
