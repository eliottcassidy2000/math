---
id: HYP-2665
title: LRC14 p1-tax envelope refinement - p1/3 is false; 3p1/8 remains viable
status: OPEN; refines HYP-2664 by refuting the p1/3 constant
source: codex-2026-06-19-S40
depends_on:
  - HYP-2664
  - HYP-2662
  - HYP-2655
related:
  - HYP-2661
  - HYP-2663
  - HYP-2648
  - HYP-2658
  - OPEN-Q-108
---

# HYP-2665 - p1-Tax Envelope Refinement

## Claim

The HYP-2664 idea is right at the level of proof currency (`p1(E')` is the
natural boundary mass), but the clean constant

```text
Delta_w^+ <= p1(E')/3
```

is false.  Exact stress data finds bounded rows with `Delta_w^+/p1 > 1/3`.
The next viable target is instead

```text
Delta_w^+ <= 3*p1(E')/8
```

or a packet-refined inequality that treats the few `>1/3` rows by their
phase-packet structure.

## Evidence

Script:

```text
04-computation/lrc14_p1_tax_envelope_codex_s40.py
```

Stored output:

```text
05-knowledge/results/lrc14_p1_tax_envelope_codex_s40.out
```

The script scans structured resonant families plus targeted bounded banks and
compares the actual ratio `Delta_w^+/p1` to an interval-length envelope.

Exact refutations of `1/3`:

```text
structured odd row:
E'=(0,1,3,5,7,9,10,11), w=13
Delta_w^+/p1 = 95943/269360 ~= 0.356189
```

```text
targeted B=13 row:
E'=(0,1,2,3,6,7,10,11), w=12
Delta_w^+/p1 = 4691/13076 ~= 0.358749
```

Two additional sampled B=13 rows also exceed `1/3`:

```text
(0,1,2,3,5,7,9,11), w=13: 320249/940758 ~= 0.340416
(0,1,2,3,6,7,8,11), w=12: 6493/19180 ~= 0.338530
```

All observed ratios remain below `3/8`.

## Interval Envelope Is Too Coarse

The elementary interval envelope asks only how long the single-missed cells are
and ignores phase cancellation.  It fails badly:

```text
max envelope/p1 = 6/7
```

In the targeted B=13 sample, every row has envelope/p1 above `1/3`, even though
only three sampled rows have actual `Delta_w^+/p1 > 1/3`.  Thus a direct
length/coarea proof cannot close the p1 tax.  The missing ingredient is
phase-packet cancellation: group endpoints by `y=frac(w*x)`, then use QR/NQR
class and packet signs before scalarizing.

## Why 3/8 Is Plausible

HYP-2664's bounded AP-window bank already computed the exact tolerance:

```text
min allowed c in p0 + (1/7+c)*p1 <= cap9
  = 388929/718718 ~= 0.541143.
```

So any theorem `Delta_w^+ <= c*p1(E')` with `c <= 3/8` clears that bounded
bank with large slack.  The S40 stress sample finds no ratio above `3/8`, with
global maximum

```text
4691/13076 ~= 0.358749 < 3/8.
```

The updated sharp target is therefore:

```text
prove Delta_w^+ <= 3*p1(E')/8,
or classify the few rows above 1/3 and prove a packet cancellation bound there.
```

## Tournament Analysis

Vertices are proof carriers:

```text
actual_p1_tax
interval_envelope
phase_packet_cancellation
endpoint_count
raw_speed
```

Pairwise observable: which carrier certifies the positive far discrepancy.

Switch/gauge: replace endpoint count by interval-length envelope, then refine
the envelope by phase packets.

Hamiltonian path:

```text
actual_p1_tax
> phase_packet_cancellation
> interval_envelope
> endpoint_count
> raw_speed
```

Challenged assumption: `p1/3` looked attractive because the bounded AP-window
bank had huge slack, but the far discrepancy itself can exceed `p1/3`.  The
right target must leave room for local packet resonance while still charging
the error to `p1`.

## Honest Status

LRC(14) is not proved.  HYP-2665 refines HYP-2664: the p1 currency remains
promising, the constant `1/3` is false, and the next plausible theorem is a
`3p1/8` tax plus packet cancellation for the near-extremal rows.
