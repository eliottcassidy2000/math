---
id: HYP-2667
title: LRC14 p1-tax packet frontier - 3p1/8 is false; B13 survives 2p1/5
status: OPEN; full B=13 exact frontier and packet anatomy
source: codex-2026-06-19-S41
depends_on:
  - HYP-2665
  - HYP-2664
  - HYP-2662
  - HYP-2655
related:
  - HYP-2666
  - HYP-2661
  - HYP-2663
  - HYP-2658
  - OPEN-Q-108
---

# HYP-2667 - p1-Tax Packet Frontier

## Claim

The `p1(E')` boundary currency remains the right way to price positive far
discrepancy, but the scalar constant must be revised again.  The full bounded
`B=13` bank refutes the S40 provisional target

```text
Delta_w^+ <= 3*p1(E')/8.
```

In that exact bank, every tested row still satisfies

```text
Delta_w^+ <= 2*p1(E')/5.
```

The two rows above `3/8` are not broad noise; they are low bounded `w=12`
rows with strong even/dyadic phase packets.  The live theorem should be either
a raw `2p1/5` tax, or a split theorem:

```text
generic phase packets <= 3*p1/8,
dyadic-even packet frontier <= 2*p1/5.
```

## Evidence

Script:

```text
04-computation/lrc14_p1_tax_packet_frontier_codex_s41.py
```

Stored output:

```text
05-knowledge/results/lrc14_p1_tax_packet_frontier_codex_s41.out
```

The full exact bank is

```text
E'={0}+7-subsets of [1,13],
w=max(E')+1,...,max(E')+8.
```

It contains `13728` primitive rows.  Exact threshold counts:

```text
Delta_w^+/p1 > 1/3: 17 rows
Delta_w^+/p1 > 3/8: 2 rows
Delta_w^+/p1 > 2/5: 0 rows
```

The worst two rows are:

```text
E'=(0,1,2,4,6,7,8,10), w=12:
  Delta_w^+/p1 = 997/2562 ~= 0.389149
  gap below 2/5 = 139/2450
  cap slack for c=2/5 = 84676/525525

E'=(0,2,4,6,7,8,9,10), w=12:
  Delta_w^+/p1 = 5347/14070 ~= 0.380028
  gap below 2/5 = 562/5145
  cap slack for c=2/5 = 1418467/8828820
```

Across the full bank:

```text
max Delta_w^+/p1 = 997/2562
min 2/5 tax gap c*w*p1 - raw = 139/2450
min cap slack for c=2/5 = 507209/14714700
```

Thus `2/5` clears the full bounded B=13 p1-tax bank with exact slack, while
`3/8` does not.

## Packet Anatomy

The worst row has full shell-1 tower `{1,2,4,8}` inside `E'`, holes
`(3,5,9,11,12,13)`, and `w=12`.  Its endpoint phase packets satisfy:

```text
positive_packets/(w*p1) = 3529/5978 ~= 0.590331
negative_packets/(w*p1) = 1804/8967 ~= 0.201182
top_positive_share      = 350/3529 ~= 0.099178
```

So the excess is not one giant endpoint.  It is a medium-sized packet stack.
The largest positive packets are concentrated at rational phases with small
denominators:

```text
y=6/7:  raw/(w*p1)=25/427,  raw=15/49
y=3/7:  raw/(w*p1)=10/183,  raw=2/7
y=2/7:  raw/(w*p1)=20/427,  raw=12/49
y=31/35,3/5,11/35: each raw/(w*p1)=58/1281
```

The second `3/8` failure also has `w=12`, mostly even/dyadic structure, and a
largest positive packet at `y=2/7` with `raw=19/49`.

This is the new motif: the interval envelope is too coarse, but packet
cancellation is also not uniformly strong enough to force `3/8`; the dyadic
bounded row can align many medium phase packets.

## Namespace Note

Concurrent mainline work currently contains two `HYP-2664` detail files:
the residual/plateau p1-tax thread and the three-tail shell-1 frontier
`T907`.  Concurrent work also claimed `HYP-2666` for the two-gate boundary
currency scout.  This HYP is therefore numbered `HYP-2667`; it depends on the
residual/plateau p1-tax thread and corrects the raw-tax constant in the
two-gate HYP-2666 route.

## Tournament Analysis

Vertices:

```text
dyadic_even_packet
phase_packet_concentration
p1_tax_constant
interval_envelope
endpoint_count
```

Pairwise observable: which carrier explains rows with `Delta_w^+/p1 > 3/8`.

Switch/gauge: first sort by exact p1 ratio, then expose endpoint packets
`y=frac(w*x)` and their QR/NQR residual ledgers.

Hamiltonian path:

```text
dyadic_even_packet
> phase_packet_concentration
> p1_tax_constant
> interval_envelope
> endpoint_count
```

Challenged assumption: after `p1/3` failed, `3p1/8` looked viable in the S40
sample.  The full B=13 bank shows a small dyadic-even frontier above `3/8`, so
constants must be stress-tested before canonization.

## Honest Status

LRC(14) is not proved.  HYP-2667 gives a sharper next obligation: prove a
`2p1/5` raw far-discrepancy tax, or split the proof into a generic `3p1/8`
packet theorem plus a finite dyadic-even packet frontier.
