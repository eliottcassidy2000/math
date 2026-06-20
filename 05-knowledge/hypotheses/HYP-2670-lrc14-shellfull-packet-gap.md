---
id: HYP-2670
title: LRC14 shell-full p1-tax packet gap beyond the B13 leader
status: OPEN; exact B=30 shell-full quotient evidence
source: codex-2026-06-20-S44
depends_on:
  - HYP-2669
  - HYP-2668
  - HYP-2667
  - HYP-2666
  - HYP-2661
related:
  - HYP-2665
  - HYP-2664
  - HYP-2662
  - HYP-2660
  - HYP-2643
  - OPEN-Q-108
---

# HYP-2670 - Shell-Full Packet Gap

## Claim

After imposing the shell-1-full quotient `{1,2,4,8} subset E'`, the p1-tax
frontier appears to split into a finite B13 packet pocket plus a genuine
new-speed decay lemma.

The tested exact quotient through `B=30`,

```text
E'={0}+{1,2,4,8}+3 extras from [1,30],
w=max(E')+1,...,max(E')+8,
```

has only one row above `1/3`, the old S41/S43 leader

```text
E'=(0,1,2,4,6,7,8,10), w=12,
Delta_w^+/p1 = 997/2562.
```

Every row introducing a speed beyond `14` stays below `1/3`, and every row
introducing a speed beyond `24` stays below `1/4` in the scan.

Thus the shell-full side of the two-gate proof should not be attacked as a
single growing scalar frontier.  The sharper target is:

```text
finite packet ledger for max(E') <= 14,
new-speed decay for max(E') > 14,
far-tail packet decay for max(E') > 24.
```

## Evidence

Script:

```text
04-computation/lrc14_shellfull_packet_gap_codex_s44.py
```

Stored output:

```text
05-knowledge/results/lrc14_shellfull_packet_gap_codex_s44.out
```

Exact scan summary:

```text
rows=20800
global shell-full max = 997/2562
new-speed max max(E')>14 = 1371/4319
tail max max(E')>24 = 932669/4085893
```

Layer table:

```text
finite <=14: rows=960,   >1/3=1, max=997/2562
new 15..24: rows=8160,   >1/3=0, max=1371/4319
tail 25..30: rows=11680, >1/3=0, max=932669/4085893
```

Exact gaps:

```text
1/3 - 1371/4319       = 206/12957
1/4 - 932669/4085893  = 355217/16343572
```

The B13 leader still has exact `2p1/5` tax gap

```text
139/2450.
```

## Packet Anatomy

The leader is not just high because it has shell 1.  It has a concentrated
small-denominator phase packet:

```text
positive_packets/(w*p1)=3529/5978
negative_packets/(w*p1)=1804/8967
small-denominator positive share, denom <= 35: 9527/10587
top positive share: 350/3529
positive denominators: 5, 7, 14, 35, 49
```

Its fold profile also keeps more high-target fold mass than the new-speed
leader:

```text
leader fold_recip = 319/420
new-speed leader fold_recip = 59/240
```

The new-speed leader

```text
E'=(0,1,2,4,8,12,16,20), w=24
```

keeps the small-denominator phase but loses enough fold/tower concentration to
fall to `1371/4319 < 1/3`.  This is the first useful explanation of HYP-2669's
stability: extending the dyadic tower creates packets, but it also dilutes the
fold target currency that made the B13 pocket dangerous.

## Proof Route

The current shell-full route should become three lemmas.

1. Finite pocket lemma: enumerate or symbolically certify the `max(E')<=14`
   shell-full packet ledger, with the B13 leader as the unique row above `1/3`
   and still below `2/5`.
2. New-speed lemma: if shell-full and `max(E')>14`, prove
   `Delta_w^+ <= p1(E')/3`.  In the B30 scan the exact margin is `206/12957`.
3. Far-tail lemma: if shell-full and `max(E')>24`, prove a stronger decay,
   plausibly `Delta_w^+ <= p1(E')/4` or a nearby packet-dependent version.

Then HYP-2666/HYP-2668 can assemble:

```text
shell-damaged rows -> shell/mouth gate,
shell-full finite pocket -> exact packet ledger,
shell-full new speeds -> p1/3 or p1/4 decay,
all shell-full rows -> 2p1/5.
```

Incoming KPS S17 work after this scan strengthens the first line of that
assembly: THM-545 proves the one-extra-hole tower-deletion layer, and the
wide k=2 tower-deletion scans report `0` sub-threshold rows through the
binding bit-4 window.  Thus the shell-damaged gate is no longer just a useful
stratification; it is close to a rigorous discharge.  HYP-2670 is the matching
post-gate obligation.

## Tournament Analysis

Vertices:

```text
finite_B13_pocket
new_speed_gap
phase_denominator_cliff
fold_target_concentration
raw_runner_vertices
```

Pairwise observable: exact `Delta_w^+/p1` after quotienting by the full shell-1
tower.

Switch/gauge: compare proof-obligation layers by `max(E')`, then expose phase
packets and fold targets.

Hamiltonian path:

```text
finite_B13_pocket
> new_speed_gap
> phase_denominator_cliff
> fold_target_concentration
> raw_runner_vertices
```

Challenged assumption: shell-full proof work is not a single growing scalar
frontier.  The high packet appears finite, while new speeds have a visible
decay gap.

Alternate vertices considered: runners, gaps, fixed sectors, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, fold
targets, Glaisher odd-carry words, and proof-obligation layers.

Preserved predicate: the quotient keeps the p1-tax inequality
`Delta^+ <= c*p1`.

Destroyed information: coarse interval-envelope placement outside the
phase-packet address.

## Honest Status

LRC(14) is not proved.  HYP-2670 is bounded exact evidence that the shell-full
half of the current proof has a finite high packet plus a new-speed decay
structure, rather than a moving scalar frontier.
