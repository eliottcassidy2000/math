---
id: HYP-2663
title: LRC14 AP-tail root-packet scout - three-tail packet clears the AP-second threshold to 35
status: OPEN; exact finite packet certificate
source: codex-2026-06-19-S38
depends_on:
  - HYP-2654
  - HYP-2659
  - HYP-2660
  - HYP-2661
  - HYP-2648
related:
  - HYP-2537
  - HYP-2658
  - HYP-2655
  - THM-541
  - THM-543
  - THM-544
  - OPEN-Q-108
---

# HYP-2663 - LRC14 AP-Tail Root-Packet Scout

## Claim

The old FKN/root-packet lesson (HYP-2537) applies directly to the LRC14
near-collar AP-tail problem: perturbations of the AP collar should be treated as
addressed packets of holes and tail insertions, not as uniform replacement
counts.  KPS HYP-2661 supplies the sharper theorem-shaped carrier: sub-second
AP-tail rows should conserve the full shell-1 dyadic tower `{1,2,4,8}` with
odd-shell carry `{1:15}`.  This scout is the next finite evidence layer for
that conservation law.

For the first AP-tail layer beyond THM-544,

```text
C = ({1,...,13} \ H) union {r,s,t},
|H|=4, 14 <= r < s < t,
```

the exact scout through `t <= 35` finds no packet below the AP one-hole second
threshold `426/35035`.  The best scanned packet is already

```text
H=(3,4,10,12), tails=(17,19,20),
meas(G_C)=4309/255255
4309/255255 - 426/35035 = 59/12495.
```

Thus the three-tail AP packet has a large finite safety margin in the same
window where the older AP-tail scout operated.

## Computation

Script:

```text
04-computation/lrc14_root_packet_ap_tail_scout_codex_s38.py
```

Stored output:

```text
05-knowledge/results/lrc14_root_packet_ap_tail_scout_codex_s38.out
```

The script reuses the sorted-arc exact interval engine from THM-544 and adds:

1. Glaisher odd-shell carry delta relative to the drop-6 collar.
2. Survival/damage of the four drop-6 mouth intervals.
3. Comb lower-bound pruning for tails already forced above `426/35035`.

Bound-35 scan:

```text
base packets: 715
r rows:       15730
s rows:       165165
exact rows:   1076482
one-tail prunes: 24618
below threshold: 0
```

The top packets all either have zero old-mouth survivor or only partial survivor
(`1/364` or smaller), so the three-tail packet is no longer in the retained
drop-6 mouth regime.  It pays above the AP-second threshold before scalar
measure becomes dangerous.  The best row damages shell-1 by `{1:-4}`, matching
HYP-2661's carry-conservation prediction that spending a dyadic-1 tower bit
pushes the row above the AP-second threshold.

## Hidden-Gem Transfer

HYP-2537 said that near a transitive tournament, the useful local coordinates
are interval-root packets, not raw Hamming shells.  The LRC translation is:

```text
replacement count              -> too coarse
hole/tail root packet           -> useful local carrier
Glaisher carry delta            -> dyadic tower address
drop-6 mouth survival/damage    -> fixed-observer geometry
```

This folds the older FKN/root-packet work into the current HYP-2654/HYP-2659
near-collar route.

## Proof Route Impact

The local AP-tail stack now reads:

```text
THM-541: one-hole AP collar has unique minimum drop 6.
THM-543: one-tail layer has one below-second exception, (6,10)->20,
         and it fully retains the drop-6 mouth.
THM-544: two-tail layer has no below-second row.
HYP-2661: shell-1 carry conservation predicts sub-second rows keep `{1:15}`.
HYP-2663: three-tail packet through tail bound 35 has no below-second row,
          with best margin 59/12495.
```

This suggests the next theorem should not be a generic replacement-count
induction.  It should be a **shell-1/mouth-damage packet theorem**:

```text
If a root packet damages the shell-1 tower or the four drop-6 mouths, or opens
genuinely new odd-shell carry, then it pays at least 426/35035.
```

Packets preserving the drop-6 mouth are finite templates feeding
HYP-2654/HYP-2659/HYP-2660.  Packets with large tails or no local mouth template
should route through HYP-2655/HYP-2658 far/core-gap recursion.

## Tournament Analysis

Vertices are proof carriers:

```text
root_packet
mouth_survival
odd_shell_carry
comb_tail_bound
raw_speed
```

Pairwise observable: which carrier preserves the predicate
`meas(G_C)<426/35035` before scalarization.

Hamiltonian path:

```text
root_packet > mouth_survival > odd_shell_carry > comb_tail_bound > raw_speed
```

Challenged assumption: AP-tail vertices are not replacement counts or raw
speeds.  They are addressed root packets, exactly as the old tournament FKN
root-shell work warned.

## Next Work

1. Promote the bound-35 packet certificate into a theorem-style finite proof
   with explicit comb cutoff for all `t>35`.
2. Search whether every packet with full drop-6 mouth survival is one of the
   already-known templates.
3. Add an interval-root incidence graph on the AP holes and tail owners, then
   prove a packet-energy lower bound analogous to the HYP-2537 shell-2 law.
4. Use HYP-2661's shell-1 carry conservation to separate harmless
   existing-shell carry from genuinely new-shell debt.
