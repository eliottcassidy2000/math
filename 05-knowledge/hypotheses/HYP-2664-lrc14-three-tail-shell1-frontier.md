---
id: HYP-2664
title: LRC14 three-tail shell-1 frontier - carry gate before root-packet comb
status: OPEN; exact cutoff-frontier atlas
source: codex-2026-06-19-S39
tangent: T907
depends_on:
  - HYP-2661
  - HYP-2663
  - HYP-2654
  - HYP-2659
  - HYP-2660
related:
  - HYP-2658
  - HYP-2655
  - THM-541
  - THM-543
  - THM-544
  - OPEN-Q-108
---

# HYP-2664 - LRC14 Three-Tail Shell-1 Frontier

## Claim

The three-tail AP-tail proof should apply KPS HYP-2661's shell-1 carry conservation before doing root-packet or comb enumeration.

HYP-2663 scanned the first three-tail AP-tail layer

```text
C = ({1,...,13} \ H) union {r,s,t},
|H|=4, 14 <= r < s < t <= 35,
```

and found no rows below the AP one-hole second threshold `426/35035`. A naive unbounded nested comb proof leaves a large exact finite residue set. The S39 frontier atlas shows why: the largest crude three-tail comb burdens are mostly packets that delete a shell-1 tower bit from `{1,2,4,8}`, exactly the packets HYP-2661 says should already pay the second threshold.

## Computation

Scripts:

```text
04-computation/lrc14_three_tail_shell1_frontier_codex_s39.py
04-computation/lrc14_three_tail_ap_tail_comb_codex_s39.py
```

Stored result:

```text
05-knowledge/results/lrc14_three_tail_shell1_frontier_codex_s39.out
```

The frontier script computes, for every four-hole AP base, the first-tail cutoff `R` from the crude three-tail comb bound:

```text
meas(G_base with three tails >= R)
  >= (4/7) meas(G_base) - 6*c_base/(7R).
```

This is not a full proof. It is a residue-burden atlas for the next proof.

## Exact Frontier Data

There are `715` four-hole AP bases. Of these, `589` damage shell 1 and `126` preserve the full tower `{1,2,4,8}`.

```text
top-40 crude comb bases damaged by HYP-2661 gate: 37/40
top-100 crude comb bases damaged by HYP-2661 gate: 87/100
```

The first-tail cutoff quantiles:

```text
all bases:       min 69, q25 124, median 151, q75 181, q90 212, q99 280, max 308
shell1 damaged:  min 69, q25 123, median 151, q75 181, q90 216, q99 284, max 308
shell1 full:     min 71, q25 125, median 150, q75 179, q90 199, q99 232, max 239
```

The global worst base:

```text
holes=(4,5,6,13), R=308, shell1 damaged by missing 4.
```

Worst shell-1-preserving base:

```text
holes=(3,5,6,13), R=239, shell1 full.
```

So shell-1 gate drops worst crude first-tail cutoff from `308` to `239`, and removes almost all top frontier. In crude tail-triple load, global top is about `4.19M` triples below first cutoff; shell-1-full top is about `1.87M`.

## Proof Route Impact

The failed broad-comb attempt was informative. Nested cutoffs are finite, but exact third-tail interval subtraction is too expensive when run before using the carry law. Proof order should be:

```text
1. shell-1 carry gate: if a packet deletes 1,2,4,or 8, use HYP-2661;
2. root-packet address: among shell-1-full packets, retain holes/tails/carry;
3. mouth-owner ledger: separate full old-mouth survival from damage;
4. nested comb finite residues only after those quotients.
```

This reframes HYP-2663 next theorem target:

```text
If a three-tail AP packet is below 426/35035, then it must be shell-1 full
and must lie in a small old-mouth/root-packet template family.
```

Then shell-1-full frontier, with max crude cutoff `239`, is the finite residue region to attack.

## Tournament Analysis

Vertices:

```text
shell1_gate
root_packet
mouth_owner
nested_comb
raw_replacement_count
```

Observable: finite residue burden for proving `meas(G_C)>=426/35035`. Orient `A -> B` when applying `A` first strictly reduces or structures finite obligations left to `B`.

Hamiltonian path:

```text
shell1_gate > root_packet > mouth_owner > nested_comb > raw_replacement_count
```

Alternate vertex sets considered: raw speeds, replacement count, AP holes, tail values, Glaisher odd shells, drop-6 mouth intervals, endpoint-owner pairs, comb cutoffs, and proof obligations. Best quotient uses shell-1 carrier gate plus addressed root packet.

Challenged assumption: first AP-tail quotient is not number of tails or maximum tail. It is whether shell-1 carrier survives.

## Next Work

1. Prove HYP-2661 shell-1 deletion theorem independently of bounded scans.
2. Run exact nested comb only on shell-1-full packets, with mouth-owner pruning.
3. Classify shell-1-full three-tail packets by old-mouth survivor value.
4. Look for direct inequality on remaining worst base `holes=(3,5,6,13)`.
