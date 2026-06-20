---
id: HYP-2673
title: LRC14 constant stack and scale-cluster budget
status: OPEN; synthesis after HYP-2671 and incoming HYP-2653c
source: codex-2026-06-20-S46
depends_on:
  - HYP-2672
  - HYP-2671
  - HYP-2670
  - HYP-2661
  - HYP-2655
  - HYP-2653
related:
  - HYP-2666
  - HYP-2664
  - HYP-2644
  - OPEN-Q-108
---

# HYP-2673 - LRC14 Constant Stack

## Claim

The "one open constant" should not be treated as one scalar.

There are two proof currencies:

```text
local shell-full currency:  Delta_w^+ / p1(E')
global far-peel currency:   w*|Delta_w|
```

HYP-2671 identifies the local shell-full constant stack after the shell-1 gate:

```text
finite pocket:  Delta_w^+/p1 <= 2/5
new-speed:      Delta_w^+/p1 <= 1/3
far tail:       Delta_w^+/p1 <= 1/4   (suggested by B30)
```

Incoming HYP-2653c corrects the global far-peel side:

```text
C(k)=sup w*|Delta_w| is not the old non-resonant ~1.95 constant.
C(k) should be O(number of scale clusters) <= O(k).
```

These are compatible, not competing.  The shell-full constants are boundary-tax
constants normalized by the available single-missed-sector mass `p1(E')`.  The
far-peel constant is a scale-cluster resonance budget normalized by `w`.

## Exact Ledger

Script:

```text
04-computation/lrc14_constant_stack_codex_s46.py
```

Output:

```text
05-knowledge/results/lrc14_constant_stack_codex_s46.out
```

Main exact constants:

```text
shell-damage threshold       = 426/35035
finite packet tax            = 2/5
new-speed packet tax         = 1/3
far-tail packet tax          = 1/4
KPS scale-count scout C(9)   ~= 293/100
```

Exact slacks from the output:

```text
drop-6 mouth floor                         = 7/858
426/35035 - 7/858                          = 841/210210

delete 1 margin above 426/35035            = 73673/2382380
delete 2 margin above 426/35035            = 3153/595595
delete 4 margin above 426/35035            = 1927/805805
delete 8 margin above 426/35035            = 1481/240240

finite B13 leader gap below 2/5            = 139/12810
new-speed dyadic-block gap below 1/3       = 206/12957
far-tail B30 leader gap below 1/4          = 355217/16343572
```

The binding local gate remains the missing-`4` tower deletion, but it is still
strictly above the threshold.  The binding shell-full new-speed row remains the
dyadic block

```text
E'=(0,1,2,4,8,12,16,20), w=24,
Delta_w^+/p1 = 1371/4319.
```

## S46 Bridge Addendum

Incoming KPS work after the first S46 push stated that the codex shell-full
`p1`-tax object and the KPS far-plateau deviation are the same exact quantity.
This was verified in:

```text
04-computation/lrc14_codex_kps_bridge_verify_codex_s46.py
05-knowledge/results/lrc14_codex_kps_bridge_verify_codex_s46.out
```

For the HYP-2671 dyadic-block extremizer:

```text
raw_wdelta(E',w)/w = p0(E' union {w}) - Phi(E') = 457/3920
Delta/p1 = 1371/4319
```

For the incoming non-shell-full warning row:

```text
E'=(0,2,3,5,6,15), w=18,
raw_wdelta(E',w)/w = p0(E' union {w}) - Phi(E') = 11/315,
Delta/p1 = 22/63 > 1/3.
```

Thus the relative `Delta <= p1/3` route is genuinely shell-full; non-shell-full
rows must be routed through the shell-damage gate or the absolute `C(k)/w`
far-peel route.

## Far-Span Reading

For the tight k=9 far-peel row, the sector margin is

```text
cap_9 - Q(8) = 129643/980980 ~= 0.132157.
```

Thus different `C` targets ask for different finite base spans:

```text
C=39/20       -> span >= 15
C=293/100     -> span >= 23
C=2804017/717360 -> span >= 30
C=9           -> span >= 69
```

This is the useful correction.  The old `1.95` target was convenient because it
dovetailed with the span-16 bank, but the theorem should instead prove a
per-scale-cluster bound and enlarge the finite base span as needed.  Even the
known multiscale examples point to bounded-span work, not a proof failure.

## Proof Route

The current route should be assembled in this order:

1. Prove the HYP-2661/HYP-2666 shell-damage gate:
   damaging `{1,2,4,8}` pays at least `426/35035`.
2. On the shell-1-full quotient, prove the finite `max(E')<=14` packet ledger,
   with the B13 leader below `2p1/5`.
3. Prove HYP-2671:
   `max(E')>14` shell-full new-speed rows satisfy `Delta_w^+ <= p1/3`, with
   the `m=4` dyadic block isolated as the sharp scout row.
4. Prove the far-tail packet lemma suggested by HYP-2670:
   `max(E')>24` rows satisfy a `p1/4`-style decay.
5. For genuinely wide/multiscale far peels, prove
   `w*|Delta_w| <= O(number of scale clusters)` jointly with plateau shrinkage,
   then finish by a finite base span around `23-30` for k=9 unless the constant
   is sharpened.

## Tournament Analysis

Vertices:

```text
shell_damage_gate
finite_packet_tax
new_speed_packet_tax
far_tail_packet_tax
scale_cluster_Ck
raw_runner_vertices
```

Pairwise observable: exact proof currency attached to each gate.

Switch/gauge: normalize endpoint discrepancy either by `p1(E')` in the local
shell-full quotient or by `w` in the global far-peel quotient.

Hamiltonian path:

```text
shell_damage_gate
> finite_packet_tax
> new_speed_packet_tax
> far_tail_packet_tax
> scale_cluster_Ck
> raw_runner_vertices
```

Challenged assumption: there is one scalar constant to prove.  The evidence
supports a stack of currencies: local boundary mass and global scale-cluster
count divided by `w`.

## Honest Status

LRC(14) is not proved.  HYP-2673 refines the target after the incoming
HYP-2653c update: prove a local packet-tax stack plus a global per-scale-cluster
resonance budget, rather than chasing the obsolete uniform `C~=1.95`.
