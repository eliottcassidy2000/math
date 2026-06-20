---
id: HYP-2669
title: LRC14 shell-full p1-tax quotient stability through B24
status: OPEN; exact shell-1-full bounded quotient evidence
source: codex-2026-06-20-S43
depends_on:
  - HYP-2668
  - HYP-2667
  - HYP-2666
  - HYP-2661
related:
  - HYP-2665
  - HYP-2664
  - HYP-2662
  - OPEN-Q-108
---

# HYP-2669 - Shell-Full p1-Tax Quotient Stability

## Claim

After the HYP-2668 shell gate is imposed, the raw p1-tax target

```text
Delta_w^+ <= 2*p1(E')/5
```

survives the corrected exact shell-1-full quotient through `B=24`.  The maximum
does not move from the S41/B13 dyadic-even leader:

```text
E'=(0,1,2,4,6,7,8,10), w=12,
Delta_w^+/p1 = 997/2562 < 2/5,
2/5 tax gap = 139/2450.
```

Thus the next proof obligation is not a growing scalar frontier in the
shell-full quotient.  It is a finite phase-packet lemma around the B13 leader,
plus the separate HYP-2661/HYP-2666 discharge of shell-damaged rows.

## Evidence

Script:

```text
04-computation/lrc14_shellfull_p1_tax_quotient_codex_s43.py
```

Stored output:

```text
05-knowledge/results/lrc14_shellfull_p1_tax_quotient_codex_s43.out
```

Exact quotient:

```text
E'={0}+{1,2,4,8}+3 extras from [1,24],
w=max(E')+1,...,max(E')+8.
```

The scan has `9120` primitive rows.  Cumulative checkpoints:

```text
B=14: rows= 960, >1/3=1, >3/8=1, >2/5=0, max=997/2562
B=16: rows=1760, >1/3=1, >3/8=1, >2/5=0, max=997/2562
B=18: rows=2912, >1/3=1, >3/8=1, >2/5=0, max=997/2562
B=20: rows=4480, >1/3=1, >3/8=1, >2/5=0, max=997/2562
B=22: rows=6528, >1/3=1, >3/8=1, >2/5=0, max=997/2562
B=24: rows=9120, >1/3=1, >3/8=1, >2/5=0, max=997/2562
```

Rows with genuinely new speeds beyond `B=14` are well below the target:

```text
new-speed max = 1371/4319 ~= 0.317435
E'=(0,1,2,4,8,12,16,20), w=24.
```

Odd-carry profiles also separate the leader from the new-speed rows.  The B13
leader keeps the full shell-1 carry and uses medium odd shells:

```text
((1, 15), (3, 2), (5, 2), (7, 1)).
```

The new-speed leader extends the dyadic-1 tower instead:

```text
((1, 31), (3, 4), (5, 4)).
```

This supports a proof split where extended-tower and far-new-speed rows decay,
while the only close shell-full packet is the finite dyadic-even row already
visible at B13.

## Interpretation

HYP-2668 refuted a global pre-gate `2p1/5` bound.  HYP-2669 says the corrected
post-gate version is still alive after a larger exact quotient scan.  This
does not prove the theorem, but it changes the shape of the remaining work:

1. Shell-damaged rows should be routed through HYP-2661/HYP-2666 first.
2. Shell-full rows with new large speeds should be controlled by a far/new-speed
   decay or fixed-observer survival argument.
3. The only close shell-full packet is the finite dyadic-even B13 leader, so
   the local packet lemma should focus there.

## Tournament Analysis

Vertices:

```text
b13_shellfull_leader
dyadic_even_packet
new_speed_decay
extended_tower
raw_ratio
```

Pairwise observable: whether forcing shell 1 lets the `2p1/5` p1-tax target
survive larger bounded quotients.

Switch/gauge: quotient first by the full dyadic-1 tower, then compare extra
odd-carry profiles.

Hamiltonian path:

```text
b13_shellfull_leader
> dyadic_even_packet
> new_speed_decay
> extended_tower
> raw_ratio
```

Challenged assumption: the S42 shell-full survival might be only a B14 artifact.
Through B24 it is not.

## Honest Status

LRC(14) is not proved.  HYP-2669 is bounded evidence for the shell-full half of
the HYP-2666/HYP-2668 two-gate proof obligation.
