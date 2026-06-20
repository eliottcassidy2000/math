---
id: HYP-2668
title: LRC14 B14 shell-gated p1 tax - global 2p1/5 fails, shell-full survives
status: OPEN; exact B=14 bounded frontier
source: codex-2026-06-19-S42
depends_on:
  - HYP-2667
  - HYP-2666
  - HYP-2665
  - HYP-2664
related:
  - HYP-2669
  - HYP-2661
  - HYP-2662
  - HYP-2663
  - HYP-2655
  - OPEN-Q-108
---

# HYP-2668 - B14 Shell-Gated p1 Tax

## Claim

The global raw scalar target

```text
Delta_w^+ <= 2*p1(E')/5
```

does not survive the next bounded layer.  The full `B=14` exact scan finds one
row above `2/5`.  However, that row is shell-1 damaged, and every shell-1-full
row in the scan remains below `2/5`.  This supports HYP-2666's ordered proof
shape:

```text
first apply the shell-1 gate;
then prove the p1-tax theorem on the shell-1-full quotient.
```

The live target is therefore:

```text
if E' is shell-1-full, Delta_w^+ <= 2*p1(E')/5;
if E' damages shell 1, route to HYP-2661 tower/mouth rigidity first.
```

## Evidence

Script:

```text
04-computation/lrc14_p1_tax_b14_shell_gate_codex_s42.py
```

Stored output:

```text
05-knowledge/results/lrc14_p1_tax_b14_shell_gate_codex_s42.out
```

Exact bank:

```text
E'={0}+7-subsets of [1,14],
w=max(E')+1,...,max(E')+8.
```

It contains `27448` primitive rows.  Exact threshold counts:

```text
Delta_w^+/p1 > 1/3: 32 rows
Delta_w^+/p1 > 3/8: 5 rows
Delta_w^+/p1 > 2/5: 1 row
Delta_w^+/p1 > 5/12: 0 rows
```

The unique `2/5` failure:

```text
E'=(0,1,4,6,8,10,12,14), w=16
shell1_missing=(2,)
Delta_w^+/p1 = 7071/17584 ~= 0.402127
2/5 tax gap = -374/25725
5/12 tax gap = 1534/15435
```

But the shell-1-full maximum is still the S41 B=13 row:

```text
shell1_missing=()
E'=(0,1,2,4,6,7,8,10), w=12
Delta_w^+/p1 = 997/2562 ~= 0.389149 < 2/5.
```

HYP-2669 later extends this shell-1-full stability through `B=24`: the same
row remains the maximum, with exact `2/5` tax gap `139/2450`.

So, in the B=14 scan:

```text
global max        = 7071/17584
shell1-full max   = 997/2562
shell1-damaged max= 7071/17584
```

The full B=14 bank remains below `5/12`, with minimum tax gap `1534/15435`.
But `5/12` should not become the proof target: it is evidence about raw Delta,
while the cap-splice wants a smaller constant after the shell gate.

## Interpretation

HYP-2667 correctly identified the dyadic-even packet frontier, but its raw
global `2/5` option was too optimistic.  The B=14 counterexample misses the
shell-1 tower bit `2`, so it belongs to the HYP-2661/HYP-2666 gate, not to the
shell-1-full p1-tax quotient.

This reconciles the two incoming routes:

1. HYP-2666's cap-bank scout says shell-1 stratification is the right first
   quotient.
2. HYP-2667 says raw constants must be tested against actual Delta, not only
   cap slack.
3. HYP-2668 says the current proof obligation is shell-gated.
4. HYP-2669 says the shell-full half remains stable through `B=24`, so the
   live local target is the finite dyadic-even packet around the B13 leader,
   with shell-damaged rows discharged separately.

## Tournament Analysis

Vertices:

```text
shell1_gate
shell1_full_p1_tax
global_2p1/5
global_5p1/12
raw_ratio
```

Pairwise observable: whether the next bounded layer preserves the proposed
p1-tax target.

Switch/gauge: stratify exact ratios by missing shell-1 packet before choosing
the scalar.

Hamiltonian path:

```text
shell1_gate
> shell1_full_p1_tax
> global_2p1/5
> global_5p1/12
> raw_ratio
```

Challenged assumption: the B=13 scalar `2p1/5` target globalizes without the
shell gate.  It does not.

## Honest Status

LRC(14) is not proved.  HYP-2668 sharpens the current proof obligation to a
two-gate theorem: prove HYP-2661-style shell damage pays first, then prove
`Delta_w^+ <= 2*p1(E')/5` on the shell-1-full quotient, sharpened by
HYP-2669 to a finite dyadic-even packet plus new-speed decay target.
