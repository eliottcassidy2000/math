---
id: HYP-2672
title: LRC14 shell-full tail stratification after the new-speed dyadic block
status: OPEN; exact B36 correction and doubled-odd tail exception
source: codex-2026-06-20-S45
depends_on:
  - HYP-2671
  - HYP-2670
  - HYP-2669
  - HYP-2668
  - HYP-2666
  - HYP-2661
related:
  - THM-545
  - HYP-2665
  - HYP-2664
  - HYP-2662
  - HYP-2643
  - OPEN-Q-108
---

# HYP-2672 - Shell-Full Tail Stratification

## Claim

HYP-2671 localizes the post-shell-gate new-speed `1/3` constant to the dyadic
block row

```text
E'=(0,1,2,4,8,12,16,20), w=24,
Delta_w^+/p1 = 1371/4319.
```

The remaining shell-full constant work does **not** have the naive HYP-2670
far-tail `p1/4` form.  Exact B36 evidence refutes

```text
max(E') > 24  ==>  Delta_w^+ <= p1(E')/4
```

by exactly one row in the scanned quotient:

```text
E'=(0,1,2,4,8,14,26,34), w=38,
Delta_w^+/p1 = 966562/3357319.
```

This row is still below `3/10`, with exact gap

```text
3/10 - 966562/3357319 = 406337/33573190.
```

Thus the corrected post-HYP-2671 split is:

```text
finite high pocket:       max(E') <= 14,
dyadic 1/3 block:         HYP-2671's m=4 row,
intermediate 1/4 ledger:  two rows in 21 <= max(E') <= 24,
doubled-odd tail row:     (14,26,34), w=38,
broad post-dyadic decay:  max(E') > 20 stays below 3/10 in B36.
```

## Evidence

Scripts:

```text
04-computation/lrc14_shellfull_tail_stratification_codex_s45.py
04-computation/lrc14_doubled_odd_tail_probe_codex_s45.py
```

Stored outputs:

```text
05-knowledge/results/lrc14_shellfull_tail_stratification_codex_s45.out
05-knowledge/results/lrc14_doubled_odd_tail_probe_codex_s45.out
```

B36 shell-full quotient:

```text
family: E'={0}+{1,2,4,8}+3 extras from [1,36], w=max(E')+1..max(E')+8
rows=39680
```

Band maxima:

```text
finite <=14:          997/2562   at (0,1,2,4,6,7,8,10), w=12
dyadic/new 15..20:    1371/4319  at (0,1,2,4,8,12,16,20), w=24
intermediate 21..24:  1880/6881  at (0,1,2,4,8,12,13,21), w=24
tail 25..30:          932669/4085893 at (0,1,2,4,5,8,27,28), w=31
tail 31..36:          966562/3357319 at (0,1,2,4,8,14,26,34), w=38
```

Threshold counts:

```text
rows above 3/10: 3
rows above 1/3:  1
post-dyadic rows max(E')>20 above 1/4: 3
far-tail rows max(E')>24 above 1/4:    1
far-tail rows max(E')>24 above 3/10:   0
```

The three post-dyadic rows above `1/4` are:

```text
(0,1,2,4,8,14,26,34), w=38, ratio=966562/3357319
(0,1,2,4,8,12,13,21), w=24, ratio=1880/6881
(0,1,2,4,8,14,20,24), w=26, ratio=1315/4858
```

The new B36 tail exception has low fold mass and a doubled-odd packet address:

```text
extras = (14,26,34) = 2*(7,13,17),
w = 38 = 2*19,
fold_recip = 1/34,
odd_carry = ((1,15),(7,2),(13,2),(17,2)).
```

## Doubled-Odd Probe

The focused doubled-odd packet scout checks

```text
E'={0,1,2,4,8,2a,2b,2c},  odd 3<=a<b<c<=29,  w=2d,  c<d<=c+8.
```

It scans `2912` exact rows.  The B36 tail exception is the unique row above
`1/4`, and no tested row exceeds `3/10`:

```text
above 1/4 = 1
above 3/10 = 0
maximum = 966562/3357319 at (a,b,c,d)=(7,13,17,19).
```

This does not prove the full tail, but it gives the next exception-ledger
address: doubled-odd packets need their own finite treatment before the broad
tail can be scalarized.

## Proof Route

HYP-2672 changes the HYP-2670/HYP-2671 route as follows.

1. Keep HYP-2671 as the new-speed `1/3` dyadic block theorem target.
2. Certify the finite high pocket `max(E')<=14`.
3. Certify the finite intermediate `21..24` rows above `1/4`.
4. Add a doubled-odd tail exception ledger, starting with
   `(2*7,2*13,2*17; 2*19)`.
5. Replace the naive far-tail `p1/4` lemma by a broader post-dyadic
   `3p1/10`-style decay target, or by a packet-dependent constant that treats
   doubled-odd rows separately.

## Tournament Analysis

Vertices:

```text
doubled_odd_tail_exception
tail_3/10_decay
intermediate_21_24_ledger
dyadic_block_exception
finite_B13_pocket
```

Pairwise observable: exact `Delta^+/p1` after the shell-1-full quotient.

Switch/gauge: remove HYP-2671's dyadic block first, then sort by max-speed band
and doubled-odd packet address.

Hamiltonian path:

```text
doubled_odd_tail_exception
> tail_3/10_decay
> intermediate_21_24_ledger
> dyadic_block_exception
> finite_B13_pocket
```

Challenged assumption: HYP-2670's B30 far-tail `p1/4` signal is not stable once
the scan reaches B36.  The stable-looking scalar in this layer is `3/10`, after
finite and doubled-odd exceptions are addressed.

## Status

LRC(14) is not proved.  HYP-2672 is an exact bounded correction to the
post-gate shell-full constant ladder.
