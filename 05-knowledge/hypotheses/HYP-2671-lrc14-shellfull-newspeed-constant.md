---
id: HYP-2671
title: LRC14 shell-full new-speed 1/3 constant and dyadic block resonance
status: OPEN; exact B=30 constant atlas
source: codex-2026-06-20-S45
depends_on:
  - HYP-2670
  - HYP-2669
  - HYP-2668
  - HYP-2661
related:
  - HYP-2667
  - HYP-2666
  - HYP-2665
  - HYP-2643
  - OPEN-Q-108
---

# HYP-2671 - Shell-Full New-Speed Constant

## Claim

The post-shell-gate open constant is the `1/3` new-speed barrier:

```text
shell-1-full and max(E') > 14  ==>  Delta_w^+ <= p1(E')/3.
```

In the exact B30 quotient, the maximum new-speed row is not diffuse.  It is a
single dyadic block resonance:

```text
E'=(0,1,2,4,8,12,16,20), w=24,
Delta_w^+/p1 = 1371/4319,
1/3 - Delta_w^+/p1 = 206/12957.
```

This row is exactly the `m=4` member of

```text
E_m={0,1,2,4,8,3m,4m,5m}, w=6m.
```

Thus the proof should isolate and certify this dyadic block first, then prove
extra packet cancellation for every other shell-full new-speed row.

## Evidence

Script:

```text
04-computation/lrc14_newspeed_constant_codex_s45.py
```

Stored output:

```text
05-knowledge/results/lrc14_newspeed_constant_codex_s45.out
```

B30 recap:

```text
rows=20800
finite max = 997/2562 at (0,1,2,4,6,7,8,10), w=12
new-speed max = 1371/4319 at (0,1,2,4,8,12,16,20), w=24
tail max = 932669/4085893 at (0,1,2,4,5,8,27,28), w=31
```

The best rows by `max(E')` show that the only new-speed layer close to `1/3`
in the B30 scan is `max(E')=20`:

```text
max(E')=20: 1371/4319
max(E')=21: 1880/6881
max(E')=24: 1315/4858
max(E')=28: 932669/4085893
```

The dyadic block family has a unique spike at `m=4`:

```text
m=3: 205/5817
m=4: 1371/4319
m=5: 131/3346
m=6: 122/2331
...
m=24: 523/6580
```

So the `1/3` constant is the rational gap above `1371/4319`, not a broad
moving frontier.

## S46 Long-Family Addendum

Script:

```text
04-computation/lrc14_dyadic_block_long_family_codex_s46.py
```

Output:

```text
05-knowledge/results/lrc14_dyadic_block_long_family_codex_s46.out
```

This extends the same exact dyadic-block family to `m=3..120`:

```text
global family top: m=4, ratio=1371/4319
tail m>24 top:    m=60, ratio=3048/25067
best m=61..120:   m=78, ratio=3632/31857
```

Thus the `m=4` spike remains isolated in the targeted family through `m=120`.
After `m=24`, the best tested row is only about `0.1216`, far below `1/3` and
far below the B30 non-family new-speed competitors.  This does not prove the
new-speed lemma, but it strengthens the interpretation: `m=4` is a finite
resonance to isolate, not the first term of a high-ratio tail.

## Fold Guardrail

HYP-2670 suggested fold-target dilution as the reason new speeds fall.  HYP-2671
sharpens that to a guardrail:

```text
fold reciprocal mass helps locate the dyadic block, but does not certify the
constant monotonically.
```

Examples from the B30 scan:

```text
new-speed maximum:
  E'=(0,1,2,4,8,12,16,20), fold_recip=59/240, ratio=1371/4319.

low-fold but still moderate:
  E'=(0,1,2,4,8,14,20,24), fold_recip=1/24, ratio=1315/4858.

high-fold but harmless:
  fold bins above 1/2 never approach 1/3 in the B30 scan.
```

Therefore the certificate likely needs both:

```text
dyadic block address + phase-packet cancellation,
not fold mass alone.
```

## Proof Route

The current shell-full proof can be sharpened to:

1. Finite B13 pocket: handle the single row above `1/3` but below `2/5`.
2. Dyadic block resonance: prove the `m=4` new-speed row is the unique
   shell-full new-speed maximizer or is safely dominant in a finite dyadic
   block class.
3. New-speed remainder: prove all non-block shell-full new-speed rows satisfy
   `Delta_w^+ <= p1/3`.
4. Tail split: HYP-2672 shows the stronger raw `p1/4` target is false at B36,
   so replace it with finite intermediate rows, doubled-odd exception ledgers,
   and broad post-dyadic `<3p1/10` decay.

Together with the incoming KPS S17/THM-545 shell-damaged gate, this turns the
two-gate route into two explicit constants:

```text
local shell damage: threshold 426/35035;
post-gate far tax: 1/3 outside finite dyadic pockets, 2/5 inside the B13 pocket.
```

## Tournament Analysis

Vertices:

```text
dyadic_block_resonance
new_speed_1/3_gap
fold_recip_bins
max_speed_layer
raw_runner_vertices
```

Pairwise observable: exact `Delta_w^+/p1` on the shell-full quotient with
`max(E')>14`.

Switch/gauge: first gate by `max(E')`, then by dyadic block family and
fold-target bins.

Hamiltonian path:

```text
dyadic_block_resonance
> new_speed_1/3_gap
> fold_recip_bins
> max_speed_layer
> raw_runner_vertices
```

Challenged assumption: fold mass alone determines the constant.  It does not.

## Honest Status

LRC(14) is not proved.  HYP-2671 identifies the exact bounded extremizer for
the shell-full new-speed `1/3` constant and gives a sharper proof target.
