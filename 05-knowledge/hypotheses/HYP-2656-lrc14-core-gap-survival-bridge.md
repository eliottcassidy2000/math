---
id: HYP-2656
title: LRC14 fixed-observer core-gap survival bridge - 6/7 far limit and AP-tail bubble ledger
status: OPEN; post-THM-544 bridge hypothesis with exact rational scout
source: codex-2026-06-19-S36
depends_on:
  - HYP-2651
  - HYP-2653
  - HYP-2655
  - HYP-2654
  - THM-541
  - THM-542
  - THM-543
  - THM-544
  - HYP-2644
  - HYP-2648
  - HYP-2652
related:
  - THM-523
  - HYP-2650
  - HYP-2638
  - HYP-2646
  - OPEN-Q-108
---

# HYP-2656 - LRC14 Core-Gap Survival Bridge

## Claim

The fixed-observer core-gap crux in THM-523 now splits into two interacting
layers:

1. The near-collar AP-tail layer proved through THM-541/542/543/544.
2. A genuinely-far survival layer governed by the same decorrelation engine as
   KPS HYP-2653 and the multiscale warning from KPS HYP-2655, but with the
   fixed-observer sign reversed.

For a positive core `C`, write

```text
G_C = {t in [0,1): ||c t|| > 1/14 for every c in C}.
```

If a speed `w` is genuinely far from a bounded base `B`, then the expected
limit is

```text
meas(G_{B union {w}}) -> (6/7) meas(G_B).
```

The reason is direct: on an already-safe interval for `B`, an independent
`w*t` survives the forbidden central sector with probability `6/7`.  This is
the fixed-observer sibling of HYP-2644's sector-cover plateau and KPS
HYP-2653's decorrelation route.  KPS HYP-2655 warns that a decoupled small
uniform constant is false in multiscale rows, so the far side should now be a
joint plateau/Delta recursion rather than a single one-shot discrepancy bound.

## Evidence

Script:

```text
04-computation/lrc14_core_gap_survival_bridge_codex_s36.py
```

Stored output:

```text
05-knowledge/results/lrc14_core_gap_survival_bridge_codex_s36.out
```

The multi-far ledger scans bounded base cores of sizes `r=1..11` and multiplies
the base minimum by `(6/7)^(12-r)`.  Every tested row remains far above the
collar.  The closest case is the 11-core

```text
B11 = (1,2,3,4,5,7,8,9,11,12,13),
meas(G_B11) = 313/9702,
(6/7)meas(G_B11) = 313/11319,
```

with margin

```text
313/11319 - 7/858 = 5737/294294 = 0.019494111...
```

over the collar.

The exact positive 12-core extension through `B=21` confirms the mainline
near-collar story: the drop-6 collar is still the unique minimum, but the
second value below the old `B<=19` threshold is the now-proved THM-543
exception

```text
C20 = (1,2,3,4,5,7,8,9,11,12,13,20),
meas(G_C20) = 3859/420420 = 7/858 + 1/980.
```

THM-543 proves the whole one-replacement AP-tail layer and identifies
`(a,b,r)=(6,10,20)` as the unique sub-`426/35035` replacement.  The contribution
of this note is the component-owner anatomy and the far-survival bridge around
that theorem.

THM-544 then proves the two-replacement AP-tail layer has no rows below
`426/35035` at all.  That pushes the bounded AP-tail obstruction beyond the
first two replacement layers, while this note keeps the exact bubble anatomy of
the one exceptional row as a template for later endpoint-owner ledgers.

## The 10->20 Bubble Ledger

The THM-543 exceptional row is the drop-6 collar with `10` replaced by `20`:

```text
(1,2,3,4,5,7,8,9,10,11,12,13)
  -> (1,2,3,4,5,7,8,9,11,12,13,20).
```

The component ledger shows that the four old collar components remain
unchanged, and the graft adds exactly two symmetric bubbles:

```text
[29/98, 83/280]   length 1/1960   owners 7 -> 20
[197/280, 69/98]  length 1/1960   owners 20 -> 7
```

Thus

```text
meas(G_C20) = 7/858 + 2*(1/1960) = 7/858 + 1/980 = 3859/420420.
```

This is the fixed-observer anti-coset signal: doubling the missing collar speed
does not create broad tail mass; it creates two addressed endpoint-owner
bubbles.  The scalar `sumset_excess` gets the order wrong here, while the
component ledger sees the theorem's exception immediately.

## Proof Route

The current OPEN-Q-108 route can be read as:

1. Use THM-541 for the AP-window single-hole collar.
2. Use THM-542 for the one-tail mouth retention layer.
3. Use THM-543 for the whole one-replacement AP-tail layer, including the
   unique `10 -> 20` graft.
4. Use THM-544 for the two-replacement AP-tail layer: no row in that layer lies
   below `426/35035`.
5. Prove a fixed-observer far-survival estimate or recursion

```text
|meas(G_{B union {w}}) - (6/7)meas(G_B)| <= C(B)/w
```

but treat KPS HYP-2655 as a guardrail: the useful statement is probably not a
small global constant, but a joint plateau/Delta recursion whose base plateaus
shrink faster than multiscale discrepancies accumulate.
6. Keep the HYP-2648/HYP-2652 address data until after scalar bounds are proved:
   endpoint owners, state words, and component lengths are the proof objects;
   raw far speed and sumset excess are routing labels.

## Tournament Analysis

Vertices are proof quotients, not runners:

```text
proved_one_replacement_tail
proved_two_replacement_tail
exact_core_gap_components
drop6_collar_ledger
owner_bubble_ledger
far_survival_multiplier
joint_plateau_delta_recursion
finite_resonance_discrepancy
sumset_excess
raw_far_speed
```

Pairwise observable: which quotient preserves the lower-bound predicate and
explains the sub-`426/35035` near-collar rows.

Switch/gauge: keep fixed-observer positivity.  Judge tails only after the
proved AP-tail layers, the `6/7` survival quotient, and the KPS HYP-2655
joint plateau/Delta recursion are invoked.

Hamiltonian path:

```text
proved_one_replacement_tail
> proved_two_replacement_tail
> exact_core_gap_components
> drop6_collar_ledger
> owner_bubble_ledger
> far_survival_multiplier
> joint_plateau_delta_recursion
> finite_resonance_discrepancy
> sumset_excess
> raw_far_speed
```

Fingerprint: transitive proof-priority tournament, score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}`, no directed 3-cycles.

## Assumption Challenge

Alternate vertex sets considered: runners, far speed `w`, missing holes, safe
components, component endpoint owners, sumset excess, state words, and proof
obligations.

The challenged assumption is that "far speed" or ordinary sumset excess is the
determinant.  The `B=20` row has exact excess `9`, while the older `B<=19`
second row has excess `1`; excess alone has the order backwards.  The
determining object is the endpoint-owner component ledger: it records that the
new row is the old collar plus two addressed bubbles, exactly matching the
THM-543 exceptional replacement.

The later S15 dyadic/apex-prime residue work is compatible with the same
lesson: dyadic richness and QR/NQR sector structure can be useful tiebreaker
addresses near the `m/7` cover walls, but HYP-2656 treats them as part of the
HYP-2648/HYP-2652 address layer, not as a replacement scalar determinant.
