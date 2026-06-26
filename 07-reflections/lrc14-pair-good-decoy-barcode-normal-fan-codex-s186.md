# LRC14 Pair-Good Decoy Barcode And Normal-Fan Refinement

codex-2026-06-26-S186

## What Changed

The useful move is to stop treating pair-good decoys as a count.  HYP-3021
gives the modular tooth generator; this pass adds where the generated false
switch sits relative to HYP-3015 barcode bars and HYP-3018 active normal-fan
supports.

A decoy is a pair-crossing time where some pair is threshold-good, but another
runner makes the full minimum fall below `1/14`.

The exact generator is:

```text
pair lane + shell + denominator + pair gap + blocker deck + barcode/normal-fan relation
```

That turns a scary multiplicity into a grammar of false switches.

## Readout

Named rows have `8889` pair-good decoys, but they split by explicit generators.
The common blockers are small core runners such as `(7,)`, `(11,)`, `(1,)`,
`(9,)`, and `(13,)`.  The bounded low-frontier one-swap atlas
`drops={10,12,13}, add<=36` has `48037` decoys and `9809` generator classes,
again dominated by exact single-blocker decks.

Every named decoy is outside a strict barcode bar.  Many still overlap real
normal-fan supports through the good pair or blocker deck, which means they can
be routed to the same support certificate.  The disjoint cases are the pure
false switches and should be discharged by lane and blocker alone.

## Assumption Challenge

The vertices are not runners or raw pair gaps.  They are generator classes and
proof carriers.  The quotient preserves the LRC predicate by retaining why a
pair-good time is not a lonely time.

## Next Pull

Add pair-good generator fields to the HYP-2963 packet bank next to HYP-3021's
tooth fields, then prove a blocker-deck lemma for the common singleton
blockers.  After that, decoy counts should be bookkeeping, not a proof risk.
