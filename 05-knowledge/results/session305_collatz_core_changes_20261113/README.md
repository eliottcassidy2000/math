# Session 305: changing cores and compiling Collatz returns

**Result:** a finite-proof reduction to four reverse ladder depths and two ready seeds; physical transport of a guarded triangle factor at unchanged count; and a matching-core move scaling both a binary and a ternary clock exactly. **Chan–Pak and both one-sign Collatz problems remain open.**

## 1. From arithmetic to a concrete physical obligation

The [elementary controller](ARITHMETIC.md) descends every allowed odd source other than 1 and 5 in at most ten signed odd steps, using valuations at most five. Unlike Session 304's endpoint-only route, every intermediate in the **reverse dyadic ladder** also avoids 7 and 21. At cap five the worst block-length requirement is exactly ten, and the best possible worst-case expanded peak factor is 195/19; both are forced by source 19. The seed 5 cannot be eliminated from this reverse arithmetic grammar: the minus step hits 7, and the plus step's reverse ladder also hits 7.

An [exact dyadic cylinder partition](CAP4.md) reduces the maximum valuation to **four**, with length at most 26. Its 1,068 disjoint cylinders cover all odd residues at precision 43, and the checked affine thresholds certify the entire infinite range n≥65. No finite census of counts implies this result. Short individual words cover the low boundary. Valuation cap two cannot give uniformly bounded descent: arbitrarily large positive integers can shadow the 2-adic alternating cycle at ±1/5, where its only allowed exits are valuation-1 expansions and the alternative sign has valuation 3. Cap three is not decided here.

The reverse direction is a target-decreasing induction at **block boundaries**. For each arithmetic edge with 3n+s=2^a m, reverse from m by one 2x−s, then a−1 moves 2x+s, and then divide the resulting 3n by 3. Having a smaller count in the induction hypothesis does not supply the bigger intermediate counts or their ports. A sufficient physical compiler must realize the selected entire reverse blocks from ready representatives, with ready seeds at 1 and 5. It does not need a universal lift of every one-step edge, or of either single-sign Collatz system.

## 2. Move an actual factor of 3 at unchanged H

[Selective shared-apex fusion](APEX_FUSION.md) gives the exact total count from seven local path and endpoint statistics per ordinal block. With a uniform apex port a block's normalized response depends only on its two-cover clock q=(2H+F₂)/H, even with selective ports elsewhere. TT₃ and C₃→TT₁ have the same clock q=8 but counts 1 and 3. Swapping their roles in two slots transports an actual peelable triangle at fixed total H and fixed order. This gives a neutral preparation when the count already carries a usable factor 3 in the wrong slot. It neither creates an additional factor nor proves repeatable reentry. A third visit to a block can detect structures that the two-visit clock cannot distinguish; the packet exhibits this failure explicitly.

## 3. Scale coupled clocks rather than fitting totals

[Marked backward matching cores](MATCHING.md), drawn from the older zero-heavy structure work, admit a signed monomial compiler. Two endpoint matchings yield the same ternary mode but different binary modes. Switching four endpoint arcs can isolate exactly a difference of 2^b at the same order. A separate one-step width exchange keeps the ternary coefficient fixed and multiplies the binary coefficient by 3/2 while extending the central gap. The exact result is 3H−2 on an infinite family of actual strong tournaments. The particular exchange cannot repeat because its prerequisite outer difference changes from one to three. This is not another fixed one-apex rearrangement; it points to a larger palette of matching cores as the next place to seek a closed return.

## 4. What sparse failure does and doesn't buy

[Finite-horizon sign repair](RARE_REPAIR.md) follows ordinary Collatz until a selected block descends or would fail the chosen K-step horizon. A minus step and a halving repair each of those exceptional starting sources. The exact source set is a dyadic prefix-survival set plus a finite boundary, and its density decays exponentially in K. This does not bound the number of sign changes along a given ordinary trajectory or solve ordinary Collatz. Its lesson for Chan–Pak is concrete: a vanishing density of failed count projections does not dispose of surviving interface states. A neutral transport or an alternative physical macro must cover those states individually.

## Replay

Run the following next to this file (all run in normal and optimized Python without dependencies):

- `python verify_arithmetic.py`: exact bounded arithmetic paths and obstructions.
- `python verify_cap4.py`: exact disjoint cover, valuations, thresholds and small boundary.
- `python check_fusion.py`: independent subset DP controls for paths, cycles and path covers.
- `python check_matching.py`: independent DP for two cores, matching-mode and exchange checks.
- `python check_rare_repair.py`: finite parity-cylinder and repair controls.
