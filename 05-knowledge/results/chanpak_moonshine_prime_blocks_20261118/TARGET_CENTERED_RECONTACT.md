# Target-centered repair slices beyond order six

## Experiment and exact physical positive windows

Using mt19937 seed351324, sample300 labelled deletion tournaments at each total order7,8,9,10. For each host compute its exact quadratic apex response from empty, singleton and pair evaluations by integer subset-DP. Enumerate every contact setting, every distinguished single repair spoke with signed increment+2, both fixed orientations of that spoke, and every move changing at most two untouched spokes whose baseline jump is≤4. Each connected component certifies its entire odd interval by the previously proved crossing lemma.

The largest intervals found had four odd values. This is a sampled-host observation, NOT a maximum-size theorem. In particular these exact windows occur well above the previous9..15 control:

|order|two physical windows|common edit mask|
|---|---|---|
|7|129/131 and133/135|17|
|8|593/595 and597/599|33|
|9|1487/1489 and1483/1485|17|
|10|1725/1727 and1721/1723|17|

Each common edit changes exactly two untouched apex spokes. `high_repair_window_replay.py` gives complete outgoing masks and replays every literal tournament without using the quadratic evaluator. Thus this is actual local exact-target transport, not predicted interpolation. The n10 interval includes1727; the data say nothing special about the taxicab number1729.

## Direct dyadic-source comparison

Take the actual five-vertex source code70, mark vertex0, and inflate it to a transitive block of k vertices. Its count is6·2^k−1 by the independently derived skeleton compiler. For each k=1,...,8, exhaust EVERY apex of this fixed tournament, all recontacts at that apex, all distinguished single-spoke repairs, and all connected small-step components above.

At k=1 the component9..15 contains both center±1 for center12. At every k=2,...,8 (orders6..12), the largest component containing both6·2^k−1 and6·2^k+1 consists of exactly those two values. The centers are24,48,96,192,384,768,1536. No growth in repair radius occurs in this exact move universe. The largest other interval in these recontact orbits is1,3,5.

This is a decisive finite boundary on the most direct proposed composition: inflating the center does not automatically inflate the local repair component. It does NOT rule out the other source of the dyadic pair, a changing deletion tournament, a multiple-spoke distinguished repair, simultaneous multiple apices, or larger permitted exchange jumps. In particular it is not a bound on all repair windows at these orders, as the preceding positive examples demonstrate.

## How the random historical pass changed the question

A second deterministic sample (seed192413 from the sorted result Markdown list) selected:

* `seam_mills_20260925.md`: its decoder must retain both quotient and shell remainder. Having a point somewhere in a shell is not placing one at its prescribed offset.
* `collatz_golden_routing_synthesis_20261007.md`: distinguish a genuinely new guarded route from merely reaching an ungrounded dependency; keep the source and terminal obligations.
* `procgen_continuous_20260925_natural_boundary_mahler_harmonic.md`: growth/transience and membership are different invariants. No Collatz analytic theorem is imported here.

These suggested measuring the component CONTAINING the dyadic target, rather than reporting the largest unrelated component. That stronger test detects the failure already at k=2.

## Next concrete direction

The positive high windows show that changing the deletion source can produce exact local transport where the fixed inflated source has none. Thus store a source's repair-slice equation and its component baseline interval as one marked certificate. A selector must choose among such certificates by absolute target, not just retain H or the exponent k. The next genuinely different test is to change the deletion core while preserving a boundary profile, or use a multiple-spoke repair D that leaves a less rigid affine gap equation. Neither is proved available universally.

## Reproduction and scope

Run `python3 high_repair_window_replay.py` for the literal witnesses, and `python3 dyadic_recontact_replay.py` for the complete recontact scan on the listed dyadic sources. Seeded random search was used to discover, not certify, the high windows. No maximality claim is made.
