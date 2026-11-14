# A count-nine cut family, Hamilton-cycle enumeration, and the metagraph

Continues the [minimum right-fiber theorem](
../chanpak_h3_right_fiber_20261114/PROOF.md) and
[selective-cut inverse](../chanpak_selective_target_20261114/README.md).
The session gives new absolute selectors and exact structural
enumerators. **Chan–Pak remains open.**

## Fixing the family

* [Two H=9 routes](H9_FIBERS.md): two triangle components still
  cannot supply three consecutive odds in a common right fiber.
  Their nonmultiple-of-three offsets require a local four-bit
  interaction with separated residue types. A single strong
  component instead gives the four-odd right window 89–95,
  which survives every common transitive continuation. A small
  atlas inverts any prescribed strong H=9 fiber without searching
  the two transitive bit strings.
* [Two-height windows and prime/7p guards](h9_bilateral/PROOF.md):
  genuine bilateral pairs and triples, including a guarded
  unbounded two-height construction and literal prime/7p pairs.
  The guards do not prove infinitely many primes on the rails.
* [Square-root target selector at fixed cycle count](
  h9_bilateral/DIAGONAL_EDGE.md): one actual apex arc reversal
  changes the Hamilton path count by two while keeping exactly
  nine Hamilton cycles. A second set of three physical cuts
  supplies a triple. From N the possible cut address is forced
  by floor(sqrt(N/84)). The family gives order sqrt(X) outputs
  through X, not a cofinite spectrum. The arc-flip operation is
  not identified with a path-complement metagraph line.
* [Exact size enumerators](CYCLE_NINE_ENUMERATOR.md): at order
  2ℓ+6 the diagonal families give 5·2^(ℓ−1) rigid classes
  with distinct H and cycle count 9, and a closed rational
  generating function for their full-path presentations.
  Independent tapes give 5(n−9)2^(n−10) rigid classes
  at cycle count 9 for every n≥10, even after converse
  merging, without requiring their path counts to differ.

* [Guarded path-to-cycle chart](h9_bilateral/GUARDED_ENUMERATOR.md):
  any realized core path count becomes the cycle count of a
  whole canonically recognizable family at all sufficiently large
  orders. Its exact ordinary, labelled, converse-merged, and
  cycle-marked series reduce to marked-core inventories. This is
  a structural lift, not a converse from all cycle counts to paths.

* [General guarded lift](h9_bilateral/GUARDED_ENUMERATOR.md):
  an exact ordinary, labelled, and converse-merged class
  series that turns any trimmed core's path count into a
  cycle count. A unique recoverable apex makes its cycle
  action free, so derivatives of the weighted series count
  unrooted cycle marks and rooted closed cuts.

## What unfolding a metagraph node exposes

* [Cycle necklaces and fans](CYCLE_NECKLACES.md): exact counts
  of closed cuts at each rotational stabilizer and of their
  blue (mirror-fixed) cuts. A marked Hamilton cycle generates
  a fan of open partners, with exact mean cyclic-triangle
  change −1 before weighting the converse merge. Its nonstrong
  partners can be read as strong bridges in the cycle-reversed
  source. This is not a monotone result for Hamilton path counts.
* [Mixed path/cycle moments](MIXED_MOMENTS.md) and
  [partial cycle markers](PARTIAL_CYCLE.md): closed counting
  formulas for joint moments and marked isomorphism classes.
  Closed endpoints have more Hamilton cycles on average in
  the weighted line graph; individual lines can go the other way.
* [Same spectra, different lines](SAME_SPECTRA_LINES.md):
  a pair of distinct rigid order-seven merged nodes agrees on
  all path and cycle counts at *every* length, yet has three
  versus six closed lines into nonstrong mates. Therefore even
  complete marginal spectra cannot serve as a stand-alone
  state for the metagraph.
* [Colored scaffold Burnside](PACKING_BURNSIDE.md): an exact
  bridge from A000568's unmarked partition sum to path, cycle,
  and packing decorations. The odd-cycle packing identity
  equates path and packing cardinalities but not their
  isomorphism orbits for any tournament with symmetry.

## Reproduce

    python replay_h9.py
    python h9_bilateral/replay.py
    python replay_enumerator.py
    python replay_cycle_necklaces.py
    python replay_mixed.py
    python replay_partial.py
    python metagraph_collision.py
    python replay_packing.py

The optional C++ subset DP in h9_bilateral independently
recounts literal path and cycle totals on the guarded witnesses.
The remaining Chan–Pak obligation is still supplying a physical
profile for an arbitrary prescribed odd target, not computing
more witnesses inside a fixed chart.
