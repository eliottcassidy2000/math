# Reservoir descent and the missing half of the wiring

2026-11-13. **Chan–Pak remains open.** This packet proves exact reservoir
changes at infinitely many absolute targets, sharpens what must be
preserved to chain such changes, and recovers a genuinely missing sign
coordinate in the merged metagraph.

## Native results

- [Dyadic switch](DYADIC_SWITCH.md): two different strong reservoirs
  realize the same absolute count 9·2^(k−2)+9 at every k>=4. Three
  different reservoir counts 5, 3, 1 realize the consecutive odd
  totals 21·2^k+15, +17, +19. But different fixed terminal counts
  synchronize on only O(r²) of the first 2^r common high prefixes, even
  allowing a bounded count error. Bounded retiming gives O(r³).
- [Changing reservoir](CHANGING_RESERVOIR.md): a pure exposed-end
  replacement survives any common binary prefix and produces every
  even gap 2r with unbounded reservoir counts. A separate cancellation
  example also gives an exact adjacent pair. The whole literal
  chain-cut grammar gives at most 1/28 of the odd targets
  asymptotically; the gap parameter is not an absolute selector.
- [Bidirectional reservoir](BIDIRECTIONAL_RESERVOIR.md): one cut word
  controls the missing starts in binary from the left and the
  split-end count in base 3/2 from the right. Independent leading
  choices yield an approximate rectangle. At exact precision the
  endpoint is forced. The split-end formula was independently
  recovered by the pure-end and bidirectional calculations.
- [Readiness check](READINESS_CHECK.md): the equal-count dyadic switch
  changes Q_2 by exactly 12−4·3^(k−2). Also fixing the low L bits of
  the missing-start count narrows the normalized end response to a
  single interval of width at most 2(2/3)^L+2(2/3)^m. This is the
  concrete failure of treating a binary guard as a separate knob.

## What merges and what can be recovered

- [Voltage lifts](VOLTAGE_LIFTS.md): retaining ordinary dart
  multiplicities restores the sign-character block. Existence of an
  involutive lift is a cycle-parity test with extra involution
  constraints. Order five still has x^8(1+x²)^2; order six changes
  from x^12(1+x²)^22 for the merged diagram to x^12(1+x^44) for the
  ordinary wiring. All 22 chiral flips must move together. A simple
  general obstruction: the half-turn of the base C4 in its connected
  C8 double cover lifts by rotations by 2 or 6; both square to the
  deck flip. Ordinary liftability need not mean involutive liftability.
- [Character channels](CHARACTER_CHANNELS.md): the actual adjacency
  family has a Fourier Gram decomposition, with exact outer-ring
  tensor growth before quotienting. Normalized class operators
  already fail to commute at order four. Cyclic substitution, rather
  than ring growth, gives a genuine wreath-product law. The marked
  stabilization parallels the [Atiyah–Bott–Shapiro Clifford model](https://www.maths.ed.ac.uk/~v1ranick/papers/abs.pdf);
  neither this character theory nor the number of fixed classes
  proves any mod-eight recurrence for tournament classes.
- [Tight suffix factors](TIGHT_SUFFIX_FACTORS.md): the graceful
  completion staircase really does split into independently
  recombinable marked factors across fixed tight flats. Two-color
  connectivity by itself does not give that product. The resulting
  fixed traces and inverse signed-cube test distinguish a genuine
  factorization from the same polynomial on a path.
- [Descent and glue](DESCENT_AND_GLUE.md): the two-step analogy is to
  test equality on a complete interface and then glue actual
  witnesses. Projecting to the target first admits false joins.
  The older inverse-stress and marked Collatz precedents explain
  why compatibility and an absolute anchor must both survive.
  An unbounded triple component anchored at 23 suffices; an
  unbounded collection of isolated triples does not.

## Next proof obligation

Construct a target-dependent completion of the marked pullback: when
changing the reservoir, choose enough free end bits to land on the
prescribed count **and** certify the new interface can enter the
next move. The present exact end switches certify those as separate
tests, and the coupled-residue and density bounds rule out treating
the current single buffer as two independent selectors.

The replay scripts are standalone Python 3 (use the packet directory):
check_wiring.py, replay_characters.py, replay_switch.py,
replay_readiness.py, replay_bidirectional.py, replay_affine.py.
The finite checks do not stand in for the all-size proofs.
