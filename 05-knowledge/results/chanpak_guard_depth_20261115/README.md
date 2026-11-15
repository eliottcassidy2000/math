# Chan–Pak: exact digits, changing roots, and the absolute scale

**Status:** new all-order results for specified construction families.
Chan–Pak remains open. This session aimed to create absolute integer
choices that survive continuation, not to raise the finite coverage limit.

## Constructive results

- **[Repeatable physical digits](REPEATABLE_CLOCK.md).** Two patches
  with the same continuation state have different responses to the
  accumulated address. They can be switched independently at every
  depth, yielding exactly (2^k) distinct counts at depth (k), with a
  greedy exact inverse. Insisting on matching their count response as
  well would erase position: the same repeated switch would yield at
  most (k+1) counts. The current digits are superincreasing, not a
  filled window.
- **[An adjacent continuation-locked pair](../chanpak_outlet_blocks_20261115/MIXED_RADIX.md).**
  Actual counts 9,971,529 and 9,971,531 have a common continuation
  profile. Every common singleton suffix preserves their difference.
  This pair matches only at its initialization; its different slot
  multipliers preclude independent repetition as a unit switch.
- **[Two exponential clocks with genuine small windows](../chanpak_backbone_20261115/MULTICONTACT_ADDENDUM.md).**
  Exact rooted transfers and an arbitrary-core readout allow a
  (2^m/3^m) clock. With fixed triangle or transitive ports one gets
  (9\cdot2^m+4\cdot3^m+\{-3,-1,1\}) for every (m\ge1), with
  an additional independent reservoir that moves both exponential
  coefficients. Preserving its first two nontrivial slot readings
  at fixed order restricts ordinary count changes to multiples of 48.
  With many ports, doubly regular cores give a growing contiguous
  third-coefficient alphabet while the two upper coefficients stay
  fixed. A lower-scale remainder is still uncontrolled.
- **[Mirrored carry decoder for the 168 small-core forms](CARRIER_WIDTHS.md).**
  Recovers all candidates for matched slopes at any width and proves
  their allowed widths are eventually periodic. The exact offsets
  still need a separate update. The independent four-bit
  three-point window survives; exhaustive bounded span matching
  through 16 finds no fourth point in the same type of universal
  affine cylinder.

## All-order limits that guide the next move

- **[Guard clock](GUARD_DEPTH.md).** A scaled 2-adic isometry reads
  depth from the count after its transient. Six towers jointly hit
  every odd multiple of 3 modulo (9\cdot2^r), at every precision
  (r), yet their absolute support has size only (O(\log X)).
  Even [changing small strong components between arbitrary guard
  wraps](GUARD_GRAMMAR.md) leaves a sublinear support bound.
- **[Arbitrary outlet modules with root changes](../chanpak_root_switch_20261115/OUTLET_MODULE.md).**
  A closed four-coordinate physical transfer separates retaining
  a root from switching into a new block. With bounded strong
  component types, (o(\log X)) root switches cannot produce
  positive density; in fact the resulting support is (X^{o(1)}).
- **[A sharp bounded-edit obstruction](../chanpak_finite_edit_20261115/FINITE_EDIT.md).**
  From bounded-count sources, almost every target up to (X) needs
  roughly (log X/\log\log X) arc reversals even when the
  source size is not bounded. At exactly two edits, only the crossing
  shape achieves cubic polylogarithmic support; its ternary
  interaction gives a literal factor-pair inverse.
- **[Whole unbounded chain families](../chanpak_backbone_20261115/BACKWARD_CHAINS.md)**
  can be locally complete modulo every power of two while
  remaining sparse. The
  [isolated-zero subclass](../chanpak_finite_edit_20261115/ISOLATED_CROSSING_CORRIDORS.md)
  has the stronger bound (o(X^{13/20})).
- **[Order-eight palette stress test](../chanpak_outlet_blocks_20261115/ORDER_EIGHT_SIEVE.md).**
  Arbitrarily many outlet SCCs of order at most eight still
  miss more than 0.5% of the targets in the fixed progression
  (1\pmod{390}) in the specified transitive-spine chart.

These results separate three properties: preserving the state
needed for continuation, making the effect of a choice depend on
its position, and placing a sufficiently wide digit range at the
prescribed absolute baseline. We can now get the first two in the
same literal tournament construction. The third remains the
proof obligation.
