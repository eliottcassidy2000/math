# LRC14 random031 owner-boundary persistence

Reserved as the owner-cobordism continuation of HYP-3494.

Known data before the executable pass:

- The max-delta seam on components `(43,54)` carries all seven owners
  `(23,45,93,113,147,169,173)`.
- The lower-delta pure bypass on the same components carries only
  `(23,93,113)` and is hit by exactly `12` phase witnesses.
- The seam-only word `(45,147,169,173)` must be treated as boundary debt, not
  as ordinary phase-flow owners.

Assumption challenge: the tournament vertices for this packet should not be
runners, arcs, or raw witnesses.  The first candidate vertex set is quotient
sidecars: seam owner word, bypass owner word, hard-component owner map,
relative-H1 boundary class, owner-current cobordism matrix, dead-island owner
union, global flow-owner union, delta route, component pair, and raw bypass
count.  The preserved predicate is exact reconstruction of the pure bypass's
owner-boundary debt.  The information destroyed by the lossy vertices is the
seam-only owner word.

The executable pass should make the loss visible as a quotient-price matrix.
