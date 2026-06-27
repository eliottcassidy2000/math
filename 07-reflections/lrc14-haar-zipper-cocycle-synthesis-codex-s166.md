# LRC14 Haar Zipper Cocycle Synthesis

This pass moved the Haar-product story one level deeper.  HYP-2989 found the
local checkerboard, HYP-2988 classified the dyadic product interactions, and
HYP-2990 gave the abstract zipper rule.  HYP-2991 says what the load-bearing
local coordinate actually is:

```text
zeta(T) = T00 - T01 - T10 + T11.
```

The useful surprise is that this is not just another coefficient.  For the
`2 x 2` fixed-margin square, `(margins,zeta)` is complete in the audited finite
range, while margins alone collide constantly.  That makes the old warning
about row/column shadows precise: a proof quotient that keeps the margins but
forgets `zeta` has forgotten the exact direction along which the local
tournament tile can flip.

The discrepancy reading is also cleaner now.  The raw dyadic product census is
mostly orthogonal zero, and every nonzero non-atom class is sign-balanced until
the LRC labels break symmetry.  So the proof is unlikely to be a scalar
positivity argument on Haar mass alone.  It should be a cocycle-routing theorem:
each mixed sign cancels by color-compatible discrepancy, stops at an AP/GW
boundary cocircuit, hands to an owner/period/certificate clock, descends to a
smaller family, or becomes a named state-lift residual.

That gives a concrete interpretation of the zipper theorem.  A zipper is not
just "two proof methods agree."  It is a promise that every forgotten local
cocycle is either reconstructed, annihilated, or recorded as debt.  In LRC14
terms, `F7` should be the residual sector of unpaired `zeta` after all the
standard teeth have failed.

The Tournament Analysis also changed shape.  Runners, arcs, dyadic rectangles,
endpoint walls, margins, switches, and proof obligations were all plausible
vertices.  The best choice for this pass was proof carriers plus the local
`zeta` obstruction.  That made the tournament transitive:

```text
haar_zipper_cocycle > certificate_handoff_atlas > exposure_kernel_audit >
tope_cocircuit_wall > color_resonance_discrepancy >
admissible_smoothing_clock > fixed_margin_tiling_shadow >
raw_component_count_K
```

This does not prove LRC14.  It does narrow the missing theorem.  The next finite
experiment should construct the labelled packet grid and try to route every
nonzero local mixed cocycle through Z0-Z4 before declaring any F7 state-lift
debt.  If that routing succeeds on the hard banks, the summit statement becomes
much sharper: no primitive counterexample can carry an unpaired mixed Haar
cocycle without also carrying a named THM-572 state-lift obstruction.
