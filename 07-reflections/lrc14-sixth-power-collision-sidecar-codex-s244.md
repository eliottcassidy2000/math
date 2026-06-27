# LRC14 sixth-power collision sidecar

The useful distinction in the two prompt equations is arity.

```text
a^6+b^6+c^6=d^6+e^6+f^6
```

is a native six-term signed relation.  It belongs in the same room as the
support-six coimage tail, the low-height wall ledger, the octahedral current,
and the route-state closure interface.  In the bounded S244 scan, the first
primitive example through base `80` is:

```text
(3,19,22) = (10,15,23)
```

meaning equal sums of sixth powers, with gcd `1`.  Its small-modulus masks are
almost too perfect: mod `7` all six unit terms collapse to `1`; mod `9` both
sides are `0^1 1^2`; mod `13` both sides are `1^2 12^1`; mod `27` both sides
are `0^1 1^1 19^1`.  Then mod `41` keeps a genuinely different phase.  That
is exactly the LRC14 pattern: local legality is cheap, but magnitude-sensitive
phase still has to be retained.

```text
a^6+b^6=d^6+e^6
```

is different.  It is a four-term equality of squares of cubes, because
`x^6=(x^3)^2`.  The bounded scan found no nontrivial hits through base `220`,
but the more important point is structural: even if a hit appears, it does not
natively occupy the support-six layer.  To use it in a six-term LRC relation,
we must add a canceling pair and mark the packet as padded/degenerate.  That
is a sidecar, not a final obstruction.

This merges cleanly with the recent HYP-3060/HYP-3074 stack.  HYP-3060's Beal
gate says primitive multi-channel collisions need common-owner information or
named debt.  HYP-3074 says proof states must be closed under legal sidecars
before medians are trusted.  The sixth-power addition is:

```text
native_3v3_support6_collision
rank_lowered_2v2_square_cube_shadow
padded_support6_canceling_pair
sixth_power_residue_phase_mask
owner_gcd_common_factor_gate
```

The proof move I would try next is small and concrete: add these fields to a
HYP-2963 packet sample already carrying `relation_lattice_covolume`,
`low_height_wall_ledger`, `cycle_class_image_status`, and
`route_state_closure_status`.  Any native `3-vs-3` wall should be forced into
finite wall/cycle image/THM-572-F7 routing.  Any `2-vs-2` equality should be
treated as a degeneracy guard unless some independent LRC sidecar makes it
native.

Tournament Analysis used proof carriers and sidecar fields, not bases or
runners.  The carrier order was transitive:

```text
labelled_packet_sheaf
> native_three_vs_three_support6_collision
> sixth_power_residue_phase_mask
> route_state_closure_sidecar
> low_height_wall_ledger
> owner_gcd_common_factor_gate
> padded_support6_canceling_pair
> rank_lowered_two_vs_two_square_cube_shadow
> raw_equal_sixth_power_scalar
```

That is the whole moral in one line: equal sixth powers are useful only after
we remember what kind of equality they are.
