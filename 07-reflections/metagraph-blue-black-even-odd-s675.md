# Metagraph Blue/Black Even-Odd Audit, S675

The answer is subtler and better than the first sentence.

If "the black part is an even graph" means the simple graph you see after
aggregating all the lines, then no: it already fails in the small audited
merged metagraphs.  In the SC-type convention, black is `SC<->NS` where `SC`
means self-converse, and it has nonzero simple boundary at `n=4,5,6`.  In the
explorer convention, pure-black complement edges also have nonzero simple
boundary at `n=4,5,6`.

But the observation has a stronger survivor.  If "black lines" means the
parallel complement lines counted with multiplicity mod `2`, then the
explorer pure-black layer has zero weighted boundary in every tested case
through `n=6`.  That feels like the real theorem candidate.  The visible graph
is the quotient shadow; the line-count chain is the actual parity object.

The blue side gives the warning.  Explorer pure-blue looks all-odd at
`n=3,4,5`, but that pattern breaks at `n=6`.  So blue is not intrinsically
the odd graph.  It is a boundary-producing color whose small cases happen to
put every incident vertex in the boundary.  At `n=6`, some blue degrees become
even and the clean all-odd story collapses.

The right language is chain-level:

```text
color layer = GF(2) 1-chain
even graph  = boundary zero
odd graph   = boundary packet present
```

This also sharpens the even/odd-number analogy.  Even numbers are closed
parity classes: no exposed endpoint.  Odd numbers are address defects.  In a
finite graph, odd-degree defects come in pairs or larger even packets by
handshaking.  That is the metagraph version of "odd carries an address and
doubling kills it."

This ties back cleanly to HYP-2245.  The tiling cube has honest side choices.
The quotient forgets address data.  What the quotient cannot remember as a
side choice returns as boundary.  So the boundary vector is a leakage detector:
where the boundary is nonzero, the quotient is missing an address coordinate.

After rebasing over S674b, this also rhymes with the signed LRC and
unit-distance trienerment addendum.  There, the observer shadow or binary
tournament edge is too coarse unless the pair-clock sign address or the
`S/U/L` equality layer is retained.  Here, the visible blue/black edge is too
coarse unless line multiplicity and boundary are retained.  Same move:

```text
keep the third/address coordinate, then collapse only after measuring boundary.
```

The Royle-even issue matters too.  The repo's equinumerosity theorem is not
about degree-even graphs.  Degree-even graph iso counts are `2,3,7,16` for
`n=3..6`, while tournament/Royle-even counts are `2,4,12,56`.  So the even
graphs that are "key to the structure" probably are not the literal degree-
even subgraphs drawn in black.  They are more likely a Royle-even target
equipped with the missing line-address side channel.

The next move I like:

```text
Find a functor F:
  tournament class + retained complement-line address -> Royle-even graph

such that:
  boundary(color layer) is preserved or killed predictably,
  H / odd-cycle packets survive as side channels,
  and the black weighted chain maps to a closed cycle.
```

That would turn "tournaments are equinumerous with even graphs" from a count
coincidence into a structure transfer.

Tournament Analysis for this session used proof routes as vertices rather
than tournament vertices or arcs.  The pairwise observable was whether a route
preserves the color-layer boundary.  The selected Hamiltonian path was:

```text
boundary_vector
> convention_split
> line_weight_parity
> royle_even_property
> degree_even_cycle_space
> raw_color_claim
```

The challenged assumption was that the visible black edge set is the even
object.  The replacement is that the weighted black line-chain may be the even
object, while the visible black graph is its address-forgetting quotient.
