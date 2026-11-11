# Independent boundary-event audit of the two-birth matching law

**PASS.** The all-order argument in
[TRIANGLE_C2_PROOF.md](../arithmetic/TRIANGLE_C2_PROOF.md) was reviewed as a
proof, separately from its finite coefficient certificates. This audit
does not claim an arbitrary-k extension of the k2 argument.

The crucial factorization is lawful when a group retains the identities
and component locations of every boundary event, and fixes the actual
endpoint vertex precisely when its orientation variable is selected.
These data enforce each bin's first and last occupied components. Interior
local portions may be empty. An old path visits components in ordinal
order, so local covers concatenate uniquely; there is no further interblock
compatibility or sign depending on an unselected endpoint label. An
incompatible event arrangement contributes zero.

An unmarked triangle's label rotation preserves every retained condition,
including unselected boundary-event type and location. It acts freely on
the labelled cover objects. This gives an independent factor3 for every
unmarked triangle, even when that triangle is split among old path bins.

At a single selected event, the root bin has a specified endpoint. With
two available bins its nonempty segment can have1,2,3 vertices, with exactly
one valid cyclic path order in each case and a uniquely ordered remainder.
Hence the local count is3. Independent exhaustive local cover enumeration
gives1,3,7 for one,two,three available bins, and3,12,33 with no fixed event.

All five valid four-event availability walks were checked. ESES alone lacks
two single-event factors3 and requires two disjoint zero-availability
intervals. Such an interval between marked triangle blocks cannot contain
even one old singleton; its endpoint triangles must be adjacent with zero
gap. The two intervals yield disjoint matching edges. This pays the nu<=1
case without assuming the same triangle pair can be used twice.

For three marked triangles and nu0, the endpoint cases have three selected
events: the only walks lacking a two-available event have a zero interval
between marked blocks, and are impossible. With four events in four blocks,
at most one two-available position can be an unmarked block; ESES again
forces a forbidden triangle pair. With events in three blocks, the doubled
block may have start/end events at the same vertex or distinct vertices.
The proof uses only its net availability change, so neither case requires
an artificial within-triangle order. If the doubled block is in the middle,
the two outer singles must be S and E and each contributes3. At either end,
at least one of the remaining singles contributes3. This checks the
potentially delicate three-marked-block step.

Finally, the constructive upper bound supplies extensions with H41 on
each selected adjacent triangle pair and H5 on each selected single
triangle. All-sink births supply3^r, so the gcd has no prime divisor other
than3. The stated lower valuation is attained by the disjoint blocks.
For r1 and r2 unused births can be sinks, and the exponent is truncated
at zero exactly as stated. Arbitrary external singleton padding introduces
no extra triangle or matching edge.

The finite companion `audit_triangle.py` checks only the local cover
numbers and complete three/four-event walks. Run it normally and with-O;
`triangle_audit.json` records these controls and the proof-file SHA256.
The general coefficient grouping and ordinal factorization were checked
analytically above, not inferred from the finite rows.
