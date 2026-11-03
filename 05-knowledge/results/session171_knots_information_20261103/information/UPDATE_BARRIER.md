# The promise gap during a protected update

This is a short consequence of the storage-code theorem, useful when interpreting the physical compiler as a dynamics rather than a static encoding.

**Theorem (no nontrivial one-arc update inside the promised memory).** Let C⊆F_2^M be any binary code of minimum distance delta≥2, embedded by the two-seed ordinal-sum compiler. No two distinct valid stored states are joined by a single arc reversal. Every literal one-arc schedule implementing a nontrivial update must leave the encoded code family. Nevertheless any pair of valid states can be joined through the larger constant-profile seed family, by changing their differing data arcs one at a time.

**Proof.** The distance isometry gives at least delta differing arcs between distinct codewords. Flipping the differing data arcs one by one is a path of length wt(h+g); every intermediate word still consists of the same seed alternatives and therefore has the same complete profile. ∎

For the holonomy code, delta=2^r≥4. Thus constant-profile microsteps do not certify that each intermediate state satisfies the depth-r holonomy promise. The subgroup promise belongs to the completed table, and a clock or latch marking completion is a necessary sidecar if intermediate states are interpreted as transport defects.

**Sharp midpoint ambiguity.** Let f be the indicator of an affine r-flat, a minimum-weight holonomy codeword of weight delta=2^r. Along a shortest update from T_0 to T_f, each intermediate word has support a subset A of that flat. Every nonempty proper A lies outside C_r, because its nonzero weight is below delta. When |A|=delta/2, the intermediate tournament is at labelled arc distance delta/2 from both T_0 and T_f. No unique decoder can extend the correction radius to delta/2 for every received tournament. When |A|<delta/2, T_0 is uniquely nearest; when |A|>delta/2, T_f is uniquely nearest. Indeed the distance to any codeword other than the nearer endpoint is at least delta minus the distance to that endpoint.

All these intermediate tournaments share the endpoints' entire profile. The invariant is blind both to the update and to the passage through ambiguous or invalid code states. This provides a concrete separation between a protected memory's code constraint, a physical single-flip trajectory, and a coarse invariant maintained along the trajectory.

This is not a no-go theorem for updates with a larger workspace, simultaneous multi-arc changes, relabellings, or representations whose promised intermediate family changes with the clock.
