# LRC14 Universal Scale Invariance Recursion Ledger

The instruction was to assume scale invariance holds universally and then use
the repo history to see the recursive patterns we have not precisely named.
That assumption does not make the problem easier by itself.  It makes the
bookkeeping honest.

The exact invariant is already in the repo in several forms: Lean's
`lonely_scale`, THM-407's `M(cS)=M(S)`, THM-532's scale-invariance of
`meas(S7)`, and the covering-bound correction that dilated AP/GW rows sit on
the floor.  The point is not merely "work WLOG primitive."  The point is that
every live proof route is a projective route.  Absolute scale is a gauge, and
each proof stage must say which coordinate survives after that gauge is
divided out.

The route has sharpened in a consistent direction:

```text
scalar route
-> proof currencies
-> scale and shell quotients
-> relation-height split
-> primitive covering residual
-> labelled packet sheaf
-> legal sidecar median
-> protected branch graph
-> Fejer / Gram / Toeplitz / Green certificate vector.
```

That is the gradual sharpening.  The object did not become a better scalar.
It became a better ledger.

The hidden recursion is:

```text
normalize scale
split the first coordinate scale cannot kill
attach a sidecar for the coordinate a quotient would forget
discharge the easy branch
recurse on the remaining primitive packet.
```

This appears at every level.  THM-407 normalizes shells but leaves the
`gcd(a,27)=1,3,9` strata.  THM-573 removes the 14-free and `>=7` multiples-of-7
branches but leaves `<=6` multiples-of-7 primitive covering rows.  THM-532
turns seven-sector mass into a high-relation-height certificate plus a finite
low-height AP-rich residual.  HYP-2963 and the HYP-3030 stack do the same with
route, topology, owner, primitive deck, Haar zeta, and residual capacitor
fields.  HYP-3229 does it analytically: order-2 Fejer/Gram positivity is
usable, but order-3 triple-overlap constants are named debt.

So the next proof interface should carry scale fields explicitly:

```text
primitive_scale_gcd
scale_orbit_representative
dilation_fixed_boundary
nonunit_residue_stratum
scale_destroyed_payload
renormalization_depth
```

These are not decorative labels.  They prevent the old error where a quotient
looks clean only because it hid the coordinate that proves the route.  A
nonprimitive tie is then recognized as a gauge copy, not as a new extremizer.
A nonunit residue stratum is recognized as real residual data, not as noise.

The tournament-analysis lesson is also sharper.  The right vertices are proof
states and scale fibers, not runners by default.  Runners, gaps, cover arcs,
Fourier modes, endpoint owners, matroid circuits, and triple overlaps are all
available vertex sets, but the chosen quotient must preserve the LRC predicate
or state exactly what it destroys.  For this synthesis, the preserved
predicate is `M(S)>=1/14` modulo nonzero scale; the destroyed data are
absolute speed size, reset-calendar complexity, finite product-bound size, and
some endpoint event order.

The proof target is now cleaner:

1. Build the scale-normal packet ledger for the post-THM-573 residual.
2. Prove scale-fiber exactness for every HYP-2963 quotient.
3. Turn S75/HYP-3229 order-3 triple-overlap debt into finite LP rows.
4. Merge high relation height with the HYP-3215/Rosenfeld wide branch.
5. Keep the induction-base flag separate from the bounded certificate.

This does not close LRC14.  It removes a layer of ambiguity about what the
route is doing.  The repo's route has sharpened from "find the invariant" to
"prove the scale-normal recursion terminates with a certificate or a named
residual."

-> HYP-3230, HYP-3229, HYP-3215, HYP-3214, HYP-3205, HYP-3162, HYP-2963,
THM-573, THM-532, THM-407, LTI-330, LTT-230, T1330, OPEN-Q-108.
