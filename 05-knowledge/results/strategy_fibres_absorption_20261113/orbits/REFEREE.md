# Independent review of neighbouring lanes

**ACCEPT, scoped proof review by arithmetic_orbits, November13,2026.**

## Five-atom fibre results F1--F5

Reviewed [fibres/PROOF.md](../fibres/PROOF.md), its exact implementation and
receipt. The normal-form proof stays in the nonnegative fibre: before the
maximum v is reached, either t and q are positive, or t=0 forces u and q
positive. The stated normal form is unique. The two-point fibres certify
indispensability, not just failure of one particular search strategy.

The temporary S simulation is valid for the declared batch rewrites and
adjacent atom swaps. The named helper can be placed in the T consumed by the
second rewrite and returned as a singleton, because vertex assignment within
the retained rewrite block is permitted. This is not an induced embedding or
an assertion about intermediate individual edge flips.

F5's Dickson argument is sound. After cancelling common coordinates of a
pair x,y with equal image, a componentwise minimal nonzero disjoint-support
relation exists below their surplus and deficit. Replacing its positive side
by its negative side stays nonnegative and reduces L1 distance by |a|+|b|.
The positive order row ensures both sides are nonempty and have equal total
vertex count. The argument proves existence of a finite move bank, without
claiming useful bounds on its largest move or new prime supply.

The P2 valuation inverse and common-helper cancellation are sound. They
correctly distinguish preserving a scalar target from preserving its next
consumer response. No unrestricted-context equivalence is claimed.

## Root growth note G1--G3

Reviewed [GROWTH.md](../GROWTH.md). G1 correctly uses persistence of literal old
cycles under births and the inherited deletion-kernel criterion; every positive
birth contributes a distinct cyclic vertex. It does not apply after rewiring.
The hard-core compiler uses at most two subtractions before each divide and
at most floor(log_3 h) divides, proving its stated vertex bound. It is not
silently lifted to a realizable tournament cycle inventory.

For G3, retaining the arcs of C_n puts all n vertices in one SCC; applying the
inherited strong floor and SCC product rule gives the exponential lower bound.
Requested the harmless explicit convention n>=3, which the root accepted.
The claim does not exclude compressed marked-response compilers.

These are proof and implementation reviews, not new numerical reruns of the
other agents' full finite universes. Own orbit tests are recorded separately
in receipt.json, receipt.optimized.json, and singular_receipt.json.
