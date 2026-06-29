# LRC14 AP84 Corridor-Splice Certificate Reflection

HYP-3462 changes the AP84 status from "three sidecars plus a splice TODO" to a
single routed bridge packet.

This is deliberately separate from HYP-3460.  HYP-3460 is the noncanonical
phase-branch color pullback for the random031 route; HYP-3462 is the canonical
AP84 carrier/splice closure.  The two packets should be composed downstream,
not conflated.

The useful correction is that HYP-3431 is a branch-union carrier.  Branch 0 and
branch 1 each have a small overlap inside the opposite corridor, but their
union is exactly `[8/49,6/35] union [29/35,41/49]`.  That is the statement the
AP endpoint and floor-count lemmas need.

The one-branch rescue split is also now explicit on the AP tower: `m=1` is the
only checked rank-`6` row, with core `(3,5,7,9,11,13)`, while every
`m=2..70` row drops to the rank-`5` core `(5,7,9,11,13)`.  This lines up with
the sidecar split: HYP-3457 handles `m=1..4`, HYP-3454 handles the endpoint
interval from `m=5`, and HYP-3456 handles the boundary count.

The challenged assumption was that "import HYP-3431" meant a pure one-branch
corridor.  It does not.  The preserved predicate is AP84 branch-union
survival; the destroyed information is arbitrary non-AP component geometry.
So the next task should not reopen AP84.  It should feed this closed splice
into HYP-3460/HYP-3453/HYP-3451 and isolate any remaining noncanonical gluing
debt, especially HYP-3455.
