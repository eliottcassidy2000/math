# LRC14 Phase-Branch Color Pullback

The useful reconnection is that the two "color" languages are not competing.
The old distance-graph coloring formulation says a witness time `t` is a
regular circular `14`-coloring.  The HYP-2593 CRT layer refines this to a phase
color `b=a mod 14` for `t=a/(14V)`.  The two-adic branch layer then pulls the
same witness through `u=2t mod 1` and records whether `t` came from the first
or second half of the circle.

That sounds formal, but it changes the random031 story.  HYP-3455 isolated
the noncanonical rank-`6` obstruction as two max-delta mirror gates.  When
`random_covering_031` is read as an S3 row with `V=173`, the colored CRT layer
finds `282` actual witnesses, exact mirror symmetry in phase and branch
counts, and zero hits on those two max-delta gates.  The same hard components
are touched only through branch-opposite, lower-delta gates.

So the max-delta mirror pair is a continuous branch-gluing proof debt, not a
finite phase-grid obstruction.  That is good news: it means the coloring
route can become a bypass test for scary local gates.  A gate that receives no
compatible phase-grid hits should not be promoted to a global obstruction
until it survives colored resonance cancellation, low-rank component escape,
endpoint-spine/wall lift, owner-current, two-adic descent, and signed-SPEC
sidecars.

The AP84 comparison is the opposite and helps calibrate the distinction.  At
`m=1`, phase colors hit the finite transient edge-singleton gates.  At `m=5`,
they hit the rank-one endpoint-phase gates and the actual deficit is already
negative.  Thus canonical AP gates are real phase-grid gates, while the
random031 max-delta gates are bypassed by the phase grid.

Post-rebase HYP-3459 makes the AP84 side stricter: the canonical splice needs
the full color-packet legality test, not just a gate color or a floor word.
HYP-3460 supplies the complementary noncanonical test: when that colored packet
is pulled through the branch atlas, do the phase-grid witnesses actually hit
the scary local gates?

The next proof task should be stated as a pullback lemma:

```text
If a branch-colored gate obstruction has zero compatible phase-grid hits,
then it must discharge through colored resonance cancellation, low-rank
component escape, endpoint/wall sidecars, owner-current imbalance, two-adic
descent, or signed-SPEC debt.
```

That lemma would merge the S359/S363 regular-coloring picture, HYP-2593/HYP-2595
colored discrepancy, HYP-2991 Haar cocycle cancellation, and HYP-3455
gate-gluing debt into one proof interface.
