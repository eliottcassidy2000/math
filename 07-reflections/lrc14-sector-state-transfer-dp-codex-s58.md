# LRC14 Sector-State Transfer DP - codex S58

The dynamic-programming lens changed the object from a row to a transfer
kernel.  THM-554 is the useful exact identity: when a speed is inserted, each
wall atom's missed-sector state either stays fixed or loses exactly one inner
sector.  Therefore the scalar increment in `p0` is exactly the mass of
one-missed-sector atoms that the new speed lands in.

This is the LRC analogue of the half-tiling address quotient.  A final row has
one `p0`, but the proof state depends on the insertion address order.  The
script shows that this order dependence is not cosmetic: greedy minimum-support
orders beat decreasing orders by a large support/area gap, and the dyadic block
specifically wants the dyadic-tower order.

The proof route I would pursue next is a transfer inequality.  At each prefix,
bound the one-missed-sector landing mass by either a finite address template
or a decorrelation estimate.  Low state support should trigger the known
finite atlases: AP, dyadic tower, AP-triple/cube-root phase, Ruzsa packet, and
residue-private ledgers.  High support should be the Weyl/BV branch: many
states mean the new speed is trying to hit a dispersed target and should pay a
decorrelation tax.

This also clarifies what HYP-2683 was seeing.  Coarse `state_mass` separated
final true-wide risk because it was already a shadow of a transition system.
The next version should classify rows by the full prefix transfer word, then
prove that each transition either belongs to a finite low-state pocket or has
small landing mass.

No LRC14 proof is finished here.  The progress is an exact DP recurrence and a
sharper proof obligation that avoids collapsing the address too early.
