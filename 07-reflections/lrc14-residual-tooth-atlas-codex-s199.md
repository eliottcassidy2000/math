# LRC14 Residual Tooth Atlas - S199 Reflection

The useful shift in this pass is that the HYP-3028 residuals stopped being a
list to fear.  The coarse ET+unit gate leaves `15` route-mixed fibers, but
HYP-3030 already showed all `38` packets in them are open.  S199 asks the
right follow-up: what is the first non-route tooth that explains each open
route collision?

Incoming S198/HYP-3034 and S197/HYP-3033 are the immediate predecessors:
HYP-3034 lifts the AP/GW equality front to explicit owner-essential closed
arc-boundary representatives, while HYP-3033 schedules the residual q-witness
versus covering certificates by topology-bucket/unit-scale teeth without route
labels.  This S199 pass asks for the first local owner-strip tooth underneath
that scheduler.

The answer is cleaner than expected.  Arc topology separates `13` of the
`15` fibers.  The two exceptions are same-topology collisions, and the coarse
largest-safe-component stalk separates both.  Thus every first repair is an
`owner_strip` in the HYP-3031 dictionary.  Exact stalk, magnitude, and explicit
q/covering certificates remain complete backups, but they are no longer the
first proof object.

This changes the proof obligation.  We do not need a theorem that every
coarse residual has a bespoke route label.  We need a residual-tooth manifest:
after status convergence, attach compact arc topology; if that is not enough,
attach the coarse safe-component stalk.  Route labels then become bookkeeping,
not structural risk.

Assumption challenged: the vertices were not runners, route labels, or raw
fibers.  They were repair teeth inside a fixed coarse residual fiber.  That
preserves the LRC predicate that matters here - strict-open status - while
destroying exact route identity until the first legal tooth is reattached.

Next step: add `first_tooth` and `residual_tooth_class` to HYP-2963 sidecars,
then prove the two owner-strip descent lemmas separately: one for the 13
topology-pure fibers, and one for the two stalk-only fibers.
