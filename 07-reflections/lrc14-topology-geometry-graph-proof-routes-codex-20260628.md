# LRC14 Topology Geometry Graph Proof Routes

This pass tries to make the pictures load-bearing.  The best current view is
not "draw the lonely runner circle and hope"; it is a small atlas of pictures
with explicit payloads.

The main update from incoming HYP-3240/HYP-3241 is that the witness side has
become much cleaner.  AP and Goddyn-Wong share the same `Phi_14` core
witnesses, and covering dilations promote the witness to `Phi_{14d}`.  The
new mac-mini HYP-3242/S78 topology packet adds the matching Cech invariant:
the cap is the measure-weighted Euler characteristic of the danger-cover
nerve.  So the hard frontier is less "find a witness in the dark" and more:

```text
prove tight-locus finiteness, prove bulk positive measure,
and glue them across the Vitali core without losing endpoint/sign data.
```

The visual proof route I now trust most is:

```text
endpoint arrangement -> topes/cocircuits -> Cech components
-> finite chamber atlas -> state lift or named residual.
```

Green conductance and algebraic connectivity still matter, but as energy
labels on the finite cells rather than as the final quotient.  A conductance
graph can show why AP is the positive/even extremal face and can discharge
finite traps.  It cannot see the odd/negative side unless negative covariance,
Hermite-Biehler, or Borsuk-Ulam payload is retained.

This also clarifies how a "visual proof" might inspire nonvisual proof routes.
The cell picture becomes an oriented-matroid theorem.  The component picture
becomes a Cech/good-cover theorem.  The electrical picture becomes a
Rayleigh/Thomson inequality on chamber labels.  The root-motion picture
becomes a Lee-Yang/ear payload theorem.  The last obstruction becomes a
state-lift theorem into the forbidden `H=7` endpoint.

The assumption challenged here is runners-as-vertices.  The useful vertices
are topes, cocircuits, endpoint owners, witness strata, conductance cuts,
normal-fan chambers, PGF root events, ear payloads, and proof obligations.
The LRC predicate they preserve is "strict open witness or known equality core
or finite forbidden residual."  Raw scalar pictures are fine as diagnostics,
but they destroy exactly the endpoint/sign/core data a proof needs.

-> HYP-3243, HYP-3242, HYP-3241, HYP-3240, HYP-3238, HYP-3237, HYP-3236, HYP-3235,
HYP-3234, HYP-3233, HYP-3232, HYP-3230, HYP-3228, HYP-3227, HYP-3225,
HYP-3224, HYP-3223, HYP-3222, HYP-3220, HYP-3201, HYP-3128, HYP-3123,
HYP-3108, THM-572, T1340, LTI-340, LTT-240, OPEN-Q-108.
