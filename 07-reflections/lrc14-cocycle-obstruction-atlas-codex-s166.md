# LRC14 Cocycle Obstruction Atlas

The cocycle pass changes the unit of bookkeeping.  HYP-2991 made the local
coordinate explicit: the fixed-margin square has a mixed Haar direction
`zeta(T)=T00-T01-T10+T11`, and margins alone forget it.  HYP-2994 says this is
not an isolated warning.  The current LRC14 proof stack is better read as a
cochain ledger.

The clean dictionary is:

```text
C0  packet labels, owner potentials, safe sections, Fejer centers
C1  handoff arrows, endpoint transfers, smoothing gauges, source pullbacks
C2  Haar switches, tope curls, color squares, boundary moments, state-lift faces
H2  unpaired mixed modes, no-hidden-kernel survivors, F7/THM-572 residuals
```

This is not meant as decorative topology.  It is a quotient discipline.  If a
coordinate is forgotten, the next proof step must say which of four things
happened: it was an exact coboundary, it was killed by a dual certificate or
boundary stop, it was a torsion/period class with labels attached, or it was
emitted as a named residual.

The tournament result was usefully non-linear.  The top layer is not simply
one best technique.  The three-carrier SCC says the local Haar `zeta`, the
certificate handoff atlas, and the exposure/Cech gluing audit are mutually
dependent.  A local mixed cocycle is not theorem-safe just because it is named;
it needs a handoff route and an exposure check proving that no hidden kernel
survives the cover.

The sparsely preserved labels are also the warning labels: exact scale, mixed
Haar sign, and phase period.  Those are exactly the coordinates scalar
summaries tend to erase.  So future LRC scripts should treat those fields as
first-class outputs, not after-the-fact annotations.

The most useful next computation is concrete: attach `zeta` and Haar-tile
signatures to HYP-2963 packet rows, then group the Fejer certificate workload by
identical cocycle signature.  If that compression works, it turns the abstract
zipper rule into an engine.  If it fails, the failure should not be anonymous;
it should produce an `F7` record with packet fiber, mixed cocycle, harmonic
sector, state-lift route, and failed exits.
