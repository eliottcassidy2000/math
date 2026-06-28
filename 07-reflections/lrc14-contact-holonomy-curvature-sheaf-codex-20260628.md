# LRC as a certificate sheaf with curvature

The useful new frame is not another scalar invariant.  It is a way to line up
the existing invariants without pretending they commute automatically.

Read the current route as a sheaf of certificates over the LRC packet space:

```text
global index          HYP-3246
calibrated boundary   unit equioscillation / Fejer autocorrelation
local curvature       shell-lag commutator
connection repair     zeta_7 contact holonomy
```

HYP-3247 already showed that lag space forgets endpoint placement.  The exact
new scout sharpens that: over the 62 mixed lag+residue fibers, the first
cyclotomic contact moment

```text
sum zeta_7^j over non-unit gap positions j
```

repairs shell magic completely.  Position sum, min/max position, support size,
gap multiset, and ordered gap values all fail somewhere.  The shell-lag
commutator is therefore not only positional; at this interface it is
cyclotomic positional.

That is why the sheaf language is useful.  HYP-3246 wants a global nonzero
index.  HYP-3245 gives a lag/autocorrelation projection.  HYP-3247 says that
projection has curvature.  HYP-3267 identifies the small holonomy coordinate
that kills that curvature on the bounded bank.

The proof target becomes concrete:

```text
every primitive residual packet has zero quotient curvature,
or its curvature is killed by zeta_7 contact holonomy,
or it lifts to an endpoint arrangement / finite chamber,
or it becomes named residual debt.
```

The caveat matters.  Holonomy alone is not a legal global quotient: empty and
full contact support both have zero first moment.  So the connection coordinate
does not replace endpoint cells.  It tells us which endpoint data must be
kept when lag space talks to the shell/index packet.

The newest HYP-3249 warning makes this sharper.  A naive odd
Borsuk-Ulam map can vanish at a runner collision with the observer, which is
the opposite of a lonely witness.  The HYP-3267 answer is not that holonomy
proves the global theorem; it says the proof map must keep cover-hole or
endpoint-cell data whenever it transports the index.  Losing that coordinate is
exactly how a topological zero can become the wrong zero.

The incoming frontier tests put a sharper boundary around this idea.  HYP-3250
turns the global proof shape into finite tight-locus rigidity plus a uniform
margin away from that locus.  HYP-3251/HYP-3252 then make the honest correction:
the index is an excellent description of the AP saddle, but the S-dependent
floor still has to do the proof work.  In that picture, contact holonomy is a
local diagnostic for the finite side of the split: if lag/residue data tries to
collapse two endpoint chambers with different shell magic, the holonomy names
the missing coordinate before the packet is allowed to descend to a margin
argument.

The same paragraph gives a natural next experiment: if the floor wants a
`Q(sqrt(-7))` basis, do not feed it raw lag/residue fibers.  Feed it fibers
already repaired by the `zeta_7` contact holonomy, then ask whether the signed
floor terms become positive or finite-chamber exact in that basis.

The creative payoff is that the topological index and the finite shell
commutator now look like the same kind of object at two scales: global degree
versus local holonomy.  Proving LRC14 from here means proving that local
curvature cannot hide from the holonomy/endpoint lift while the global odd
index remains nonzero.
