# LRC14 Packet Migration Gauntlet

HYP-2955 is the first check that made the recent picture feel less like a pile
of analogies and more like one proof object.

The object is not a runner tournament and not the scalar `M(S)`.  The object is
the labelled packet:

```text
qdiv gate
regular-open Haar witness U(S)
finite boundary owner skeleton
C27 transfer label
K33/state-lift flag
wide/decorrelation escape class
```

The finite result is simple:

```text
one-swap add<=420: only GW is boundary-only
two-swap add<=60: no boundary-only rows
three-swap add<=30: no boundary-only rows
```

There are also no covered hard rows.  Every non-AP/GW hard row in this gauntlet
has positive strict Haar interior.

## Missing Picture

The prior route sometimes phrased K33 as the endpoint.  This run suggests a
more careful order:

```text
first ask whether the packet migrates to positive Haar interior;
only if it refuses to migrate should K33/HYP-2908 be asked to close it.
```

For the known near-miss rows, K33 is a label on a positive-open front.  That
already proves those rows are safe.  The state-lift theorem is reserved for a
hypothetical zero-interior packet that keeps the K33/nonunit labels without
opening a witness interval.

So the proof may split as:

```text
AP/GW source core          zero regular-open witness
finite AP-neighborhood    forced migration to open front
wide/unbounded packets    decorrelation migration
non-migrating labelled    HYP-2908/THM-572 state lift
```

This connects the local AP/GW census with the older wide gK8/decorrelation
route: both are migration theorems.  One is topological near the boundary; the
other is analytic in the tail.

The concurrent HYP-2952/HYP-2953 work makes this cleaner.  HYP-2952 supplies
the derived boundary tournament filter for AP/GW-kind atoms; HYP-2953 supplies
the source-spectrum pullback that can remember where a packet's deepest source
node lives.  HYP-2955 is the observed local transition rule between them:
outside AP/GW, the source address migrates to positive Haar interior in every
tested AP-neighborhood packet.

## New Hypotheses

1. **Boundary-source collapse.**  In the reduced bounded AP-facing atlas, the
   only zero-regular-open source packets are AP and Goddyn-Wong.

2. **Open-front migration.**  Every finite non-AP/GW AP mutation with qdiv>=14
   has a positive Haar-open front unless it carries an explicit non-migrating
   state-lift obstruction.

3. **K33 reserve principle.**  K33 should not be used as the first proof
   endpoint when a positive Haar front is already present.  It is the endpoint
   only for non-migrating zero-interior packets.

4. **Wide analogue.**  The unbounded proof should imitate the local result:
   show that any packet outside the finite AP source atlas loses its source-core
   address by decorrelation before it can become a counterexample.

## Guardrail

The qdiv/Haar scout does not compute exact `M(S)` for every positive row.  That
is intentional: positive strict Haar interior is already stronger than needed
for LRC safety at threshold `1/14`.  Exact `M` remains necessary when ranking
low-frontier values, but not when discharging a row by open witness mass.

The quotient destroys unbounded row identity and analytic tail behavior.  That
is the remaining theorem, not something to hide.
