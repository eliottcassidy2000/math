---
source: codex-2026-06-01-S549
status: reflection + hypothesis claim
tags:
  - lonely-runner
  - conditional-clearance
  - cascade
  - tournament-analysis
  - transitivity
  - wedge-debt
---

# LRC conditional clearances and the hidden transitive wedge

The useful move in the prompt is to keep the S548 product formula but stop
thinking of its factors as naked probabilities.  A cascade is a product of
conditional clearances:

```text
F_k = times surviving runners 1..k,
P_k = mu(F_k) / mu(F_{k-1}),
mu(F_m) = prod_k P_k.
```

For LRC, `F_m` is the observer-lonely set.  This already connects global
spread to local emptiness: each runner clears conditionally, and the product is
the simultaneous clearance.

But the tournament reminder adds a second payload to every factor.

## The hidden fact

The usual transitivity rule is:

```text
X -> Y and Y -> Z imply X -> Z.
```

The edge-local companion is:

```text
if X -> Y, then not (Z -> X and Y -> Z).
```

That is just "no directed triangle" seen from one accepted edge, but it is the
right operational form.  Once `X -> Y` is admitted, no future vertex may sit in
the backward wedge

```text
Y -> Z -> X.
```

So transitivity does not merely close paths forward.  It also clears forbidden
wedge positions behind every accepted edge.

For a tournament `T`, define

```text
W_T(X,Y) = { Z : Z -> X and Y -> Z }.
```

Then

```text
T transitive  <=>  W_T(X,Y)=empty for every edge X -> Y,
sum_edges |W_T(X,Y)| = 3*c_3(T).
```

THM-395 now records this as a canon theorem.  This packages the hidden fact as
a measurable obstruction.  The first directed 3-cycle is exactly the first
unpaid backward wedge.

The later S545/S545o notes use complementary language: HYP-2041 calls the same
edge-local obstruction a no-return or 3-term resonance debt, while HYP-2042
marks the limitation that the Helly-3 layer is necessary but not sufficient for
the full LRC cascade.  So the theorem is the exact transitive ledger; the LRC
proof still needs the correct order-`n` lift.

## LRC reading

In the observer-source tournament, a clearance is an edge saying the observer
has beaten a danger obligation for one runner/event.  In endpoint-pressure
language, a clearance is an owner or wall event being peeled.  In either lift,
the accepted edge should export no-wedge obligations to the remaining objects.

So the cascade factor should be mentally refined from

```text
P(next runner safe | prefix safe)
```

to

```text
P(next runner safe and no new backward wedge | prefix safe).
```

This is not a request to force the LRC movie itself to be transitive.  The LRC
trienerment is allowed to have ties and cycles.  The point is that any proof
route that uses a transitive or hierarchical cascade also gets the no-wedge
propagation for free, and should account for it as part of the dependence term.

The S548 last-runner bottleneck becomes cleaner in this language.  Early
clearance factors look like the outside credit `(n-2)/n`.  The AP/regular row
does not lose the battle gradually; it aligns the backward-wedge debt so that
the last runner or compact wall pays almost all of it.  Generic/non-AP rows
decorrelate that debt across the cascade, leaving positive open clearance.

## Which vertices?

The dangerous habit would be to put runners in the tournament by reflex.  For
this particular idea I considered:

```text
runners,
gaps,
fixed circle sections,
section boundaries,
danger arcs,
clearance obligations,
wall-crossing events,
endpoint owners,
residues,
p-adic zero branches,
cover arcs,
Fourier modes,
Gabor zero columns,
matroid circuits,
proof obligations.
```

Runners alone are probably too lossy: many LRC runner-level gauges collapse to
transitive shadows.  Arcs alone also risk forgetting the observer target.  The
more promising vertices are clearance obligations, endpoint owners, or
wall-crossing events, optionally decorated by the p-adic/Gabor labels already
used in HYP-2036 and HYP-2032.

Predicate preserved:

```text
observer tie-degree zero / observer-source target / zero danger occupancy.
```

Information destroyed:

```text
exact gap lengths, internal safe-sector positions, and arithmetic resonance
labels unless retained as decorations.
```

Challenged assumption:

```text
the product cascade is scalar.
```

The better object is scalar product plus edge-local wedge ledger.

After THM-395, the exact theorem-level content is no longer speculative:
transitive clearance is equivalent to zero backward-wedge debt.  The open LRC
choice is the lift in which that debt is certificate-bearing rather than only
descriptive.

## Proof shape this suggests

A future proof route could try to build a finite lifted tournament at each
cascade prefix:

```text
vertices: current clearance obligations and remaining wall/owner events
edge: forced-before / clears-before / protects-before
statistics: P_k, backward-wedge count, c_3, SCCs, HP count
```

The S550 verifier covers the exact base layer of this plan: it exhaustively
checks the backward-wedge identity for all tournaments on up to six vertices
and prints exact `P_k` rows for AP, generic, and lacunary LRC cascades.

If the lifted tournament stays transitive after enough clearances, the no-wedge
rule should peel the remaining endpoint core.  If it develops cycles, those
cycles name the exact dependence debt instead of leaving it as a mysterious
epsilon in the product.  Either outcome is useful: transitive means clearance
propagates; cyclic means the script has found the real obstruction layer.

The proof obligation can therefore be written as a cascade of ordinary
conditional clearances multiplied by THM-395 anti-wedge clearances.  A zero in
the second ledger is not noise; it is exactly a directed triangle in the lifted
object.

The slogan I would keep:

```text
conditional clearance is multiplicative;
transitive clearance is also anti-wedge.
```

For LRC, the product asks whether all danger arcs clear together.  The hidden
transitivity fact asks whether a cleared relation has accidentally created the
only shape that can block a hierarchical proof: a third object wrapping behind
it.
