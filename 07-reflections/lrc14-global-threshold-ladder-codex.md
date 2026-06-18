# LRC14 Global Threshold Ladder

Codex 2026-06-18.

The important move today was not computational bravado; it was changing the
subscript.

The old shorthand `rho*(P,E)` hid a threshold:

`rho_alpha(P,E)=meas(G_P cap {x : maxgap(E*x)>alpha})`.

Once that is explicit, the apparent contradiction in the recent frontier
disappears.  The via-max threshold `alpha=2/7` is false as a uniform floor, but
the global LRC witness threshold is only `alpha=1/7`.  These are different
objects.

## What The Ladder Says

The exact scan tested `132005` cases across consecutive clusters, structured
perforated near-AP shapes, relation-lattice tail stresses, and random
bounded-spread rows.

The summary is stark:

- `rho_{2/7}` has exact zeros (`54` in this bank).
- `rho_{4/15}` has no zeros, minimum `2/525`.
- `rho_{1/4}` has no zeros, minimum `1/140`.
- `rho_{1/7}`, `rho_{1/6}`, and `rho_{3/14}` have no zeros; their aggregate
  minimum is the `k=3` small-part floor `14249/252252`.

The explicit anti-correlation examples make the geometry visible.  The case

`P=(1,2,3)`, `E=(0,2,3,4,5,6,7,8,9,10)`

has `rho_{2/7}=0`, but `rho_{1/4}=1/140` and `rho_{4/15}=2/525`.  So the
zero is not a global-witness obstruction.  It is a boundary phenomenon at the
via-max certificate.

## Why This Helps The Denominator Route

The small-remainder language `M=j/D`, `D=14j-r`, becomes useful only after the
other runners clear.  A positive `rho_{1/7}` reservoir says a global witness
exists in the limiting phase space, but it gives no slack.  A positive
`rho_alpha` reservoir with `alpha>1/7` gives clearance margin:

`alpha - 1/7`.

At `alpha=4/15`, the margin is `13/105`.  At `alpha=1/4`, it is `3/28`.
Those are large compared with a `1/Vmax` Weyl or CRT placement error.  This
suggests a cleaner two-step proof:

1. Prove a uniform density floor at an intermediate threshold.
2. Use the slack to show finite denominator grids cannot be fully covered by
   conditional blockers.

This reframes "covering forbids small binding denominator remainder" as a
placement theorem with slack, rather than as a private pair inequality.

## The New Conjectural Lever

The critical-threshold probe computed

`sup_{x in G_P} maxgap(E*x)`.

Among the exact minima and zero examples checked, the worst was `51/182`, still
above `4/15`.  That is not a uniform theorem, but it suggests a useful split:

- First prove a **critical slack** floor:
  `sup_{x in G_P} maxgap(E*x) >= beta > 1/7`.
- Then prove an actual **density** floor below that `beta`.

The density step is harder, but now it has a visible target below the false
`2/7` boundary.

## Assumption Challenge

I deliberately did not use runners as tournament vertices.  Runner-level
blocking has been transitive in too many of the hard families.  Here the
vertices were proof obligations: threshold levels.  That preserves the
implication chain and shows exactly which target collapses (`2/7`) and which
ones remain plausible (`1/4`, `4/15`).  It destroys the local interval
topology, so if this route stalls the next tournament should use safe
components, denominator events, or short relation-lattice patterns.

The lesson is small but sharp: the proof probably does not want the strongest
sufficient certificate.  It wants the weakest certificate with enough slack to
survive discretization.
