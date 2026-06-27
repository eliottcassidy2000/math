# LRC14 Arc-Cech Nerve Carrier - S174 Reflection

This pass tries to make the geometric carrier closer to the actual LRC14
covering predicate.

The old endpoint/topes language was good at saying "strict safe interval" or
"boundary equality wall."  The old runner nerve and sequence-shadow languages
were good at naming rows or proof routes.  Neither directly stored the exact
topology of the danger cover.

The closed arc-Cech nerve does.

## What Changed

For each row, the script builds danger endpoints at threshold `1/14`, samples
open cells, then records two nerves:

```text
open arc nerve:   only arcs active on open cells
closed arc nerve: open-cell facets plus endpoint equality facets
```

That distinction matters.  AP and Goddyn-Wong look disconnected if one sees
only open cells:

```text
open_arc_betti=(6,0)
```

But the endpoint equality facets glue those six pieces into the full circular
cover:

```text
closed_arc_betti=(1,1)
```

So the equality atom is not "six independent boundary accidents."  It is a
single closed good-cover cycle whose open pieces touch only at cocircuits.

## Concrete Output

Named rows:

```text
AP, GW: closed arc beta1 = 1, safe_mu = 0
K33 12->36: closed arc beta1 = 0, safe_mu = 1/1260
petal 10->20: closed arc beta1 = 0, safe_mu = 1/980
petal 13->26: closed arc beta1 = 0, safe_mu = 1/182
covering 12->84: closed arc beta1 = 0, safe_mu = 563/105105
fibbinary first13: closed arc beta1 = 0, safe_mu = 66077/399840
Moser first13: closed arc beta1 = 0, safe_mu = 4264747/40348854
```

AP and GW also have six boundary-safe points with owner-pair sums all zero
modulo `14`.  That recovers the taut-current law, but now as boundary facets
of a closed Cech cycle.

The one-swap AP scan through `add<=160` gives exactly one zero-open row:

```text
drop=12 add=24
```

The smallest positive row is the K33 near miss:

```text
drop=12 add=36, safe_mu=1/1260
```

This confirms the known frontier, but with a cleaner object above it.

## Assumption Challenge

Vertices considered:

```text
runners
individual danger arcs
endpoint cells / topes
boundary cocircuits
safe components
nerve homology classes
runner quotient defects
Fejer atoms
automaton states
proof obligations
```

The chosen proof-carrier tournament uses proof carriers as vertices.  Inside
the row audit, the actual nerve vertices are individual danger arcs, not
runners.  This preserves the circle-cover predicate and destroys the runner
collapse.  The destruction is recorded as `runner_quotient_betti_defect`.

## Quotient Lesson

The runner-level nerve is not false, but it is a quotient.  In the named audit
it often has connected closed runner nerve even when the individual arc nerve
has several closed components.  It can also erase the cover cycle by merging
many disjoint arcs owned by one speed.

So any runner-level proof must carry a side channel:

```text
closed_arc_cech_beta
runner_quotient_betti_defect
boundary_cocircuit_facet_word
```

Without those fields, it is not a theorem-safe quotient.

## Proof Direction

The new theorem target is:

```text
Every primitive zero-open LRC14 packet either has the AP/GW closed arc-Cech
cycle and zero owner-current boundary law, exits through an existing K33 /
Fejer / Ramanujan / Haar / state-lift certificate, or becomes the first genuine
F7 good-cover quotient defect.
```

This also reorders the older proof carriers.  Automata, Moser/fibbinary,
Fermat-Catalan, and sequence shadows should be read as sidecars over the exact
arc topology, not as substitutes for it.  The proof summit should make the
closed cover nerve primary, then let dual certificates and finite-state labels
discharge the non-cover cases.

## Next Pulls

1. Add arc-Cech fields to the HYP-2963 packet bank or a sidecar manifest.
2. Run the closed nerve audit over the full HYP-2963 bank.
3. Build a family theorem for the K33 and petal beta1-zero exits.
4. Compare Fejer certificate degree against `runner_quotient_betti_defect`.
5. Define F7 as a good-cover quotient-defect class before adding more scalar
   filters.
