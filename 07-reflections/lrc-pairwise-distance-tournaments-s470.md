# LRC Pairwise Distance Tournaments and Two-Neighbor Pressure

**Session:** codex-2026-05-31-S470

The repo already had a boundary-side tournament bridge:

```text
LRC endpoint protection ~= tournament good-cut protection on a wrapped path.
```

This session adds a time-slice bridge.  Instead of only asking whether the
stationary origin is far from every moving runner, include the stationary
runner and look at all pairwise distances among the `n` runners.

## Two Tournament Objects

The first object is geometric and direct.  At a time `t`, orient each pair
`i,j` by the open semicircle:

```text
i -> j iff j lies clockwise from i by less than 1/2.
```

Collisions and antipodal pairs are missing arcs.  So every time slice gives an
incomplete circular tournament.  This is useful because tight even-denominator
unit skeletons literally become incomplete tournaments:

```text
n=14, t=1/14: 7 antipodal ties
n=16, t=1/16: 8 antipodal ties
```

The odd baseline `n=15, t=1/15` has no antipodal ties.  This makes the
even-composite anomaly visible as tournament-completion debt, not just as a
scalar denominator fact.

The second object uses the user's suggested extra distance data.  For each
runner `i`, keep its two nearest runners with multiplicity.  For an unordered
pair `{i,j}`, define the relief

```text
relief_i(j) = nearest distance from i after deleting j - nearest distance from i.
```

Then orient

```text
j -> i
```

when `relief_i(j) > relief_j(i)`.  In words: the stronger blocker points to
the more blocked runner.  Ties stay missing.

This is the blocker-pressure graph.  It is not just a nearest-neighbor graph;
it asks whether a nearest neighbor is irreplaceable.  That question cannot be
answered without the second-nearest runner.

## What The Probe Found

The script is:

```text
04-computation/lrc_pairwise_tournament_s470.py
```

and the stored output is:

```text
05-knowledge/results/lrc_pairwise_tournament_s470.out
```

It samples exact endpoint and gap times for:

```text
initial n=14
n14 seven-ladder
n14 S380 gate ladder
initial n=15
initial n=16
```

The strongest fact is negative but useful:

```text
all selected strict pressure graphs have largest SCC = 1
all selected strict pressure graphs have 0 directed 3-cycles
```

So the standard near-misses do not have a mobile blocker-pressure core.  They
look peelable.

For the seven-ladder origin gap midpoint:

```text
t = 3053/25872
origin d1,d2 = 29/24, 27/22 thresholds
pressure strict/tie arcs = 7/84
pressure SCC = 1
```

For the S380 gate ladder origin gap midpoint:

```text
t = 4339/51744
origin d1,d2 = 4339/3696, 29/24 thresholds
pressure strict/tie arcs = 6/85
pressure SCC = 1
```

Even at the best sampled S380 origin time:

```text
t = 6749/51744
origin d1,d2 = 29/24, 27/22 thresholds
pressure SCC = 1
```

The pairwise pressure still does not close.

## Why This Matters

THM-365 says an LRC counterexample needs labelled endpoint cycles.  The new
pressure graph suggests a mobile analogue:

```text
no pressure cycle -> peel a private runner/blocker relation
pressure cycle    -> check arithmetic labels, as with endpoint cycles
```

This gives a possible extra certificate layer between scalar gap size and the
full endpoint incidence matrix.  It is cheaper than endpoint protection but
more structured than min-distance alone.

## Reframed Proof Route

For a branch in the n=14 proof search:

1. Build the exact endpoint rows as before.
2. Sample the active endpoint/gap times for the branch.
3. Compute the blocker-pressure graph on `{0} union V`.
4. If the strict pressure graph is acyclic, try to turn a source/sink into a
   private endpoint row or a Hall-dual leaf.
5. If a strict SCC survives, label its arcs by speeds and endpoint residues
   and try to kill it with the THM-365 labelled-cycle obstruction.

This pairs well with the `14`-gate fan tax.  HYP-1920 says the local gate row
forces the odd fan plus one even bridge.  The pressure graph asks whether that
bridge creates a real cyclic blocker core or just a peelable local moat.

## Caveats

Pairwise mobile loneliness is not the same as the reduced origin LRC condition
for a fixed speed set.  For example, the S380 row at `t=1/2` has a very lonely
moving runner while the origin collides with the gate-heavy pile.  This is
still meaningful, but it belongs to translated difference sets rather than to
the origin row directly.

Also, the current pressure graph is unlabeled except for speeds.  To become a
proof tool it needs the endpoint labels:

```text
owner speed
protector speed
endpoint residue
strict slack
```

That is the bridge back to THM-365.

## Next Computations

1. Add pressure metrics to the n=14 branch-and-bound ledger.
2. For each forced odd-fan plus even-bridge branch, compute endpoint leaves
   and pressure leaves together.
3. Search deliberately for speed sets whose critical sampled pressure graph
   has a nontrivial SCC; if none appear, formulate a pressure-peeling lemma.
4. If SCCs appear, label their arcs and compare them with the existing
   endpoint-protection cycles.
