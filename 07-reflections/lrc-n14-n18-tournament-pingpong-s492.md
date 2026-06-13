---
source: codex-2026-06-01-S492
status: supplemental synthesis to HYP-1950
tags:
  - lonely-runner
  - tournament-analysis
  - n14
  - n18
  - pressure-peeling
  - nearest-neighbor
---

# n=14 / n=18 Tournament Ping-Pong: k=2 and Deficit Relief

S490 already ran the broad feedback loop between `n=14` and `n=18`.  This
supplement pushed one specific stuck point: maybe the old S470/S490 pressure
graph is too one-neighbor.  If the real obstruction needs the two nearest
neighbors, then nearest-only deletion relief could falsely look peelable.

The noise source was the nearest-neighbor graph literature.  In a metric set,
the nearest-neighbor digraph is directed and strongly constrained; in the
standard unique-neighbor case, directed cycles are only 2-cycles.  The useful
translation was not to import that theorem directly, but to ask:

```text
what if LRC pressure needs k-nearest or threshold-deficit relief,
not just one-neighbor relief?
```

References used as noise:

- Nearest neighbor graph overview:
  <https://en.wikipedia.org/wiki/Nearest_neighbor_graph>
- Mutual k-nearest-neighbor graph language:
  <https://rdrr.io/cran/cccd/man/nng.html>
- Local transitivity as a tournament constraint:
  <https://users.cecs.anu.edu.au/~bdm/data/digraphs.html>

## What S492 Added

`lrc_n14_n18_tournament_pingpong_s492.py` compares eight rows:

```text
initial n=14
initial n=18
n=14 lpd ladder, scale 7 skip 6
n=18 lpd ladder, scale 9 skip 8
n=14 gate ladder, scale 14 skip 6
n=18 gate ladder, scale 18 skip 8
n=14 single-gate repair: replace 6 by 14*16
n=18 single-gate repair: replace 8 by 18*18
```

At selected exact endpoint/gap times it computes three incomplete tournament
shadows:

```text
k1 relief:      deleting j improves i's nearest distance
k2 relief:      deleting j improves sum of i's two nearest distances
deficit relief: deleting j reduces i's two-neighbor threshold deficit
```

The stratified scan keeps every witness/boundary time and caps large endpoint
families at `640` scanned times, with the result reporting
`sampled_times=scanned/total`.

## Main Result

No selected row produced a strict pressure SCC or a directed 3-cycle under any
of the three pressure lifts.

```text
n=14 selected rows:
  k1      cyclic-or-SCC rows 0/23
  k2      cyclic-or-SCC rows 0/23
  deficit cyclic-or-SCC rows 0/23

n=18 selected rows:
  k1      cyclic-or-SCC rows 0/23
  k2      cyclic-or-SCC rows 0/23
  deficit cyclic-or-SCC rows 0/23
```

This strengthens HYP-1950.  The old nearest-only graph was not merely too weak
on these rows: even the two-neighbor and threshold-deficit variants still see a
peelable pressure DAG.

## n=14 Notes

The two classic hard rows behave as expected:

```text
n=14 scale 7:  gap/th=5/924,  unprotected=84
n=14 scale 14: gap/th=5/1848, unprotected=168
```

The single-gate repair row is a useful opposite stress test:

```text
(1,2,3,4,5,7,8,9,10,11,12,13,224)
gap/th=11/672
unprotected=12
```

It has much lower endpoint exposure than the ladders but a larger visible gap.
Still, it also has no k2 or deficit pressure core.  That supports the split
already seen in HYP-1828:

```text
small visible gap       -> quotient ladder, high endpoint debt
small endpoint exposure -> gate repair, larger visible gap
```

Neither side has yet produced a mobile cyclic dependency.

## n=18 Notes

The mixed `2*3^2` rows follow the same pattern:

```text
n=18 scale 9:  gap/th=1/176, unprotected=176
n=18 scale 18: gap/th=1/352, unprotected=352
```

The single-gate analogue

```text
(1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,324)
gap/th=97/3564
unprotected=24
```

again lowers endpoint exposure relative to the big ladders but does not create
a pressure SCC.  Thus n=18 still looks like a mixed-torsion proof laboratory,
not a disproof machine, unless perturbations can finally create a nontrivial
labelled SCC.

## Methodological Update

The next search should explicitly separate three layers:

```text
endpoint handoff graph  = arithmetic ownership/protection labels
mobile pressure graph   = geometric removability/leafiness
safe-gap mask           = cycle edge-cover obstruction
```

A counterexample-like branch should have to survive all three peelings.  If it
has endpoint witnesses, it is not a cover.  If it has a pressure leaf, remove
or charge that private blocker.  If safe gaps have adjacent pairs, the static
configuration has mobile lonely runners and should be translated back through
difference-speed labels before being treated as an origin obstruction.

The concrete next test is no longer "find a smaller gap."  It is:

```text
find a perturbation with k2_largest_scc > 1 or deficit_largest_scc > 1.
```

If bounded perturbations around `n=18` scale `9` and `18` still fail that test,
then HYP-1950 should be promoted from search heuristic to a candidate peeling
lemma.
