---
id: HYP-1930
status: OPEN
source: codex-2026-05-31-S470
related:
  - THM-357
  - THM-365
  - HYP-1853
  - HYP-1900
  - HYP-1910
  - HYP-1920
  - HYP-1921
---

# HYP-1930: LRC pairwise moat pressure is an incomplete tournament peel certificate

## Statement

For a speed set `V` with `n-1` nonzero speeds, include the stationary runner
and study the time-slice positions of speeds `{0} union V`.  The pairwise
distance data carries two tournament-like structures:

```text
semicircle orientation:
  i -> j if j lies in the clockwise open half-circle from i;
  collisions and antipodal pairs are missing arcs.

blocker-pressure orientation:
  i -> j if deleting i improves j's nearest-neighbor moat more than
  deleting j improves i's moat.
```

The second construction uses the first two nearest-neighbor distances with
multiplicity.  The conjectural proof use is:

```text
Every all-protected LRC endpoint core either has a pressure leaf, giving a
private endpoint/runner peel, or it has a strict pressure SCC whose labelled
arcs satisfy an arithmetic endpoint-cycle constraint of the THM-365 type.
```

Thus nearest-neighbor data alone is not enough.  The second-nearest runner is
the datum that measures deletion relief.

## Evidence

`lrc_pairwise_tournament_s470.py` sampled exact endpoint and gap times for
canonical rows:

```text
initial n=14
n14 seven-ladder
n14 S380 gate ladder
initial n=15
initial n=16
```

The strict pressure graph had no directed 3-cycles and no nontrivial SCC in
every selected critical row.  The known near-disproofs therefore look
pressure-peelable in this mobile pairwise sense.

At the tight initial rows, the two-nearest correction is essential:

```text
initial n=14, t=1/14: origin d1,d2 = 1,1 thresholds
initial n=15, t=1/15: origin d1,d2 = 1,1 thresholds
initial n=16, t=1/16: origin d1,d2 = 1,1 thresholds
```

The first nearest neighbor alone would miss that each witness is a two-sided
tight moat.  In the semicircle tournament, the even rows expose the expected
antipodal missing arcs:

```text
n=14: 7 antipodal ties
n=15: 0 antipodal ties
n=16: 8 antipodal ties
```

For the n=14 near-disproof rows, the origin gap times are not close to a
pressure core:

```text
seven-ladder gap midpoint:
  t = 3053/25872
  origin d1,d2 = 29/24, 27/22 thresholds
  pressure strict/tie arcs = 7/84
  strict pressure triangles = 0
  largest pressure SCC = 1

S380 gate-ladder gap midpoint:
  t = 4339/51744
  origin d1,d2 = 4339/3696, 29/24 thresholds
  pressure strict/tie arcs = 6/85
  strict pressure triangles = 0
  largest pressure SCC = 1
```

The S380 row also shows a useful warning: at `t=1/2`, one moving runner is
very lonely relative to the other positions, but the origin collides with the
gate-heavy pile.  Mobile pairwise loneliness and the reduced origin LRC
condition are related by translated difference sets, not identical for a
fixed `V`.

## Predictions

1. Add pressure-SCC and pressure-leaf features to the endpoint/IP search
   ledger.  Near-counterexamples with empty pressure SCCs should be easier to
   peel than rows with cyclic pressure cores.
2. A genuine open-cover counterexample should force a labelled pressure cycle
   after all private endpoint leaves are removed.
3. Even denominators should keep producing incomplete circular tournaments at
   their unit skeletons, with antipodal ties as completion debt.  Odd
   denominators should more often produce complete circular tournament
   skeletons.
4. The `n=14` gate fan tax from HYP-1910/HYP-1920 should be combined with a
   pressure-leaf rule: if the forced odd fan plus even bridge has no pressure
   SCC, the bridge branch should expose a private endpoint or product-depth
   row.

## Sources

- `04-computation/lrc_pairwise_tournament_s470.py`
- `05-knowledge/results/lrc_pairwise_tournament_s470.out`
- `07-reflections/lrc-pairwise-distance-tournaments-s470.md`
