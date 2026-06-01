---
id: HYP-1942
status: OPEN
source: codex-2026-06-01-S481
related:
  - HYP-1840
  - HYP-1866
  - HYP-1891
  - HYP-1905
  - HYP-1910
  - HYP-1930
  - HYP-1940
  - HYP-1941
---

# HYP-1942: n=14 and n=18 first-even LRC rows are bridge-fiber export problems

## Statement

The first-even LRC rows `14 = 2*7` and `18 = 2*9` should be attacked by the
same bridge-fiber proof lemma:

```text
unit/half fan invoice -> one forced bridge fiber -> exported endpoint frontier.
```

At a local `n`-gate, lower columns must cover the `2n` endpoints owned by the
`n`-speed.  In both `n=14` and `n=18`, the exact lower-cover size is `8`:
six unit residue columns, the half-gate, and one bridge.  The difference is
the bridge fiber:

```text
n=14: bridge in {2,4,6,8,10,12}
n=18: bridge in {6,12}
```

Thus `n=18` is not locally larger at the gate invoice.  It is locally more
rigid: the square `3`-torsion row collapses the bridge fiber to the two
residues `+-6 mod 18`.  Its difficulty is not the local gate cover but the
endpoint debt exported into `3^2` and `3^3` layers.

## Evidence

`lrc_n14_n18_tournament_pingpong_s481.py` alternates exact checks between the
two rows.

Local gate cover:

```text
n=14 forced=(1,3,5,7,9,11,13), exact=8, covers=6
n=18 forced=(1,5,7,9,11,13,17), exact=8, covers=2
```

Row-parent/gate/double-gate ladders preserve exact gap-debt products:

```text
n=14 gaps=(5/924, 5/1848, 5/3696)
     debts=(84, 168, 336)
     products=(5/11, 5/11, 5/11)

n=18 gaps=(1/176, 1/352, 1/704)
     debts=(176, 352, 704)
     products=(1, 1, 1)
```

The frontier mass is likewise stable by row:

```text
n=14 frontier_mass = 66/7
n=18 frontier_mass = 316/27
```

Tournament Analysis at selected times separates rank-like gauges from
edge-local gauges.  Mixed two-neighbor slack collapses to transitive orders,
while close-threshold and deletion-pressure switches retain cyclic/SCC
structure.  For example, at the `n=18` gate midpoint the close-threshold gauge
has `60` cyclic triples and SCC sizes `17,1`, while pressure has `22` cyclic
triples and SCC sizes `16,1,1`.

## Predictions

1. A general first-even lemma should prove: after the unit/half fan invoice,
   every bridge fiber either reopens a small-denominator row or exports
   positive endpoint frontier mass.
2. `n=18` should be easier locally than expected because its bridge fiber has
   only two options, but harder globally because the exported rows occupy
   several square-torsion depths: `{2:+r,3:+2}` and `{2:+r,3:+3}`.
3. Rank-style two-neighbor summaries are diagnostic controls, not proof
   objects.  The proof object should be an edge-local deletion-pressure
   tournament or a close-threshold switch tournament.
4. A useful dual certificate should weight bridge-fiber rows and endpoint
   frontier rows simultaneously, rather than charging speeds one by one.

## Sources

- `04-computation/lrc_n14_n18_tournament_pingpong_s481.py`
- `05-knowledge/results/lrc_n14_n18_tournament_pingpong_s481.out`
- `07-reflections/lrc-n14-n18-tournament-pingpong-s481.md`
- HYP-1840
- HYP-1866
- HYP-1891
- HYP-1905
- HYP-1910
- HYP-1930
- HYP-1940
- HYP-1941
