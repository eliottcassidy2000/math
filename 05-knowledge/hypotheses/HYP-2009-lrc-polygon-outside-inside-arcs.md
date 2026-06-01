---
id: HYP-2009
status: OPEN
source: codex-2026-06-01-S529
related:
  - HYP-1995
  - HYP-1998
  - HYP-1999
  - HYP-2001
  - HYP-2003
  - HYP-2007
  - HYP-2008
  - THM-374
  - THM-381
---

# HYP-2009: LRC is an outside-clasp gap problem with hidden inside chord debt

## Statement

A tournament can be read as a binary relation on the edges of a simplex, but the
LRC-realizable tournaments come from a stricter polygonal source.  The visible
data is the outside cyclic gap necklace of points on a circle.  The hidden
inside arcs are not independent simplex coordinates; each chord orientation is
the half-turn threshold of a consecutive outside-gap sum.  The dihedral group
acts on the outside necklace, and the inside arcs split into chord-length
channels.

For LRC, the stationary observer is a distinguished clasp vertex.  A source
time is an outside statement: the two clasp-adjacent gaps are at least `1/n`.
The runner-runner tournament inside that polygon is shadow data of the outside
necklace, except on boundary walls where diameter and endpoint ties require a
compactifying Hamiltonian tie path.

At the AP/regular-polygon witness for even `n`, the observer is the missing
clasp of the regular `n`-gon.  The clasp gaps are exactly `1/n`, so the source
is closed-wall only, and the runner interior contains hidden diameter ties.  At
`n=14` there are `78` runner-runner inside chords determined by `14` outside
gaps, split into `D_14` channels
`12,12,12,12,12,12,6`; the last channel is the `6` hidden diameter ties.

## Evidence

`lrc_polygon_outside_inside_s529.py` audits clasp-deleted regular `n`-gons for
`n=4..18`.

Selected rows:

```text
n=6:  hidden=10  outside=6   diameter_ties=2
n=8:  hidden=21  outside=8   diameter_ties=3
n=14: hidden=78  outside=14  diameter_ties=6
n=18: hidden=136 outside=18  diameter_ties=8
```

For `n=14`, the AP witness `t=1/14` has clasp gaps `(1/14,1/14)` and channel
inventory `((1,12),...,(6,12),(7,6))`.  The tie-resolved inside tournament has
score histogram `((5,6),(6,1),(7,6))`, `c3=85`, and `H=2641713`, confirming that
the wall compactification carries substantial hidden interior structure even
though the outside polygon is maximally simple.

Tournament Analysis over `n=4..18` by hidden boundary burden is transitive
(`c3=0`, singleton SCCs, one Hamiltonian path), with even frontiers ranked by
their diameter-tie burden before odd open-body rows.

## Predictions

1. The n=14 proof should separate outside clasp-gap forcing from inside
   diameter/endpoint debt.  A non-AP sweep should be forced to open the clasp
   before it can cycle all hidden boundary debt.
2. The HYP-2003 condition "wall-only => AP" should be restatable as: only the
   regular outside necklace can keep the clasp gap closed while its inside
   chord-channel debts remain dihedrally balanced.
3. Boundary-compactified source classes beyond the round/A000016 body should be
   classified by which chord channels hit half-turn or endpoint walls, with even
   `n` diameter channels forming the first obstruction.

## Sources

- `04-computation/lrc_polygon_outside_inside_s529.py`
- `05-knowledge/results/lrc_polygon_outside_inside_s529.out`
- `07-reflections/lrc-polygon-outside-inside-arcs-s529.md`
