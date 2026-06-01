---
id: HYP-1932
status: OPEN
source: codex-2026-06-01-S454
related:
  - THM-372
  - HYP-1931
  - HYP-1930
  - HYP-1904
  - HYP-1895
  - HYP-1900
  - THM-002
---

# HYP-1932: Tournament Analysis switchboards are the primary relation layer

## Statement

Tournament Analysis should treat the pairwise comparator family itself as a
mathematical object.

A **switchboard** is:

```text
one bit-channel c_ij(x) for every unordered pair {i,j},
a fixed Hamiltonian-path tie-break for zero channels,
and the resulting path T(x) in {0,1}^{binom(n,2)}.
```

Rank lifts assign one scalar to each object.  Switchboards assign one channel
to each relation.  If the problem lives in pairwise tension, the switchboard is
closer to the object than any scalar rank shadow.

S500 canonized the construction as THM-372: a switch value for every unordered
pair, together with a fixed Hamiltonian tie path, determines a unique
tournament.  Finite wall sets therefore produce piecewise-constant tournament
movies.

## Evidence

`tournament_analysis_switchboard_s454.py` extends S23 by measuring live edges,
flip counts, Hamming diameter, Hamiltonian-path ranges, and SCC ranges across
basketball, circle-runner, cuboid, sphere, simplex, and LRC switchboards.

The rank/analyzer split sharpened:

```text
rank shadows:     transitive states 172/172
analyzer shadows: transitive states 31/672
analyzer mean live edges per path: 22.95
```

The basketball example shows why the comparator is the creative layer.  The
same synthetic five-starter data supports a pass-flux tournament, an
assist-rank tournament, a reciprocity switch, a two-hop lens, and a pressure
switch.  The rank lift is always transitive, while flux and switch lifts move
through cyclic states with Hamiltonian-path ranges such as `1..13` and `3..15`.

The continuous examples show that arbitrary pair metrics can be made meaningful
by declaring the switch:

```text
circle chord-annulus switch:   H=33..571
cube Linf-resonance switch:    H=43..609
sphere greatcircle switch:     H=51..651
simplex KL-skew lens:          H=9..585
```

For LRC snapshots, the safe-distance switch separates initial skeletons from
structured ladders:

```text
n14 initial: safe-switch cycles 0
n14 seven-ladder: safe-switch cycles 93
n16 initial: safe-switch cycles 0
n16 dyadic ladder: safe-switch cycles 139
```

Thus the marked lonely-runner bracket is only one visible constraint.  The
surrounding pairwise moat switchboard carries additional tournament structure.

## Predictions

1. Every Tournament Analysis artifact should declare its switchboard type:
   flux, rank, symmetric switch, lens, phase, volume, or hybrid.
2. Feature extractors should record `live_edges`, wall count, Hamming diameter,
   Hamiltonian-path range, SCC range, and exact-path clusters.
3. Two metrics that induce the same tournament path should be identified as the
   same lens for proof purposes, even if their geometry looks different.
4. LRC proof search should keep both the marked two-neighbor lonely bracket and
   at least one pairwise switchboard shadow, because scalar brackets can miss
   cyclic crowding.
5. Human systems such as basketball lineups are legitimate low-dimensional
   laboratories for Tournament Analysis: the base labelled Hamiltonian path is
   a playbook/tie-break layer, not bookkeeping.

## Sources

- `04-computation/tournament_analysis_switchboard_s454.py`
- `05-knowledge/results/tournament_analysis_switchboard_s454.out`
- `07-reflections/tournament-analysis-switchboard-s454.md`
- `04-computation/tournament_analysis_metric_lifts_s23.py`
- `07-reflections/tournament-analysis-metric-lifts-s23.md`
- THM-372
- HYP-1931
