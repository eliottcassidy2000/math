---
id: HYP-1903
status: OPEN
source: codex-2026-05-31-S431
related:
  - HYP-1895
  - HYP-1900
  - HYP-1853
  - THM-354
  - THM-357
  - THM-359
  - HYP-1836
  - HYP-1842
---

# HYP-1903: LRC has a two-nearest-neighbor tournament lift

## Statement

For a reduced LRC speed set `V={v_1,...,v_k}` and time `t`, do not only record

```text
min_i ||v_i t||.
```

Record the full circular configuration

```text
P_t = {0} union {v_i t mod 1}.
```

This configuration carries a natural incomplete round tournament:

```text
i -> j  iff  0 < clockwise(P_i,P_j) < 1/2.
```

Coincident runners create identified vertices, and antipodal pairs create
missing/tied arcs.  It also carries a two-nearest-neighbor digraph, where each
vertex points to the two closest vertices on the circle, retaining ties at the
second order statistic.

The hypothesis is that a minimal LRC counterexample must induce a leafless,
arithmetic-realizable nearest-neighbor protection core in this lift.  In
particular, the scalar endpoint-protection core from THM-357/THM-359 should be
a projection of a richer circular cut-protection core:

```text
zero-side bracket pair       = circular cut around the stationary runner
two-nearest contacts         = local protector arcs
pairwise distance collisions = quotient identifications
antipodal ties               = incomplete-tournament defects
```

## Evidence

`04-computation/lrc_neighbor_tournament_lift_s431.py` audits known tight and
near-tight examples at exact safe/tight times.

The scalar LRC data collapses several examples that the lift separates.

- `initial k=4` and `sporadic n=5` both become the same regular `5`-gon at
  `t=1/5`; the half-distance tournament has regular score sequence
  `(2,2,2,2,2)`.
- `initial k=5` is the regular `6`-gon at `t=1/6`, but the half-distance
  relation is incomplete because of three antipodal pairs.
- `sporadic n=6` has the same lonely scalar threshold at `t=1/6`, but the lift
  exposes a moving-runner collision: `v3` and `v9` occupy the same point.  Thus
  `pairwise_min_gap=0` even though the stationary runner is lonely.
- The near-tight `k=6` example and the two `n=14` ladder examples have the
  stationary runner bracketed by one nearest runner on each side, but their
  pairwise minimum distance is far below the lonely gap:

```text
near-tight k=6:      pair/th = 0.535714
n14 seven-ladder:    pair/th = 0.169913
n14 exported-debt:   pair/th = 0.053301
```

So tiny LRC gaps are not the same as globally regular packings.  They can be
steep valleys inside crowded pairwise configurations.

All audited half-distance tournaments are strongly connected after removing
only antipodal/collision ties, and the two-nearest graph is connected in every
example.  That is not a proof object yet, but it says the right next invariant
is not connectivity alone.  It should be a private-leaf or cut-protection
invariant inside the two-nearest lift.

## Predictions

1. A useful LRC near-counterexample score should include the vector
   `(lonely_gap, second_zero_gap, pairwise_min_gap, antipodal_count,
   collision_count, two_nearest_private_leaf_count)`.
2. Known tight examples should quotient to small cyclic tournaments with
   multiplicities.  Sporadic tightness is likely residue-polygon equality plus
   quotient identifications, not a fundamentally new scalar shape.
3. Any genuine open-cover counterexample should leave a nonempty
   two-nearest-neighbor protection core after deleting private leaves and
   nonreciprocal zero-bracket contacts.
4. The difference-speed set `V-V` should expose false near-disproofs: if
   pairwise distances collapse at every scalar near-witness, then the set is
   not close to a global circular packing.
5. The tournament good-cut theorem THM-354 should have a circular round
   analogue in which unprotected zero-bracket cuts are the LRC boundary
   witnesses.

## Sources

- `04-computation/lrc_neighbor_tournament_lift_s431.py`
- `05-knowledge/results/lrc_neighbor_tournament_lift_s431.out`
- `07-reflections/lrc-neighbor-tournament-lift-s431.md`
- HYP-1895 and `07-reflections/lrc-distance-tournament-two-neighbor-s22.md`
- HYP-1900
- HYP-1853
- THM-354
- THM-357
- THM-359
