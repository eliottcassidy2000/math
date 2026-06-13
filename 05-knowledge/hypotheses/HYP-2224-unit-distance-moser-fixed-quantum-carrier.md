# HYP-2224: Unit Distance Has a Moser Fixed-Quantum Carrier

**Status:** OPEN proof-use hypothesis with exact THM-408 carrier scout.

## Claim

The unit-distance analogue of the aliquot/Pillai fixed-point carrier is not a
new scalar equation on `u(n)`.  It is the fixed edge-channel increment of the
THM-408 Moser spine ladder.

For both Moser families `P_m^+` and `P_m^-`, after the cap transient, adding
one full slab has the same direction-pair carrier increment:

```text
Delta_total = (0,1,8,4,0,4,4,2,4)
sum Delta_total = 27.
```

The fixed quantum splits as

```text
Delta_spine = (0,0,0,1,0,1,4,0,2), sum 8
Delta_bulk  = (0,1,8,3,0,3,0,2,2), sum 19
27 = 8 + 19.
```

Thus the unit-distance `27` is a fixed carrier under the add-one-slab
operation:

```text
one slab -> 8 new traceable-section edges + 19 hidden bulk edges.
```

This is the unit-distance face of S646/HYP-2222's `C=27` fixed-clock lesson.
The raw scalar `E` is not enough; a proof route must retain the
spine/bulk/direction side channel.

## Evidence

S648 imports the THM-408 symbolic Moser ladder and computes total, spine, and
bulk direction counts for both families through `m=5`.  In both `P_m^+` and
`P_m^-`, every transition from `m>=1` to `m+1` has the same stable carrier:

```text
Delta_total = (0,1,8,4,0,4,4,2,4)
Delta_spine = (0,0,0,1,0,1,4,0,2)
Delta_bulk  = (0,1,8,3,0,3,0,2,2)
```

There are `9` antipodal direction pairs in the rank-4 Moser unit shell, so

```text
27 = 3 * 9.
```

But the quantum is not uniformly distributed across directions:

```text
Delta_total - (3,3,3,3,3,3,3,3,3)
= (-3,-2,5,1,-3,1,1,-1,1).
```

The scalar `27` hides a large direction defect.

## Frontier Reading

For the exact `n=21` Moser row `P_2^-`:

```text
vertices=21
edges=57
spine_edges=20
bulk_edges=37
pure_bulk_edges=20
```

For the `n=22` Moser lane `P_2^+`:

```text
vertices=22
edges=60
spine_edges=21
bulk_edges=39
pure_bulk_edges=21
```

So at the current `n=21/22` frontier, pure-bulk direction mass equals the unit
spine length.  This equality is not a general `m` identity; it singles out the
`m=2` frontier rows as especially balanced section/bulk carriers.

The `n=22` row also has

```text
degree_hist={3:1,4:3,5:8,6:5,7:5}.
```

By the S614 deletion reduction, a hypothetical `61`-edge `22`-point row cannot
have degree `<=3`, because deleting such a vertex would leave more than
`u(21)=57` edges.  Therefore any `61`-edge improvement inside this carrier must
repair the degree-`3` cap endpoint channel, or else leave the fixed Moser
carrier entirely.

## Proposed No-Leak Target

The proof-facing version is:

```text
fixed 27 quantum + unrepaired cap endpoint -> at most 60 edges
61 edges -> cap endpoint repaired or different carrier
```

Equivalently, the next search should not ask only whether a graph has `60` or
`61` edges.  It should retain:

- direction-pair counts;
- spine vs bulk labels;
- cap endpoint degree;
- deletion-core edge counts;
- endpoint-compatible ear labels;
- totally-unfaithful obstruction labels.

This is the unit-distance counterpart of the S646 LRC no-leak target:

```text
mass-changing row -> loose
mass-fixed row -> AP/Vstar or killed by side labels.
```

## Tournament Analysis

S648 uses proof-carrier lenses as vertices.  Alternate vertices considered:
points, unit edges, direction pairs, slabs, caps, ears, deletion cores,
obstruction labels, and proof obligations.  The chosen quotient preserves
section/bulk proof data and destroys raw point order.

The majority tournament is transitive with one Hamiltonian path, ranking

```text
fixed_27_edge_quantum
> spine_bulk_direction_split
> cap_endpoint_repair_channel
> pure_bulk_direction_jackknife
> degree_deletion_core_ledger
> traceable_section_word
> triangular_lattice_baseline
> raw_edge_count_only
```

Fingerprints: `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`,
`directed_3_cycles=0`, singleton SCCs, and `hamiltonian_paths=1`.

**See also:** HYP-2222, HYP-2221, HYP-2217, HYP-2204, HYP-2203,
HYP-2188, HYP-2176, THM-408, `04-computation/unit_distance_moser_fixed_quantum_s648.py`,
`05-knowledge/results/unit_distance_moser_fixed_quantum_s648.out`,
`07-reflections/unit-distance-moser-fixed-quantum-s648.md`.
