# Moser Unit Spines and Flip-Gauge Tournaments (S626)

The user's proposed rule is natural: start with a transitive tiling of point
pairs, flip a tile when the distance is exactly `1`, and ask whether the
tournament's mandatory Hamiltonian path is made of unit pairs or non-unit
pairs.  But the computation quickly separates two questions that should not be
merged.

Incoming HYP-2201 gives the lattice traceability theorem, and HYP-2202 gives
the carrier-scout reduction: the base Hamiltonian path can be all-unit exactly
when the unit-distance graph is traceable, while canonical tiling orders can
flop earlier as a gauge artifact.  This S626 packet is the complementary
Moser-carrier check: it probes the non-lattice `n=22` lane where a genuine
true-optimum flop could first live, and it isolates the lexicographic
point-order flop as gauge noise.

The intrinsic question is graph-theoretic:

```text
Does the unit-distance graph contain a spanning path?
```

Call such a path a **unit spine**.  If a unit spine exists, then either flip
convention can be made to carry an all-unit mandatory path by choosing the tie
order from that spine.  Unit-unflipped gives the spine order; unit-flipped gives
the reverse spine order.  So the flip convention is not the invariant.  The
invariant is traceability plus the retained tie-order witness.

S626 tested the Moser-carrier witnesses produced by the S622/S617 beam.  The
width-`1200` run recovered exact values through `n=14`, exact `n=21` with `57`
edges, and the known `n=22` Moser lane with `60` edges.  Every one of these
witnesses has a unit spine.

The exact small witnesses are not barely traceable.  They have many unit
Hamiltonian paths:

```text
n=10: 672
n=12: 11976
n=14: 55160
```

This makes the answer to "where does the Hamiltonian path flop?" two-layered:

- **Intrinsic flop:** not observed through the exact small Moser witnesses or
  the retained `n=22` `60`-edge lane.
- **Gauge flop:** already visible at `n=7` if the point order is fixed
  lexicographically.  The unit graph is traceable and has `60` unit Hamiltonian
  paths, but the lexicographic flip tournament's best directed Hamiltonian path
  uses only `5` of `6` unit arcs.

That is the clean warning.  A point-order tournament can claim the path has
"flopped" while the geometry still has a perfectly good all-unit spine.  The
flop belongs to the chosen gauge, not the point set.

The recursive pattern is striking.  The `n=14` spine is a two-layer
Moser-coordinate ladder, and the `n=22` spine is the same ladder with one extra
outer layer:

```text
n=14:
(1,1,1,-1) ... (1,0,1,-1)
(0,1,1,-1) ... (0,0,0,0)

n=22:
(2,1,1,-1) ... (2,0,1,-1)
(1,1,1,-1) ... (1,0,1,-1)
(0,1,1,-1) ... (0,0,0,0)
```

This suggests a spine-ladder recursion for dense Moser clusters.  The right
proof statement may not be "all optimal point sets have a point-order
tournament whose Redei path is unit."  A better statement is:

```text
Dense Moser-carrier extremizers admit recursive boundary-ear deletion
that preserves a unit spine.
```

For `n=22`, this should be attached to HYP-2176/HYP-2188: any putative `61`
edge set extends a `56/57`-edge `21`-core by a degree-`5/4` vertex.  The next
spine question is whether all such proof-relevant core extensions preserve
traceability, or whether the first genuine flop would have to be a three-ear
extension obstruction.

Tournament Analysis should therefore use multiple vertex sets:

- points for diagnostics;
- unit-direction pairs for carrier support;
- dense cores and frontier ears for recursion;
- traceability obligations for the actual predicate;
- obstruction filters for the proof route.

The S626 route tournament ranks these as:

```text
unit graph traceability
> frontier-gain recursion
> direction-pair quotient
> spine-order flip gauge
> proof-obligation tournament
> lexicographic point flip
> raw pair-distance tournament
```

The result is transitive, with no directed 3-cycles and one Hamiltonian path.
That is almost too neat, but it says the right thing: make the unit spine the
base object, and let tournaments organize the side channels rather than replace
them.
