# HYP-2203: Moser Unit Spines Separate True-Optimum Traceability from Flip-Gauge Artifacts

**Status:** SUPPORTED by S626 Moser-carrier computation; OPEN for all
non-lattice extremal planar unit-distance configurations.  This supplements
HYP-2201, which verifies that Harborth triangular-lattice optima never flop
through `n=28`, and HYP-2202, which separates graph-level unit HPs from
canonical tiling-order flops in triangular and Moser carrier scouts.

## Claim

For a point set `P`, the intrinsic version of the user's tournament question is:

> Does the unit-distance graph `G_1(P)` contain a Hamiltonian path?

If yes, `P` has a **unit spine**.  Once a unit spine exists, either flip
convention can make the tournament's fixed/tie Hamiltonian path all-unit by
choosing the tie order from that spine:

- if unit edges remain unflipped, the spine order itself is all-unit;
- if unit edges are flipped, the reverse spine order is all-unit.

Thus the flip convention is not the invariant.  The invariant is traceability
of the unit graph plus retention of the spine/tie-order side channel.

## Evidence

`04-computation/unit_distance_unit_spine_tournament_s626.py` builds on the S622
Moser-carrier beam and stores output in
`05-knowledge/results/unit_distance_unit_spine_tournament_s626.out`.

The width-`1200` Moser beam recovers:

- exact values through `n=14`;
- the exact `n=21` value `57`;
- the known `n=22` Moser lane with `60` unit edges.

Every checked witness has a unit spine:

```text
n=2..14: unit spine True for every exact beam witness
n=21:    57 edges, unit spine True
n=22:    60 edges, unit spine True
```

The exact small witnesses also have many unit Hamiltonian paths, not just one:

```text
n=10: 672
n=12: 11976
n=14: 55160
```

This supports the working conjecture that dense Moser-carrier extremizers are
recursively traceable by unit edges.

## The First Apparent Flop Is a Gauge Artifact

Under the lexicographic coordinate order, the point-flip tournament already
stops seeing an all-unit directed Hamiltonian path at `n=7` in the S626 beam:

```text
n=7: unit graph traceable, but lex unit-flip max unit arcs = 5 of 6
```

The same row still has `60` undirected unit Hamiltonian paths.  Therefore the
"Hamiltonian path flops to non-unit pairs" phenomenon can happen immediately
inside a fixed external point order while the underlying geometry still has a
unit spine.  This is not a geometric transition; it is a gauge/tie-order
transition.

## Recursive Pattern

The `n=14` and `n=22` beam spines share an explicit Moser-coordinate ladder.
The `n=14` path is a two-layer motif in the first coordinate:

```text
(1,1,1,-1) -> ... -> (1,0,1,-1)
(0,1,1,-1) -> ... -> (0,0,0,0)
```

The `n=22` path is the same motif with an added outer layer:

```text
(2,1,1,-1) -> ... -> (2,0,1,-1)
(1,1,1,-1) -> ... -> (1,0,1,-1)
(0,1,1,-1) -> ... -> (0,0,0,0)
```

S628 promotes this motif to a proved subtheorem.  THM-408 defines two explicit
rank-4 Moser slab/cap families:

```text
P_m^+ : |P_m^+| = 8m+6,  E(P_m^+) = 27m+6 for m>=1
P_m^- : |P_m^-| = 8m+5,  E(P_m^-) = 27m+3 for m>=1
```

Both families have Hamiltonian paths made entirely of Moser unit-shell steps.
The S626 `n=14`, `n=21`, and `n=22` witnesses are exactly
`P_1^+`, `P_2^-`, and `P_2^+`.  Thus the ladder is no longer just a beam
pattern for these rows; it is an infinite traceable-section theorem.  The
remaining open part is whether every proof-relevant dense Moser extremizer
admits a compatible slab/ear decomposition, or whether a true extremizer can
force incompatible ears.

## What a Genuine Flop Would Require

A real intrinsic flop requires an extremal unit-distance graph with no
Hamiltonian path.  Graph-theoretically this would need a forced traceability
obstruction such as:

- more than two essential leaves or ears;
- a cut vertex whose removal separates three or more path-obligatory branches;
- a dense core plus frontier attachments that cannot be linearly threaded.

S626 sees none of these in the exact small Moser witnesses or the `n=22`
`60`-edge lane.  The likely proof route is an ear/deletion theorem:

> Every dense Moser-carrier extremizer has a removable boundary vertex whose
> deletion leaves a dense extremizer with a unit spine, and the removed vertex
> can be attached to one end of a spine.

The `n=22` frontier version should be phrased over `21`-core extensions:
classify which degree-`4/5` extension vertices preserve a unit spine, then
combine with HYP-2176/HYP-2188's `56/57`-edge core ledger.

## Tournament Analysis

S626 deliberately compares vertex choices:

- **points:** useful diagnostic, but order-sensitive and too lossy;
- **unit-direction pairs:** keeps the Moser shell and direction support;
- **frontier ears / dense cores:** keeps the recursive construction route;
- **traceability obligations:** preserves the actual unit-spine predicate;
- **proof-obligation routes:** keeps obstruction and deletion side channels.

The route tournament is transitive:

```text
unit graph traceability
> frontier-gain recursion
> direction-pair quotient
> spine-order flip gauge
> proof-obligation tournament
> lexicographic point flip
> raw pair-distance tournament
```

It has score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, no directed
3-cycles, singleton SCCs, and one Hamiltonian path.

The challenged assumption is that the tournament's vertices should be only
points and its tie path should be an arbitrary coordinate order.  For the unit
spine problem, the point tournament is downstream of the traceability witness,
not the primary object.

## Next Problems

1. Count unit spines in the five exact `n=21`, `57`-edge graph6 cores from
   HYP-2176, not only the Moser beam representative.
2. Add a spine-preservation score to the S617 `21`-core extension ledger.
3. Run an exact deletion/ear audit through `n=14`: which boundary deletions
   preserve both optimal edge count and traceability?
4. Test whether direction-dropout witnesses in HYP-2194 stay traceable after
   losing one edge at `n=14`.
5. Promote the spine-ladder motif to a canonical Moser-coordinate normal form,
   or refute it by finding a different exact witness with no slab recursion.
