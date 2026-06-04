# HYP-2202: Unit-Distance Carrier Scouts Separate Geometric Unit Spines From Canonical Tiling-Order Flops

**Status:** OPEN, supported by S625 carrier scout; complements HYP-2201.

**Claim.** In dense small unit-distance constructions, the meaningful
Hamiltonian-path question has at least three layers:

1. Does the undirected unit-distance graph have a Hamiltonian path?
2. Does the nonunit/complement graph have a Hamiltonian path?
3. Under a fixed points-to-tournament tiling order, does the guaranteed
   directed Hamiltonian path use only unit pairs?

S625 suggests that the geometric unit path persists in the tested optimal or
beam-optimal carrier representatives, while the fixed canonical tiling order
loses an all-unit directed Hamiltonian path at `n=7`.

This is the carrier-scout companion to HYP-2201: HYP-2201 records the stronger
triangular-lattice traceability theorem/lower-bound pattern, while HYP-2202
keeps the points-to-tournament mapping diagnostics and the Moser-carrier check.

## Evidence

S625 adds `04-computation/unit_distance_tournament_hampath_s625.py` and stores
`05-knowledge/results/unit_distance_tournament_hampath_s625.out`.

The script tests two points-to-tournament rules:

- **unit-flip:** start from a transitive base order and flip exactly unit pairs;
- **nonunit-flip:** start from a transitive base order and flip nonunit pairs.

It also records the more geometric **snake-base** variant: choose the base order
from a unit Hamiltonian path whenever one exists.

### Carrier Rows

Tested rows:

- triangular carrier through `n=22`, width `300`;
- Moser carrier through `n=14`, width `260`, matching known exact values through
  `n=14` in this scout.

In every tested row, the undirected unit-distance graph has a Hamiltonian path.
No graph-level unit-path flop was observed.

The nonunit/complement graph first has a Hamiltonian path at `n=6`, fails at
`n=7`, then reappears from `n=8` onward in the tested compact rows. The `n=7`
exception is structural: the compact hexagon has a center joined by unit edges
to all six ring points, so the center is isolated in the nonunit complement.

The fixed canonical unit-flip tiling loses an all-unit directed Hamiltonian path
at `n=7` in both triangular and Moser carrier rows. Thus `n=7` is an
order/tiling flop, not evidence that the unit graph itself has lost a
Hamiltonian path.

## Pattern

The recursive pattern is a boundary snake. Compact triangular-lattice animals
grow by attaching perimeter pieces; a Hamiltonian snake can be extended around
the new boundary. This explains why the unit graph can keep a Hamiltonian path
even when the canonical lexicographic tiling order no longer follows it.

The complement has the opposite recursion: it is sparse in the first complete
unit cases, gains a path at `n=6`, loses it at the `n=7` center-isolated
hexagon, and then becomes dense enough to carry paths again once a new point
breaks the center isolation.

Directed Hamiltonian-path histograms under the canonical unit-flip tournament
do not collapse to all-unit or all-nonunit. Their modes sit in the middle. For
example:

| Carrier | n | Mode unit-steps | Max unit-steps | All-unit directed HPs |
|---------|---|-----------------|----------------|-----------------------|
| triangular | `7` | `4` | `5` | `0` |
| triangular | `10` | `6` | `8` | `0` |
| Moser | `7` | `4` | `5` | `0` |
| Moser | `10` | `5` | `7` | `0` |

So the richer invariant is a unit-step profile distribution, not only a binary
"unit path vs nonunit path" label.

## Tournament Analysis

Vertices are mapping/proof lenses, not points:

- unit-graph Hamiltonian path
- nonunit-complement Hamiltonian path
- canonical unit-flip tiling
- canonical nonunit-flip tiling
- snake-base tiling
- boundary-shell recursion
- distance-profile HP histogram
- raw edge optimum

Pairwise observable: which lens better preserves the question "does the
Hamiltonian path live on unit pairs, nonunit pairs, or an order artifact?"

The geometry gauge and tiling gauge are both transitive with score histogram
`{0:1,...,7:1}`, zero directed `3`-cycles, and one Hamiltonian path. They flip
`14/28` edges. Geometry ranks unit-graph HP and boundary recursion highest;
tiling ranks the distance-profile histogram and canonical unit-flip rule
highest.

## Assumption Challenge

Candidate vertices include points, unit pairs, nonunit pairs, distance classes,
Hamiltonian-path obligations, boundary-shell additions, base-order tiles, and
recursive construction moves. S625 uses mapping/proof lenses because the
preserved predicate is whether the guaranteed Hamiltonian path is geometric or
an artifact of the chosen tiling order.

The quotient preserves unit graph HP, complement HP, directed HP unit-step
histograms, and the canonical Rédei insertion path's unit-step count. It
destroys classification of all optimal point sets, continuous deformation
families, and exact less-than-one versus greater-than-one nonunit distances
outside the tested finite carriers.

## Next Test

Search the known `n=11..14` Moser exact rows with larger width, then test
multiple optimal representatives rather than one canonical beam output. The
next real theorem would prove a boundary-snake extension lemma for the
triangular-lattice lower-bound family and then ask which non-lattice optimal
families preserve it.
