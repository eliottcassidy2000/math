# Proper induced unit-distance realizations of countable pseudoforests

**Session140 refinement, 2026-11-02. Status: PROVED elementary.** This is a
separate strengthening of the checkpoint1 existence theorem. No frozen
checkpoint1 file is altered. No literature-priority claim is made.

## Theorem

Every countable locally finite simple pseudoforest has an injective **proper
induced unit-distance realization** in R²: adjacent vertices are exactly the
pairs at Euclidean distance one, and every compact subset of R² contains
only finitely many image vertices. The straight edges may cross.

Here “proper” concerns the discrete vertex set. Because every edge has length
one and the graph is locally finite, every compact region also meets only
finitely many drawn edges. This is not an embedding into a fixed lattice
and does not imply a uniform positive lower bound on all pair distances.

## Proof

Each connected component is a tree or has one finite cycle with rooted trees
attached to its vertices. Use one root for an acyclic component. For a cyclic
component, first draw its cycle as a regular polygon of side one. Such a
polygon has no additional unit chords. A generic rotation makes the cycle
vertices' x-coordinates pairwise distinct: equal x-coordinates for any fixed
pair exclude only finitely many rotation directions. The cyclic seed is
finite, so all these exclusions are finite.

Give each component its own open vertical strip of finite width. A cyclic
strip is chosen wide enough to contain its regular polygon with a positive
margin; an acyclic component can use width one. Place the strips successively
from left to right with gaps greater than one and left endpoints tending to
positive infinity. This is possible even when the widths are unbounded.
Thus distinct components automatically have no collisions or unit distances,
and any bounded region meets only finitely many component strips.

Inside a component strip, give each initial root an open horizontal interval
containing its x-coordinate, with closure inside the strip. In the cyclic
case these root intervals may be chosen pairwise disjoint, using the generic
rotation above. The polygon edges themselves need not stay inside
these smaller intervals; the entire seed stays in its component strip.

Maintain the invariant that an already placed vertex v has an assigned open
interval I_v containing x_v. For its finitely many children, choose pairwise
disjoint nonempty open intervals J_w with closures inside

```
I_v ∩ (x_v−1/4, x_v+1/4).
```

There is room for any finite number of such intervals. For each child w,
choose x_w in J_w and put

```
y_w = y_v + sqrt(1−(x_w−x_v)^2),    I_w=J_w.
```

The available positions form a nonempty open upward arc of the unit circle
about v. Avoid all already placed vertices and all unit circles about old
vertices other than v. At a finite stage these exclusions remove only
finitely many points: distinct unit circles meet in at most two points.
They cannot exhaust the available open arc. Thus the new vertex has exactly
its prescribed unit adjacency to the old set and no collision. In a rooted
tree it has no other prescribed old neighbor; all cycle vertices were placed
together initially.

Process the finitely many vertices at each successive rooted-tree depth.
Local finiteness and the finite seed ensure that every level is finite, so
at every placement stage the previously placed set in that component is
finite. Each vertex is eventually processed. Alternatively, over all
components use a fair finite-stage interleaving; the strip gaps already
remove all cross-component constraints.

Every parent-child edge raises the y-coordinate by strictly more than

```
sqrt(1−1/16) = sqrt(15)/4.
```

If b is the smallest y-coordinate of the finite seed, a vertex at rooted
depth d therefore has `y > b+d*sqrt(15)/4` for d>0. A bounded upper height
allows only finitely many depths; their levels are finite. Thus a compact
set meets only finitely many vertices in each component. It meets only
finitely many component strips, proving global properness.

Finally, each pair's correct edge or nonedge status was fixed when its later
endpoint was placed, or when its finite cycle seed was placed. Therefore
the realization is induced. ∎

## Consequence for the geometry tradeoff

Apply the theorem to the simple loopless underlying standard Collatz graph.
There is an unconditional proper induced unit-distance realization, even
though no bounded-dilation injective placement in a fixed-spacing honeycomb
or triangular lattice exists.

The session's independent Collatz lane proves at least `2^r` predecessors
within graph radius `5r` of any root not divisible by three. A uniformly
δ-separated unit-edge placement would then obey the disk-packing bound

```
2^r ≤ (1+10r/δ)^2,
```

which is impossible for fixed δ>0. Hence the proper realization necessarily
has arbitrarily close **nonadjacent** pairs. Properness shows that these close
pairs can be chosen escaping every compact set: a finite set of image
vertices has positive distance from every other image vertex, since each
vertex is isolated and every bounded neighborhood contains finitely many
vertices. Thus finite-point accumulation is **not necessary**; crowding can
occur only farther and farther out in the plane.

The distinctions are:

- **Possible:** proper, induced, unit-edge, injective realization in R².
- **Not supplied:** noncrossing straight edges, a lattice placement, a canonical
  layout, or an effective algorithm deciding unknown components/cycles.
- **Impossible for Collatz:** adding a uniform positive separation of all
  vertices while retaining bounded edge lengths.

The upward coordinate is a construction-depth device relative to chosen
roots, **not** a global height decreasing along the Collatz map. In a possible
acyclic nonconvergent component, an arbitrary chosen root need not lie in the
forward direction of every edge. No convergence conclusion follows.
