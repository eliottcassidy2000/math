# Scope clarifications verified before publication

## Tangled arcs versus pre-existing closed components

The original Section C stated that endpoint matching determined the closure component count without explicitly repeating Section B's arc-only premise. Under the standard convention permitting closed tangle components, this wording was too broad: one arc and the same arc plus a disjoint circle have identical endpoint matching but different closure component counts.

The final publication revision now explicitly assumes exactly r arc components and no pre-existing closed components for the matching-only formula. It also states the general extension m+c(P,Q), retaining the number m of closed components. This repairs the issue without changing the matching recurrence, contextual quotient, braid example, or Paley results.

The amended independent discovery checker adds actual three-edge circles to matching graphs for r=1,2,3 and m=0,1,2,3. All 76 equality checks and 57 nonzero-circle hostiles passed under normal and optimized Python. The publication audit's separate binary-incidence-rank checker extends this to r≤4, giving 496 cases. The original snapshot is preserved; the issue is resolved in the final scope.

The m term refers to the input tangle's closed components. The allowed connector Q remains a matching realized by arc connectors, not a more general tangle carrying additional closed circles.

## Original fence intersections

The original fence proof used the inherited model in which two original fences meet in at most one point. Its checker already rejected overlapping segments. The publication revision now explicitly excludes collinear overlap in its opening definitions.

This is an accurate clarification of the inherited THM-4505/4509 model. It makes the minimum-degree and original-segment accounting premises explicit. No numerical result, graph class, theorem, or geometry checker changed.

`SOURCE_EVIDENCE.json` records original and final proof hashes, plus the amended tangle-checker hash. These are reviewed revisions, not edits to the preserved discovery snapshots.
