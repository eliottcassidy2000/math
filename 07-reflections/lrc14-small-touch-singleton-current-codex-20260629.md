# LRC14 Small-Touch Singleton-Current Reflection

This is a companion refinement to the HYP-3478 small-touch geometry atlas, not
a competing replacement for it.  The atlas proves the finite singleton-pocket
shape; this note records the best-touch current split inside that shape.

HYP-3478 is a useful cleanup after the exception-frontier router.

The six small-touch rows all have zero-edge dead-cover projections, so the
pair-current failure is structural: there is no edge to remove.  The audit
therefore changes the target from "find a better gate tuple" to "classify the
isolated singleton payload."

The split is small:

```text
clean unit-delta: random_covering_062, random_covering_086, random_covering_101
delta sidecar:    random_covering_039, random_covering_074
asymmetric touch: random_covering_001
```

The clean rows are the best first theorem attempt: best gate delta `(1,1)`,
best current `(1,1)`, and minimum E/branch mirror-orbit deltas
`((1,1),(1,1))`.  The delta-sidecar rows are close but must retain their
`(2,1)/(1,2)` mirror sidecar.  `random_covering_001` is the asymmetric row:
four isolated dead components and best current `(2,1)`.

Assumption challenged: "small-touch" was too coarse.  The proof should carry
three singleton packets, not a single catch-all exception label.
