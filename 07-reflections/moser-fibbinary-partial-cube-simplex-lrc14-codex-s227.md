# Moser/Fibbinary Partial-Cube Simplex Carrier

codex-2026-06-26-S227

User request: work on LRC forum posts and merge Moser-de Bruijn, fibbinary,
partial cubes, simplices, doubled triangular numbers, and old `5,6,7` geometry
motifs.

The old automaton results are still the right base:

```text
Moser-de Bruijn = base-4 digits 0/1 = binary 1s only in even positions.
Fibbinary = binary strings with no adjacent 1s.
Moser closes under x -> 4x.
Fibbinary closes under x -> 2x.
Moser is a subset of fibbinary, but the native transition differs.
```

The new carrier is geometric.  Fibbinary words of fixed length form the
Fibonacci cube, an isometric partial-cube language of independent sets in a
path.  Moser values form an even-coordinate Boolean cube inside that language,
but only if the bit-position phase is retained.  That gives LRC a concrete
sidecar: automaton state, native transition, bit-position phase, hypercube
dimension, Theta-class word, gated-subcube status, and median/interval status.

The doubled triangular numbers

```text
2, 6, 12, 20, 30, 42 = n(n+1) = 2*T_n
```

are the directed edge counts of simplex 1-skeleta and the same `u=n(n+1)`
carrier used in the Faulhaber/odd-moment work.  They should be attached as
`directed_simplex_edge_count` or `doubled_triangular_layer`, not used as a
standalone proof scalar.

The `5,6,7` geometry motif remains valuable only with the S225 correction:
`{3,5}/{3,6}/{3,7}` or `(2,3,5)/(2,3,6)/(2,3,7)` give
spherical/Euclidean/hyperbolic, while tournament sizes `n=5/n=6/n=7` are a
separate boundary/pivot/obstruction axis.

Practical next step: annotate a small HYP-2963 packet sample containing AP, GW,
K33, C27 petals, covering, fibbinary, and Moser controls with the
partial-cube/simplex sidecar and test whether any sequence/cube quotient keeps
boundary/open status and route pure after exact `M`, endpoint owners,
magnitude cocycle, closed-arc topology, and geometry-regime fields are present.
