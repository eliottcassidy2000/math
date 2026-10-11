# Adaptive cycle search through signed boundary fragments

2026-10-11. **PROVED:** exact substitution of the boundary oracle below.
**FINITE-EXACT:** the declared comparisons. No plane six/seven decision or
algorithmic speed improvement is claimed.

This integrates the [adaptive odd-cycle search](plane_odd_cycle_channels_20261011b.md)
with the [signed boundary compiler](plane_six_reframes_parity_20261011b.md).
It changes the oracle used on each candidate, retaining the candidate
labels, original graph, complete search quantifier, and actual witnesses.
The earlier complete palette-relation and signed-parity proofs are inherited.

## Exact source-preserving split and lift

Sort the supplied graph's edges and put alternating edges in fragments
\(F,H\). Each fragment receives the endpoints of its edges; assign all
original isolated vertices to \(F\). Their ordered boundary is precisely
\(B=V(F)\cap V(H)\), in increasing original vertex order. Consequently
\(F\cup H=G\), the interiors are disjoint, and no original vertex is lost.
Restrict the supplied ternary labels to the two vertex lists.

For each candidate, `split_oracle` calls `compile_boundary` on both
fragments and `glue` on the authenticated packets. An internal obstruction
or contradictory seam returns an actual monochromatic odd cycle. Remove
the repeated closing vertex to match the search certificate's convention;
optional chord shortening retains a monochromatic odd cycle in the same
supplied graph. A successful seam produces the full original colouring.
One or two channels simply leave some ternary labels unused; the adapter
checks that its lift introduces no such unused channel. Its declared
domain is exactly one, two, or three channels.

Thus this oracle succeeds exactly when the direct residual-graph oracle
succeeds. Either oracle may return a different valid cycle or colouring.
Different cycles can cause later search candidates to differ; equality of
their full search traces is not asserted.

## Complete adaptive certificate

The scheduler takes the first canonical label assignment satisfying all
retained cycle constraints. A failed candidate adds its actual odd cycle,
which excludes that candidate while preserving every proper colouring.
There are finitely many assignments, so this process terminates on each
finite graph. A positive result is its checked full colour witness.
A negative result passes the original `SearchCertificate` verifier:
**every** canonical assignment must make at least one retained cycle
monochromatic. Checking only the candidates previously visited would not
be sufficient.

The integration does not alter either imported module, assume a graph
is geometrically realizable, or infer any universal plane bound.
The exact coordinates for the pentagon and spindle controls come from
their inherited algebraic graph constructors.

## Reproduction and hostile controls

The universe is all 75 labelled graphs of orders one through four, with
all labels in one, two, and three channels: 6,594 direct-oracle comparisons
and 225 full adaptive-search comparisons (201 positive, 24 negative).
Additional controls are the unit pentagon with one channel, the unit
Moser spindle with two, and abstract \(K_7\) with three. Learned cycle
counts for integrated/direct search are respectively \(1/1,3/4,35/35\).
These finite counts are not a runtime theorem.

Deleting a necessary cycle from either negative special certificate is
rejected. A single triangle clause cannot refute all ternary assignments.
Empty graphs, empty fragments, retained isolates, invalid label types,
unsupported channel counts, and a triangle whose obstruction crosses the
seam are explicitly checked.

```text
python -B -X utf8 04-computation/experiments/plane_cycle_boundary_integration_20261011b.py
python -B -O -X utf8 04-computation/experiments/plane_cycle_boundary_integration_20261011b.py
```

[Script](../../04-computation/experiments/plane_cycle_boundary_integration_20261011b.py);
[saved output](plane_cycle_boundary_integration_20261011b.out).
Normal, optimized, and saved output agree: **13,811 explicit checks**.
LF SHA256: `99148b46bbefa1c5975c2511529310a15758cc54ff99446e625e3a90a38276f5`.
