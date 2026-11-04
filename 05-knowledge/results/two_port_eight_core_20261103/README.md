# Prime twelve-vertex residuals and the eight-private-point barrier

2026-11-03. **PROVED scoped, with FINITE-EXACT dependencies and independent
audits.** Universal two-port reconstruction remains **OPEN**. This packet
finds neither an exact response collision nor a six-chromatic unit graph.

## What advances

**A hidden twelve-vertex module reversal must have a prime core.**
[C1–C4, algebra/PROOF.md](algebra/PROOF.md) combine the inherited necessary
conditions `S=E mod2` and `C=C^t mod4` with a new exact twin-orientation sign.
For twins `a->b`, contraction to `z`, and an exterior vertex `w`,

    C_bw-C_aw+C_wa-C_wb = C_zw(Q)+C_wz(Q) > 0.

Reversing the same named twin pair changes this sign even if the rest of
the tournament changes arbitrarily. Positive scaling of the full response
cannot hide it either. Congruence inheritance treats modules and quotients
in different directions; a separate directed-triangle identity handles
the ten-vertex submodule case. The conclusion applies to a twelve-vertex
module inside an ambient tournament of any size. It does not classify
arbitrary pairs of tournaments or settle prime twelve-vertex cores.

The [specified scout](algebra/normal.json) tests **311,872 labelled
constructions**, with duplicates permitted: **146** pass endpoint parity,
and **none** also pass cover symmetry modulo four. These are fixed families
and deterministic streams, not a complete order-twelve census. The unique
six-vertex cover-mod-four hostile is retained; its endpoint imbalance
demonstrates why the two necessary conditions cannot be conflated.

**Every proper five-colouring of the fixed Moser spindle extends across
any eight private planar points.** [E1–E6, geometry/PROOF.md](geometry/PROOF.md)
retain all unit contacts and every boundary assignment. The complete
`2^28` labelled graph census yields **38** minimum-degree-three classes
passing necessary metric filters. **Thirty** have arbitrary-three-list
polynomial certificates. The remaining eight cannot be realized injectively:
six require a point with three contacts to a rigid Moser spindle, and two
contain a six-vertex circle fan forcing two distinct vertices to coincide.
An independent deletion-based census reproduces all **724,080** labelled
survivors and their full permutation orbits.

The same classification proves an exact small-order boundary: a planar
unit-distance graph on at most eight distinct vertices is three-choosable
**if and only if** it contains no Moser-spindle subgraph. A graph containing
one has choice number exactly four. These are list-colouring assertions
for arbitrary palettes, not tests restricted to five named colours.

An actual connected fifteen-point control has a cubic private cube, exactly
25 unit edges, and all **42** five-colour boundary partitions survive, with
**17,796** extensions per old assignment. The full four-, five-, and
six-colour [boundary tables](geometry/boundary_tables.json) and exact
coordinates are retained. A blocker on this fixed carrier therefore needs
**at least nine private vertices**; existence at nine remains **OPEN**.

## Full observations, not endpoint marginals

[response/RESULT.md](response/RESULT.md) retains the complete named
`S,E,M,C,U,V,D` tensors for a specified prime twelve-core and its converse.
Each four-index three-path tensor has **20,736** entries, including
singleton middle paths. Literal permutations and cuts agree with an
independent squarefree-support product construction. The complete
**6,241-query** interpolation battery for each core, plus nine high-weight
queries each, is checked against **12,500** direct fourteen-vertex
Hamilton-path DPs. Converse covariance is checked on all tensor entries.

The two saved batteries differ in **6,240** positions. The first actual
same-mask query gives **627,597 versus 1,139,039**. This is an exact hostile
to endpoint-parity sufficiency, not a near collision of full endpoint data;
that query already reduces to a one-port observation by universal-source
padding. No equality of a truncated sidecar is substituted for full
two-port equality.

## Connections, losses, and next decisive objects

The inheritance pass and live concepts are in [BOARD.md](BOARD.md). The
closest preceding packet is [two_port_seven_core_20261103](../two_port_seven_core_20261103/README.md),
including its complete order-ten modular sieve, order-eleven self-converse
census, and recovered square-sum, graceful-tree and Petersen maps. Their
earlier finite scopes are preserved, not silently promoted to all orders.

The common principle is to retain the coordinate that composition uses.
Module contraction preserves a particular congruence, while named twin
contrasts retain a sign that unknown positive scale cannot destroy.
Geometric abstraction preserves edges and lists but loses simultaneous
realizability; rigid-circle incidence restores that missing condition.
These are precise maps with different targets, not an asserted equivalence
between tournaments and unit-distance colouring.

The [incoming-work integration](recovery/INCOMING.md) retains physical
names, joint ports, accepted source words and operation scopes. A finite
return is not reusable closure. The [latest integration](recovery/LATEST_INCOMING.md)
adds full ambient modulo-four descent, an exact 1,753-query converse
battery checked against our complete twelve-core data, and the incoming
physical missing44 return with its all-order isolation. It also records
the exact transitive-chart escape and fixed-support arithmetic boundaries.
[The independent reviews](recovery/ALGEBRA_RESPONSE_REVIEW.md)
and [geometry review](recovery/GEOMETRY_REVIEW.md) state the accepted
quantifiers and computational limits.

The next reconstruction target is a prime, even, nonselfconverse
twelve-vertex module with the necessary modular data, followed by its full
response. The next geometric target is the [nine-private-point
residual](geometry/FRONTIER_NINE.md): if its core contains a Moser spindle,
two adjacent hubs must each have exactly two private-spindle contacts.
For every legal private-spindle colouring, a blocker must leave those hubs
the same singleton list. The exact residual pair factor is
`|A_x||A_y|-|A_x intersect A_y|`. Realizing this simultaneously with the old
boundary is still open.

## Reproduction

From the repository root:

```text
python 05-knowledge/results/two_port_eight_core_20261103/verify.py
python 05-knowledge/results/two_port_eight_core_20261103/verify.py --replay
```

The first command checks the frozen artifact inventory, inherited source
hashes, paired receipts, and explicit scope checks. The second also reruns
every registered engine in normal and optimized Python modes, including
the complete eight-graph censuses and full twelve-core tensors. It never
rewrites its expected results. See [VERIFICATION.md](VERIFICATION.md) for
the recorded close-out checks. No literature-priority claim is made.
