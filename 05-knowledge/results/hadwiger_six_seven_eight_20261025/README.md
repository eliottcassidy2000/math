# Hadwiger–Nelson through six, seven, and eight

2026-10-25. **PROVED scoped elementary results / CITED classical and
inherited results / FINITE-EXACT controls.** The investigation and its
limits are in the [synthesis](SYNTHESIS.md).

This continues [six/seven geometry](../six_seven_geometry_20261025/README.md),
[Gale and branch cubes](../six_seven_gale_cube_20261025/README.md),
[Session61's HN connections](../session61_hadwiger_six_seven_20261025/SYNTHESIS.md),
and the concurrent [Session62 bridge packet](../session62_agent_theme_bridges_20261025/README.md).
Those results are recovered and credited, not rediscovered under new IDs.

| Lane | Concrete outcome | Proof / independent audit |
|---|---|---|
| Metric fold | Every planar unit cube with exactly seven distinct images is the A2 flower; exact loss of translations and adjacency reflection | [Result](metric_spindle/RESULT.md) / [audit](audits/METRIC_AUDIT.md) |
| Projective eight | Classical equianharmonic S(3,4,8) transported to the marked conic labels; chirality, basis and affine-origin hostiles | [Result](projective_eight/RESULT.md) / [audit](audits/PROJECTIVE_AUDIT.md) |
| Six paired colors | Eight endpoint frames, six facets, seven affine directions; C2 x S4 residual and common-fixed-point obstruction | [Result](fano_eight/RESULT.md) / [audit](audits/PALETTE_AUDIT.md) |
| Boundary relation | Complete35-relation pair-shadow classification; one missing cubic count; two exact prescribed six-palette unit graphs | [Result](cube_boundary/RESULT.md) / [audit](audits/BOUNDARY_AUDIT.md) |

The same marked Walsh basis now separates cube coordinates, divisor
relations, pair interactions and completion parity. It does not identify
their target predicates. In particular the prescribed boundary examples
are not unprescribed six-color obstructions. The plane six/seven decision,
integer-square hourglass existence and LRC(14) are not closed here.

## Reproduce

From repository root, with Python3 and the standard library:

```sh
python3 -B 05-knowledge/results/hadwiger_six_seven_eight_20261025/verify_packet.py
python3 -B 05-knowledge/results/hadwiger_six_seven_eight_20261025/replay.py
```

The first checks the complete file set and exact bytes. The second runs
four producer and four independent audit scripts, each normally and with
Python optimization, then compares their JSON data with frozen receipts.
Raw normal/optimized outputs must match one another byte for byte. JSON
comparison to saved receipts permits platform-specific line endings only.
The optional figure generator uses matplotlib; it is not a mathematical
verification dependency. Its SVG and PNG were visually inspected.

Universes are stated per lane. The largest include all256 Boolean relations,
all65,536 Boolean-table gluing pairs, all6,561 ternary-weight tables, all8!
permutations of the affine eight, all336 projective maps, all28 collision
pairs, and221,184 independently audited transport systems. These are finite
controls supporting proofs with explicitly stated broader quantifiers.

`MANIFEST.json` freezes the packet. `SOURCE_PINS.json` identifies inherited
sources and the inspected Git states. [Publication checks](PUBLICATION_CHECKS.md)
record verification and concurrency handling. [Concept board](CONCEPT_BOARD.md)
records the inheritance pass, Anchor/Niche/Wildcard portfolio, decisive
tests and final comparisons. No Lean build or axiom audit is claimed.
