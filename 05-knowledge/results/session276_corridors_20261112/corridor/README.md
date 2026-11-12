# Corridor lane: preserve the escape, retain its receiving defect

See [PROOF.md](PROOF.md) for complete statements, proofs, boundaries and scope.

**Positive results:** an actual41->43->45 framed corridor; an exact all-order
classification of left/right h1 contexts preserving the deletion-changing
43->45 gate; an unbounded family of three-value transported corridors; the
more general same-receiving-mass transport law H'=H+2^tV. A finite zero-count
internal-edge repair restores V for H159, and a three-tile version of the
move carries the repaired same-V pair.

**Exact obstructions:** the41/43 C5 cut island survives every neutral tail
length under one/two apex-incidence flips; its sharp boundary costs are8
(core exit) and14 (tail activation). All preserving one-sided h1 transports
have explicitly lacunary shifts. No common same-mark h1 fusion sandwich
increases counts while restoring the core gate's common receiving mass.

**Additional structure:** transitive-deletion seam-charge fibers are exactly
ternary-addressed reflected-pair cubes. Low-jump widths through order14 are
finite evidence, not an all-order theorem.

Chan–Pak and unbounded continuation of the initialized23 corridor remain OPEN.
The three-value gate is not mistaken for an interval covering all heights;
the three-tile repair is explicitly outside the previous one/two-tile graph.

## Reproduce from repository root

```
python3 05-knowledge/results/session276_corridors_20261112/corridor/check.py
python3 05-knowledge/results/session276_corridors_20261112/corridor/escape_check.py
python3 05-knowledge/results/session276_corridors_20261112/corridor/charge_probe.py
python3 05-knowledge/results/session276_corridors_20261112/corridor/c5probe.py 10
python3 05-knowledge/results/session276_corridors_20261112/corridor/ve_corridor.py
python3 05-knowledge/results/session276_corridors_20261112/corridor/reroot_probe.py
python3 05-knowledge/results/session276_corridors_20261112/corridor/repair_probe.py
```

The first two commands are the focused independent literal-graph controls;
their expected outputs are CHECKS.json and ESCAPE_CHECKS.json. They use no
fusion-matrix evaluator for new scalar identities. The charge and fixed-core
probes have complete finite universes stated in their source. The V-level
search fixes six vertices and root0; it finds480 states with values41,43,45.
The reroot probe checks complete root menus only for six named transported
triples. The repair probe checks one/two edge flips only around three named
inputs. None is used as an all-order nonexistence proof.

`escape_search.py`/`escape_search.json` retain the first failed trial: fixed
powers of one h1 context through deletion order8 did not grow counts while
retaining gap2. The successful construction instead varies a single-emission
context's tail; its exact classification is in PROOF.md§6. This experimental
failure was not promoted to an impossibility theorem.

The independent session review is in the sibling `review/` directory.
