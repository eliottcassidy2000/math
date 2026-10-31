# Tournament counting engine

Read [THEORY.md](THEORY.md) for the exact formulas, proofs, retained-state
recursions, prior-work boundaries, and metagraph interpretation.

The central arithmetic family is

`F_n(q,x) = sum_{odd partitions lambda of n} q^E(lambda) x^length(lambda)/z_lambda`.

It computes A000568 as `F_n(2,1)`, self-converse tournaments as
`F_floor(n/2)(4,2^(n mod 2))`, and converse-merged nodes as their average.
The integer cycle-insertion recursion is independently checked and is a
concrete target for the session's separate FRACTRAN compiler.

The [converse-folding packet](../converse_folding/) supplies the actual
half-size class bijection. A marked anti-involution must be retained when
comparing automorphism groups, not just isomorphism classes.

Run from this folder:

```sh
python check_counts.py > normal.json
python -O check_counts.py > optimized.json
cmp normal.json optimized.json
python verify_hostile.py > hostile.normal.json
python -O verify_hostile.py > hostile.optimized.json
cmp hostile.normal.json hostile.optimized.json
```

Main results: [normal.json](normal.json), [counts.csv](counts.csv).
Sources: [SOURCE_PINS.json](SOURCE_PINS.json) and byte-preserved repository
sources in [sources/](sources/). No external source is needed to rerun tests.

Coverage: formulas through order 25; literal signed pair-action Burnside
through order 8; every labelled tournament through order 5; every pair of
four-vertex tournaments with optional vertex coloring. External OEIS
comparisons have separately recorded finite ranges and are not inputs to
the counting algorithms. No metagraph wiring or scalar-node recurrence is
claimed from these counts alone.
