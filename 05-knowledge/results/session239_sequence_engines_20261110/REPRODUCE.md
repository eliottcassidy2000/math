# Reproduce Session239

Run from the repository root. All correctness scripts use standard Python3.
The historical C/GMP comparison additionally requires a C compiler and GMP
headers/library; that benchmark is not required to use the new Python engine.

Fast final correctness replay:

    python3 05-knowledge/results/session239_sequence_engines_20261110/replay.py
    python3 -O 05-knowledge/results/session239_sequence_engines_20261110/replay.py

Individual checks:

    python3 05-knowledge/results/session239_sequence_engines_20261110/audit/verify.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/audit/strong_gap_audit.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/audit/palette_audit.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/verify.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/prime_gap.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/c7_low_fibres.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/precision/precision.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/verify_cli.py

The audit scripts use explicit runtime checks that remain active under `-O`.
Normal/optimized arithmetic conclusions agree; timing fields are not expected
to be byte-identical. Archived benchmark values are one-machine observations.

Slower performance replay (writes fresh local timing JSON next to each script):

    python3 05-knowledge/results/session239_sequence_engines_20261110/burnside/benchmark.py --orders 40 60 80 --repeats 3
    python3 05-knowledge/results/session239_sequence_engines_20261110/burnside/scale.py
    python3 05-knowledge/results/session239_sequence_engines_20261110/burnside/validate_extended.py

The benchmarks import the inherited response engine from
`05-knowledge/results/tournament_future_boundary_20261020/verify.py`.
The extended historical comparison extracts tracked C source with `git show`.
These inherited repository dependencies are not optional for those comparisons.
They are unnecessary for `generate.py` and the delivered new engines.

`SHA256SUMS` pins all deliverable files except itself and the top-level manifest.
After deliberately rerunning a benchmark or changing proof prose, refresh
hashes rather than interpreting obsolete timing hashes as a failed theorem.
