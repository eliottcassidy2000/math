# Verify primorial offsets and supersingular predicates

See [THEOREM.md](THEOREM.md) for definitions, elementary proofs, exact membership through 71, inherited boundaries, and the comparison to conditional hole supports.

- `check.py`: standalone exact verifier; optional SymPy-backed certificate generation.
- `certificates.json`: recursive Lucas certificates and exact compositeness witnesses.
- `results.json`: verified sequences, exact complete sets, class-number/point-count rows, and sharp 30030 offset examples.
- `check.out`, `check_optimized.out`: identical normal and `-O` outputs.
- `test_rejections.py`, `rejections.out`: seven deliberately false certificate variants rejected.
- `SOURCES.json`: public literature and inherited repository source pins.
- `MANIFEST.json`: file hashes, reproduction instructions, scope, and status.

Reproduce with `python3 check.py`. No publication was performed by this worker.
