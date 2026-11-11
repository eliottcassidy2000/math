# Reproduction and verification scope

From the repository root, run:

```sh
python3 05-knowledge/results/session272_yang_mills_bridges_20261111/reproduce.py
python3 05-knowledge/results/session272_yang_mills_bridges_20261111/reproduce.py --with-uv
python3 05-knowledge/results/session272_yang_mills_bridges_20261111/package_checks.py
```

The first invocation uses Python's standard library and exact rational/integer
controls. The optional slower UV diagnostic requires mpmath (tested version
1.3.0); it compares 35- and 50-decimal quadrature and is not a directed-interval
certificate. The analytic U1 proof, not the numerical trend, proves the limit.

Exact finite universes: 100 two-clock controls, a two-point reflection-positivity
hostile, two full-support three-point laws, one Gaussian product identity;
4,096 disk and 256 torus link assignments with all 16 face subsets each;
4,088 word/phase combinations and 21 parity-transfer powers. Floating heat
values printed after the exact spectral/loop controls are diagnostics only.

The source construction and inherited repository theorems were trusted and
read, not rebuilt. The new UV argument and spectral statements received
independent analytic reviews, with precise repairs recorded. The remaining
bridge lemmas have proofs in the notes and targeted finite controls, not Lean
formalizations or exhaustive verification of arbitrary hypotheses.

The initial repository-wide `agents/check_docs.py` found a pre-existing
startup-size failure. Concurrent upstream commit 0e63417e46a2 repaired it; after
integrating that work, the full documentation check PASSED (179987 bytes).
Scoped Python syntax, links, JSON, source pins and new code controls also pass.
No general mathematical proof audit or CI success is claimed from these tests.
