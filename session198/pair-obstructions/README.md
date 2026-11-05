# Source realizability barriers

This packet proves architecture-specific obstructions, not noncoverage of the full C1/C2/C3 spectrum.

- `PROOF.md`: exact dyadic one-path source slice; nonsemilinearity even for primitive sources; limits of finite congruence tests; a false formal C2 value67 and its actual C3 repair; persistent source-content divisors and necessary primitive inputs for prime targets; zero-density output from a finite bank with only sink padding.
- `verify.py`: literal gadget enumeration, independent bounded source inventory, source graph checks, exhaustive small C2 exclusion, and direct Hamilton-word count of the C3 repair.
- `VALIDATION.json`: normal/optimized replay result. The checks use explicit exceptions, so optimization does not remove them.
- `audit/`: independent review and an independently coded hostile-control verifier.

Reproduce from this directory:

    python3 verify.py
    python3 -O verify.py
    python3 audit/verify_guardrail.py
    python3 -O audit/verify_guardrail.py

The verifier uses only the Python standard library. `core.py`, `source_oracle.py`, and `c3.py` are unchanged Session197 helpers used as comparison paths. New source constants are independently obtained by counting literal gadget words and covers.

Frozen structural inputs:

- `dependencies/BILINEAR_C2.md`: Session197 `session197/bilinear-c2-inverse/PROOF.md`.
- `dependencies/C3_TWO_SOURCE.md`: Session197 `session197/c3-two-source/PROOF.md`.
- `dependencies/ROOTED.md`: the ordinary rooted source grammar and gadget compiler.

The source-gcd example intentionally records divisibility preservation, not equality of gcd: gcd can increase. The barriers concern a presentation or restricted source language. Replacing a trapped presentation by unrelated sources remains allowed in the unrestricted coverage question.
