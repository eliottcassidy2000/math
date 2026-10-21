# Independent pentagon anti-fixed audit

- `PROOF.md`: complete five-case lattice proof, reduced denominator, recurrence, exact quasiperiod, asymptotic and reciprocity.
- `check.py`: standard-library exact checker; no coefficient fitting or external dependencies.
- `CHECK.json`: normal Python output.
- `CHECK-OPTIMIZED.json`: byte-identical optimized Python output.

Run:

    python check.py > CHECK.json
    python -O check.py > CHECK-OPTIMIZED.json
    cmp CHECK.json CHECK-OPTIMIZED.json

No primary-packet files were edited and no Git operation was performed.
