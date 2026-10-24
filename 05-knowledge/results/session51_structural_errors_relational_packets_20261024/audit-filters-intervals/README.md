# Verify the independent prime/interval audit

Read `REPORT.md` for separate conclusions and `SOURCE-REVIEW.md` for cited-source scope.

Reproduce with standard-library Python:

    python3 verify_inputs.py
    python3 verify_manifest.py
    python3 independent_check.py
    python3 -O independent_check.py
    python3 hostile_certificate_checks.py
    python3 portable_audit.py
    python3 copies/prime-offsets/check.py
    python3 -O copies/prime-offsets/check.py
    python3 copies/implication-intervals/verify_manifest.py
    python3 copies/implication-intervals/check.py
    python3 -O copies/implication-intervals/check.py
    python3 copies/implication-intervals/portable_check.py

`verify_inputs.py` verifies all copied payloads and locally archived source pins without an external checkout. To additionally reproduce the source-to-copy and exact Git-blob comparison used in `input-integrity.json`, run it with `--source-workspace /path/to/original/session/workspace`. That optional check expects the original sibling producer directories and session50-publication/clean-repo. The copied producer checks, independent mathematics, and hostile mutations do not require that checkout. Normal and optimized independent output is byte-identical. The audit's manifest hashes all delivered files except itself and Python bytecode caches.
