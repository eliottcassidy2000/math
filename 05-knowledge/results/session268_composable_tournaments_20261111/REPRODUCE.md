# Reproduce the session

From the repository root, using Python3 and only the standard library:

    python3 -B 05-knowledge/results/session268_composable_tournaments_20261111/reproduce.py

This checks the saved full-census hash, runs all seven new programs in normal
and optimized Python, compares deterministic output with the committed records,
and checks the spatial transfer recurrence through N=20. It works in a temporary
copy, so the committed receipts are not rewritten. Timing is excluded from
output equality. `FINAL_CHECKS.json` records the resulting controls.

The complete order 16 modular census is stored losslessly as
`modular/full_receipt.json.gz.b64`; `modular/receipt.json` is its readable summary.
Recover the complete JSON with:

    python3 -c "import base64,gzip,pathlib; p=pathlib.Path('05-knowledge/results/session268_composable_tournaments_20261111/modular'); print(gzip.decompress(base64.b64decode((p/'full_receipt.json.gz.b64').read_bytes())).decode(),end='')" > /tmp/session268-full-census.json

The full record contains every H histogram through 16 and full response buckets
through 12. The hash in the summary checks its exact original bytes.

Construct a target witness in the declared language:

    python3 05-knowledge/results/session268_composable_tournaments_20261111/spectrum/compile_target.py 495 --output /tmp/witness495.json

Rows are integer bitmasks: bit j of row i means i→j. A negative answer applies
only to this compiler's explicit atom/seam/lift language, not all tournaments.
The recorded minimum for495 is9 vertices in that language; it is not asserted
to be the globally smallest tournament realizer.

No inherited Lean build, repository-wide test suite, or complete-prime catalogue
is part of this replay. The measured modular runtime is environment dependent.
