# Session145 attachment carrier extension

The main proof is `POOLED_PARITY_CARRIER.md`. It strengthens the inherited
endpoint-block carrier by pooling internal vertices across path components
while retaining every component's parity. The new closed coordinate counts
start `1,1,2,7,28,126,571,...`. The order-six minimum is exactly **571**,
by a matching independently audited exterior-coefficient rank certificate,
improving the old upper bound 661. Higher-order sharpness remains open.
It also classifies exactly which signed forest monomials have identical
exterior-coefficient response columns.

## Reproduce the discovery checks

    python -B check_carrier.py --output normal.json --reps reps6.json
    python -B -O check_carrier.py --output optimized.json
    cmp normal.json optimized.json

The checker is self-contained and uses only Python's standard library. It
contains no removable asserts and imports no inherited or auditor code.
It independently checks:

- the full Hamilton polynomial coefficient formula by direct multiplication
  over every vertex permutation at orders 1 through 6;
- all graph edge sets through order 6, both carrier counts, and the first
  strict old-to-new quotient;
- constructive zero/nonzero separating exterior contexts for every pair of
  distinct signatures through order 6;
- all exterior monomials on two new vertices through old order 4;
- the exact strictly positive six-vertex indistinguishable laws and their
  distinct old-coordinate expectations;
- cycle, star, odd-path, zero-pair, and optimized-execution hostile controls.

These discovery tests complement, rather than replace, the independent
audit's pullback checks and literal Boolean tournament counts.
The source files are frozen under `sources/`; `SOURCE_PINS.json` records
commit, exact repository paths, and byte hashes. The archived source bytes use `.md.txt` filenames so their original relative
Markdown links are not presented as links at a different location. Original
repository paths and exact hashes remain recorded in the source pins.

## Independent exact minimum certificate

    python -B independent/verify_exterior_rank.py
    python -B -O independent/verify_exterior_rank.py

The independent subdirectory contains its frozen context bank, receipts,
actual-graph checks, closure tests, hostile controls, and proof review.
The proof explains precisely why exterior-coefficient rows are rational
linear combinations of actual fixed-cut response queries.

Discovery normal and optimized runs agree byte-for-byte: **598,066 active
checks**, including 33,867 full-Hamilton coefficients, 171,010 signature-pair
separators, 18,376 exterior congruence comparisons, and all 32,768 positive
law weights. No all-order response-rank theorem is inferred from these checks.
