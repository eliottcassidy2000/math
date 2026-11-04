# Independent tribonacci and finite coverage audit

REVIEW.md states the theorem, finite coverage, corrected bookkeeping, and exact limits. This packet retains one strong construction witness for every odd value from 69 through 702189, not merely a count of the covered set.

To reproduce, from this directory:

    python3 verify.py > normal.stdout
    cp receipt.json receipt.normal.json
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp receipt.json receipt.normal.json
    sha256sum -c SHA256SUMS

Only Python's standard library and installed g++ are required. The script recompiles the independent predecessor coverage engine, regenerates the complete witnesses from the retained physical base codes, validates every witness using a successor endpoint engine, and repeats the actual tower/near-join and readiness controls. It does not download anything or run the incoming scripts.

Input code parts and witness parts are each smaller than 1MB. WITNESSES.json describes the ordered line-to-value correspondence and bit encoding. Base codes use the pinned deletion census's column-major upper-triangle convention, documented in the C++ readers. Both actual arc orientations and the root name are retained in every witness.

The source-defined union uses near-join parameters1–8 and tower heights2–10. Height1 adds five low values but does not increase the verified upper interval. The coverage is finite and does not supply an all-order selector, a joint-port grammar return, or an H-neutral transport with a fixed number of names.
