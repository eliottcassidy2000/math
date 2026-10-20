# Publication layout and independent verification

The final reviewer approved two editorial clarifications in the eight-file
scientific packet: “six non-top squarefree divisors” in the prime note and
the standard-path row/column of the 25×25 pair grid in the README. The proofs,
maps, code, and numerical certificates are otherwise unchanged; checksums are
refreshed. Original texts and the exact diff are retained in audit/provenance.
The separate surface addendum is not included.

The historical audit report and approval retain their original preparation
paths as provenance. From this published directory, the portable commands are:

    python3 -B verify.py > /tmp/prime-divisor.json
    cmp /tmp/prime-divisor.json receipt.json
    python3 -B -O verify.py > /tmp/prime-divisor.optimized.json
    cmp /tmp/prime-divisor.optimized.json receipt.json
    python3 -B verify_field.py > /tmp/prime-field.json
    cmp /tmp/prime-field.json field_receipt.json
    python3 -B -O verify_field.py > /tmp/prime-field.optimized.json
    cmp /tmp/prime-field.optimized.json field_receipt.json
    python3 -B audit/independent_check.py > /tmp/prime-independent.json
    cmp /tmp/prime-independent.json audit/independent-normal.json
    python3 -B -O audit/independent_check.py > /tmp/prime-independent.optimized.json
    cmp /tmp/prime-independent.optimized.json audit/independent-optimized.json
    sha256sum -c SHA256SUMS

Only four input-file locators in the independent checker were changed for this
layout; the original script and exact diff are preserved under audit/provenance.
The audit verifies the revised publication_inputs.json inventory; the original
exact_manifest.json remains unchanged as audit provenance. The marked
twenty-five-path table is copied unchanged as a finite audit input.
All six replays match the original receipts byte-for-byte. No producer code or
mathematical claim was modified. The inherited source hashes remain dated
provenance, not a claim that historical scratch paths are repository paths.

The incoming trijection_reset_20261020 packet independently recovers the
prime-shape classification and a planar divisor model. This packet retains
its different oriented tetrahedral fibre correspondence and explicitly
marked determinant/icosahedral graph. No literature-priority claim is made.
