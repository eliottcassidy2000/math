# Two-cycle source comparison

COMPARISON.md separates the new verified digit-two 12875 witness, the corrected marked/unmarked correspondence, source-reported finite cycle counts, and conjectural all-order structure.

Reproduce the narrow checks with:

    python3 verify.py > normal.stdout
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    sha256sum -c SHA256SUMS

The incoming order-nine census retains SINGLE-RUN source status. This packet does not rerun it. The new named order-thirteen witness, its actual strong child, and the weak order-fifteen padding boundary are independently checked.
