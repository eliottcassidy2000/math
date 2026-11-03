# Independent extra audit

Run from this directory with the sibling `../machines/check.py` present:

    python3 audit_extra.py > normal.json
    python3 -O audit_extra.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

`audit_extra.py` dynamically imports the production `check.py` and records its SHA-256 in each receipt. It checks 3,198 single-head input traces (N=1..6, every Boolean vector, targets 0..4N+4), including each intermediate head position, prefix state, translated evolution and independently calculated final output. It checks translation consistency for all 6,400 pairs of local cells in the complete N=2 tape-symbol/head-state alphabet, including multihead collisions. Finally, it compares both PRIMEGAME streams for N=6..10 with an independent proper-divisor primality test. It does not independently simulate the entire FRACTRAN lists or prove the general phase law. SHA256SUMS covers the script and matching receipts, which were generated under normal and optimized Python.
