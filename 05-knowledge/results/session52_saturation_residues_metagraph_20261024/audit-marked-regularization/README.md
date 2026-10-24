# Marked regularization: independent audit packet

Read `AUDIT.md` for the proof review, independent finite scope, and limits.
The complete frozen candidate is in `producer/`; `audit.py` imports none of
its code. `normal.json` and `optimized.json` contain the independent results.
`compare.py` checks agreement and rejects corrupted result near misses.

Reproduce using Python's standard library:

    python3 verify_manifest.py
    python3 audit.py --output replay-independent.json
    python3 -O audit.py --output replay-independent-optimized.json
    python3 compare.py --output replay-comparison.json
    python3 replay.py --output replay-portability.json

The final command creates a fresh temporary directory containing only
manifest-declared actual files, reruns both independent and producer checks
normally and under -O, and verifies their byte-level results. It also checks
that a corrupted relocated producer receipt is rejected. Nothing is published.

Expected finite scope: 1,099 labelled sources through order five; 10,650
single-arc reversals; 7,887 paths; 1,100 recenterings; all three complete
collision witnesses. Ordinary source/image counts are 1/1,1/1,2/2,4/3,12/9;
converse-merged counts are 1/1,1/1,2/2,3/3,10/9. The unmarked image is not all
regular tournaments and these values are not a new recurrence for A000568.
