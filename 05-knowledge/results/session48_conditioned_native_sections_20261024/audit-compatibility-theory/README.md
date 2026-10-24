# Session48 independent compatibility/theory audit

**PASS.** Read `AUDIT.md` for the independent proofs, exact scope, and hostile controls. `normal.json` and `optimized.json` are identical independent-check receipts. `PRODUCER_REPLAY.json` records unchanged producer replay. `INPUT_PINS.json` pins every frozen input, and `PORTABILITY.json` records replay from an isolated temporary copy.

The package is self-contained and standard-library-only. No Git operations are used.

    python -B independent_check.py --output /tmp/session48-independent.json
    cmp normal.json /tmp/session48-independent.json
    python -O -B independent_check.py --output /tmp/session48-independent-O.json
    cmp optimized.json /tmp/session48-independent-O.json
    python -B replay_producers.py
    python -B verify_manifest.py

Do not promote these scoped results into a common global K-host claim for several prefix-key observers, an ordinary pairwise helper construction, a fixed-host forest sufficiency claim, or an unmarked/isomorphism-class count.
