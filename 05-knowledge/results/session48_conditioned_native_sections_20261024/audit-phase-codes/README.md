# Independent native affine phase-code audit

**PASS in the stated scope, after precision corrections separating raw owner values from canonical phases and distinguishing designed from full-budget edges.** Read `AUDIT.md` for independent proofs, finite universes, and the exact correction. The packet has no Git operations, network dependencies, or third-party Python imports.

Reproduce from any working directory:

    python -B independent_check.py --output /tmp/phase-audit.json
    cmp normal.json /tmp/phase-audit.json
    python -O -B independent_check.py --output /tmp/phase-audit-O.json
    cmp optimized.json /tmp/phase-audit-O.json
    python -B replay_producer.py
    python -B verify_manifest.py
    python -B portable_check.py

`normal.json` and `optimized.json` are deterministic independent results. `PRODUCER_REPLAY.json` is separate producer evidence. `INPUT_PINS.json` freezes the source inputs. Historical CUBES snapshots use byte-preserving `.txt` archive names; original source paths and digests are retained. `PORTABILITY.json` records manifest-only temporary-copy validation. `SUMMARY.json` is the concise integration record. Raw owner equations use A o=b+(8r)A1, separately from the converse-head phase shift2A1. The main scope distinctions are higher-arity versus pairwise guards, original low-weight generators versus reduced rows, full-budget versus designed graphs, and continuously held resources versus full-record non-return.
