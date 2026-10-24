# Disjoint compatible charge sections

This packet separates one-demand hole coverage from simultaneous full-record allocation. Its main result is an exact, cancellation-free square-zero occupancy compiler, plus a precise Hall specialization and sufficient depletion tests for whole-resource hyperedges. It includes a genuine native two-demand obstruction and a six-state native example in which the exact 7-clock boundary is not enough to overcome 12-clock internal phase components.

Read `THEOREM.md` for proofs and scope. The native supply keeps three C3 kernel lifts per normalized unit, complete K/converse closures, literal shared-owner relations, class-dependent guards and exact per-slot clock costs. It does not transport the original numerical dyadic current bank, solve arbitrary heterogeneous source/target demands, or manufacture one global observer.

Run:

    python -B check.py
    python -O -B check.py
    python -B verify_manifest.py
    python -B reproduce.py --output /tmp/disjoint-relocation.json

`normal.json` and `optimized.json` are identical mathematical receipts. `PORTABILITY_REPLAY.json` records a real manifest-only copy into a new temporary directory and execution there, not merely a relative-path inspection. The final manifest adds that replay receipt after its successful run, avoiding a circular checksum: the replayed payload hashes and the manifest hash used at that time are explicit in the receipt.

Historical Markdown is copied byte-for-byte to `.txt` files, with original repository paths, commit and SHA256 digests in `SOURCE_PINS.json`. Those historical links are archival text, not live links from this packet. No Git mutation or upstream edit is performed here.
