# Mixed(1,3), c=3: complete exact exclusion

Read `mixed13_c3_exclusion_20261019.md` for the self-contained theorem, full-branch corollary, attribution and scope.

Terminal result: all 82,987 degree-valid cores covered, with 39,575 analytic exclusions and 43,412 exact cone-certified residual cores. Every one of the 401,867 residual projected orders was independently replayed. No cone obstruction occurred.

Run `python -O audit_union.py` for quick deterministic coverage/provenance verification. Run `python -O replay_full.py` for a fresh full mathematical replay using only the standard library; this is substantially slower. The already-recorded independent replay was performed for the complete serial prefix and separately before every new shard record was committed.

`union_receipt.json` is the terminal combined receipt. The root `progress.json` and `independent_replay.json` refer only to the serial prefix, which was intentionally paused before sharding. Their partial/pause labels do not describe the combined terminal result. Each shard's `progress.json` has terminal independently-replayed success, and the union audit checks disjoint complete coverage.

No full seven-thin balance theorem or novelty claim is made. Other two-sided root multiplicities remain separate.
