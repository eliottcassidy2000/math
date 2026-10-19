# Mixed(1,4), c=2: complete exact exclusion

Read `mixed14_c2_exclusion_20261019.md`. All 87,137 valid cores are covered, including 54,910 residual cores and 519,617 residual orders independently replayed. No cone obstruction occurred.

`union_receipt.json` is the terminal combined receipt. The root progress/replay files refer to an intentionally empty serial prefix; their pause/partial labels do not describe the completed four-shard union.

Run `python -O audit_union.py` for a quick coverage/provenance check, or `python -O replay_full.py` for a fresh complete mathematical replay using only the standard library. The latter is slower. All mathematical checks were already independently performed before committing the shard records.

No full seven-thin theorem or novelty claim is made.
