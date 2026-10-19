# Mixed(1,5), c=1: complete exact exclusion

Read `mixed15_c1_exclusion_20261019.md`. All 91,604 valid cores are covered, with 63,000 residual cores and 630,808 residual orders independently replayed. No cone obstruction occurred.

`union_receipt.json` is the terminal combined receipt; root progress/replay files describe an intentionally empty serial prefix. Run `python -O audit_union.py` for quick coverage/provenance verification or `python -O replay_full.py` for fresh full mathematical replay. The latter is slower and needs only the standard library.

No full seven-thin theorem or novelty claim is made.
