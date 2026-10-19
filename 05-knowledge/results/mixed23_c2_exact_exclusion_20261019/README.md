# Mixed (2,3), c=2: exact exclusion

Read `mixed23_c2_exact_exclusion_20261019.md`. All 87,137 valid cores are covered, including 57,686 residual cores and 476,837 residual orders independently replayed. The proof explicitly uses Peczarski's known six-thin theorem for 8,450 short-span orders. No cone obstruction occurred.

`union_receipt.json` is the terminal combined receipt. Root progress/replay pause labels describe an empty serial prefix, not the completed union.

Run `python -O audit_union.py` for quick coverage/provenance verification, or `python -O replay_full.py` for fresh complete standard-library mathematical replay. Every record was already independently replayed before commit.

The combined c≤2 statement links the separate c=0,1 result. The c=3 case and the full seven-thin conjecture remain outside scope. No novelty claim is made.
