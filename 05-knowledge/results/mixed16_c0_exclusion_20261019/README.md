# Mixed(1,6), c=0 and the two-sided reduction

Read `mixed16_c0_exclusion_20261019.md` and `two_sided_two_exclusive_reduction_20261019.md`. All 96,428 cores are covered, with 64,402 residual cores and 696,230 residual orders independently replayed. No cone obstruction occurred.

`union_receipt.json` is the terminal combined receipt; the root progress/replay files refer to an intentionally empty serial prefix. Run `python -O audit_union.py` for quick coverage/provenance verification, or `python -O replay_full.py` for fresh full mathematical replay using only the standard library. The latter is slower.

The combined corollary leaves at least two exclusive neighbors on both sides of any maximum-span critical pair in a hypothetical seven-thin counterexample. It does not prove the full seven-thin balance conjecture. No novelty claim is made.
