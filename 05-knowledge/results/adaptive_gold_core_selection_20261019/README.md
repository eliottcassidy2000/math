# Adaptive first-core-comparison selection

Read `ADAPTIVE_GOLD_CORE_SELECTION.md`. Two actual seven-thin exteriors of the same ten-vertex core have disjoint balanced-core-pair lists, excluding one universal first core comparison. Their entire two-law segment nevertheless has strict Gold by ordering the same two questions adaptively, with normalized margin greater than 1/25.

The note also proves the general more-balanced-first ordering lemma and the exact weak scalar criterion, with a strictness counterexample. Run `python verify_adaptive_gold_selection.py` or `python -O verify_adaptive_gold_selection.py`; both give the frozen 424-guard receipt. Every ambient count is independently computed by integer ideal dynamic programming.
