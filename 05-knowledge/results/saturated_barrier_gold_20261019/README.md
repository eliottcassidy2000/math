# Saturated barriers and two free ports

Read [SATURATED_BARRIER_GOLD.md](SATURATED_BARRIER_GOLD.md).

A connected dominating set of saturated incomparability vertices forces a two-sided exterior decomposition from induced containment alone. Applying this to the previously known twelve-point no-good-pair seed proves that every seven-thin poset containing it induced has one fixed strict-Gold tree with all margins greater than E/9.

A separate affine-gap/deletion lemma removes the degree bounds at two attachment ports. A 144-generator exact cone certificate gives the sharp fixed-tree margin infimum 386/6245 for that unbounded attachment family.

Run `python verify_saturated_barrier_gold.py` and `python -O verify_saturated_barrier_gold.py`; the standard-library receipts agree. General seven-thin balance and Gold remain open. No historical novelty is claimed.

## Independent reproduction

Read the [independent review](review/REVIEW.md). From this directory, run `python -B review/independent_controls.py` and repeat with `python -B -O review/independent_controls.py`. Compare stdout with the respective `review/independent-stdout.json` and `review/independent-optimized-stdout.json` receipts. The checker writes its full result and count-table JSON files beside itself. It reads the included `review/original-normal.json` author receipt for the two full-table hash comparisons and imports no author code.

`SOURCE_PROVENANCE.json` preserves reviewed hashes and documents the status and presentation updates. Both constants are sharp for the stated fixed tree in their respective families.
