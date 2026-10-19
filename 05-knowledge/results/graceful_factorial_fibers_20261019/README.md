# Factorial graceful-tree fibers and exact reconfiguration distance

Read [the proof](GRACEFUL_FACTORIAL_FIBERS.md) and [independent algebra review](review-algebra/REVIEW.md). The explicit upper-tail forest has exactly k! one-color components, each a k-cube; its width-two graph has the stated signed-permutation distance and linear diameter. General graceful-tree existence and global one-color exchange remain open.

From this directory, run:

    python3 -B verify_factorial_fibers.py
    python3 -B -O verify_factorial_fibers.py
    python3 -B review-algebra/independent_metric_check.py
    python3 -B -O review-algebra/independent_metric_check.py

Compare the first outputs with factorial_fibers_receipt*.json; both independent outputs must match review-algebra/independent_metric_receipt.json. The author checks 222,652 arithmetic and graph conditions. The independent checker verifies 8,122 abstract graph distances. Exploratory component probes are not part of this proved packet. SOURCE_PROVENANCE.json records reviewed source hashes and presentation-only edits.

The original pre-publication author receipts remain under reviewed-original-receipts/. Current receipts differ only in the self-reported paper hash, regenerated after adding reviewed-status and dependency links. All mathematical values are unchanged.
