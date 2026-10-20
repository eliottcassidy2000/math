# One-swap mechanical Collatz exclusion

Read [the proof](ONE_SWAP_EXCLUSION.md). The specified one-swap mechanical family contains no positive integer cycle except the ordinary (q,p)=(2,1) case. General Collatz is not resolved.

The proof uses the [published sharp norm theorem](../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md), [critical square-root-seven gap](../collatz_one_swap_critical_gap_20261020/CRITICAL_CLOCK_SQRT7_GAP.md), [Matveev input](PRIMARY_SOURCES.md), exact rational continued fractions and a complete small-word certificate. The [repetition appendix](REPETITION_FACTOR_OBSTRUCTION.md) removes coprimality, with a [separate audit](audit/NONCOPRIME_REPETITION_AUDIT.md). [Independent review](audit/DIOPHANTINE_CLOSURE_AUDIT.md).

Run from this directory:

    python3 certify_one_swap_closure.py > /tmp/closure.replay.json
    python3 -O certify_one_swap_closure.py > /tmp/closure.optimized.replay.json
    cmp closure.json /tmp/closure.replay.json
    cmp closure.json /tmp/closure.optimized.replay.json
    python3 audit/verify_diophantine_independent.py --closure-report closure.json > /tmp/closure.independent.json
    cmp audit/diophantine.independent.json /tmp/closure.independent.json
    python3 verify_repeated_word_obstruction.py > /tmp/repetition.replay.json
    cmp repetition-receipt.json /tmp/repetition.replay.json
    python3 audit/verify_nonprimitive_repetition.py > /tmp/nonprimitive.replay.json
    cmp audit/nonprimitive.receipt.json /tmp/nonprimitive.replay.json
    sha256sum -c SHA256SUMS

All checking code is Python standard-library-only. Published Collatz convergence scans, remembered irrationality constants, and exploratory obsolete cutoffs are not dependencies.
