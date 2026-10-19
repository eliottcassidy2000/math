# Frozen c=4 saturation/frontier packet

This packet proves the saturation zero-increment lemma and verifies one limited four-comparison screen obstruction together with its repair. It does not establish the full c=4 case, all seven-thin posets, or the general 1/3–2/3 conjecture.

Read mixed12-c4-frontier-stage/seven_thin_saturation_increment_lemma_20261019.md and mixed12-c4-frontier-stage/c4_saturation_and_migrating_screen_frontier_20261019.md. Both rigorous audits use the standard library:

    python -O mixed12-c4-frontier-stage/verify_saturation_increment_lemma_20261019.py
    python -O mixed12-c4-frontier-stage/verify_c4_restricted_screen_obstruction_20261019.py

The exploratory sample scripts and original receipts are retained separately from the proofs. Their directory-relative dependency is included; NumPy and SciPy are needed only for exploratory LP searches. The exact final audits above do not import them or that dependency.

This is a fixed snapshot taken before the complete c=4 scan. No claim that the sampled cores or orders exhaust c=4 is made. Manifest hashes exclude manifest.json itself. The original c≤3 frozen packet is unchanged.
