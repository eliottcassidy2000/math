# Independent resource audit

Verdict: PASS for Theorems 1–2, Corollary 3, and their stated restricted constructor scope.

The exact core filter correctly discards quotients with two or more nontrivial SCCs instead of improperly making their components available as separate generators. Both directions of prime-spectrum equality hold: padding contributes one, while two strong components contribute a product with both factors at least three. Reversal and ordered sums create no new prime values. Positive transitive inflation preserves the SCC factors used in the proof.

The inherited triangle bound applies separately to every module of every retained strong core, so every size is at most floor(log2(N−1)); deduplicating actual eligible cores by isomorphism/reversal is safe. Counting parameter tuples proves sum L^q, not realizability of arbitrary numerical capacity values.

The central-binomial estimate is valid: each p-adic exponent has at most floor(log_p(2m)) unit terms, hence its prime power is at most 2m. Removing primes 2 and 7 loses at most two values. Its asymptotic lower bound has leading constant log 2, sufficient for the stated liminf and vanishing-fraction consequences. No prime number theorem is required. These are bounds only for the specified constructor language, not unrestricted B2 circuits or global tournament attainment.

The Session 175 clone-obstruction proof was directly inspected as the inherited source. The other historical bridge paragraphs are explicitly proposals/qualified inherited summaries; this audit does not separately rederive those earlier circuit theorems.

`independent_audit.py` uses its own arc-set representation, graph-search SCC algorithm, endpoint recurrence, literal permutation checks, and prime sieve. It checks 1,099 base tournaments, 285 independently counted single-module inflations (all cross-checked by permutations), and 1,001 prime-bound cutoffs including 10,000. Normal and optimized JSON receipts agree. The author's complete checker was also rerun normally and with optimization; receipts agree with each other and with the submitted certificate.
