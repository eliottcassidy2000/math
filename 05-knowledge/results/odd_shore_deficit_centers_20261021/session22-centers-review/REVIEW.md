# Independent Session22 odd-shore audit

Current verdict: PASS after the separately documented repair review in SECOND_REVIEW.md and strengthened degree-47 proof review in THIRD_REVIEW.md. The original consumer defect below is retained as historical audit lineage, not an outstanding defect. The initial review did not cover subsequently added proof section 8.

Reviewed source SHA256: 8c8423841ec49af438e5599d3740a68378e50ef58c12496b0382cfbabfec9d82. Review performed 2026-10-21 UTC. The theorem review takes the inherited q <= Delta+2 and standard strong-tournament pancyclicity as its explicit inputs; it does not reauthenticate historical census/algebra certificates.

## Analytic checks

- The half-supported prefix argument correctly requires n >= 2r. Nonidentity powers of a power-of-two r-cycle retain even order; two disjoint conjugate supports preserve that order. Oddness of a tournament automorphism group rules these out. Reversal/converse mates must be counted in pairs, so t qualifying cuts force only 2(t-1), not 2t. In particular t=4 is compatible with Delta=6, whereas t=5 is excluded.
- Independent permutation calculations at r=2,4,8,16 and n=2r-1,2r,2r+1 demonstrate the genuine boundary: the squared reversal product has order 2r-1 just below half-support, but order r at and above half-support. These controls support the precise hypothesis, not an extension below it.
- Rotation shores are cyclically vertex-transitive and hence regular. A regular odd tournament is strong: a nontrivial SCC ordering contradicts equal degrees. The matching complement argument works for q=2 as well as q>2. Size-one shores cause no application gap because a nonempty proper deleted subset is unavailable there, and the obstruction only uses m>=9.
- An eight-cycle in a proper odd shore has strong complement. A_0-to-B row sums are equal by the same literal rotation, without any exact-period assumption on a separately constructed finite frame. The total 8k cannot be 4m for odd m; hence at least five exits or entries. Applying an ambient anti-automorphism handles the entry case and preserves complement strongness and support sizes.
- The regular connectivity inequality is correct: for surviving SCC blocks X->Y, averaging outdegrees gives y<=s and x<=s, so n<=3s. Thus n>=49 permits deletion of a 16-cycle, giving Delta>=14. Odd n yields the stated cutoff 47.
- Combining q<=8 and m<=7 gives n<=56. No assertion that all periods shown by a frame are exact cut indices is made or needed for finite deletion controls. The general bound uses r_s>d and r_1>d+2 strictly, and its odd-order endpoint 3r_1-1 is correct.
- The recurrence-onset conclusion is correctly labelled conditional on the inherited numerator-degree/algebra input; this audit adds no new verification of that algebra.

## Fresh replay and independent controls

Copied producer executions under normal Python and -O both pass exactly 10,351 checks. Source digest and every nonvolatile receipt field agree. Start/end timestamps and wall durations differ and must not be described as byte-identical runs.

All 16 original prefix receipts pass an independent strict validator checking the actual vertex-label universe, exact cycle rotation prefixes, common suffix, terminal, directed adjacency, distinct words, tournament structure, and half-support. All 16 saved representatives happen to be entry cases. Their converse matrices with the unchanged literal words give 16 independent valid exit cases.

A separate bitset reachability implementation and independent adjacency constructor replayed exactly the documented scoped universes: 136 small frames, 3,888 nonempty proper shore-deletion tests, and 2,048 restricted order-18 frames. No full order-18, all-self-converse, exact-period, or all-isomorphism-class census is claimed.

## Consumer finding

Original audit_receipt accepts a genuine receipt after replacing every occurrence of vertex 17 in its words with -1 (and terminal if applicable). Python negative indexing preserves the adjacency accesses, while the current distinct-count condition does not enforce labels in range(n). The independent validator rejects this malformed receipt. This does not invalidate the original receipts or analytic theorem, but the original consumer is not a hostile-input verifier.

Recommended fix: before any indexing, enforce direction in {'exit','entry'}, all C labels and terminal in range(n), terminal outside C, and each word's labels exactly set(range(n)); check integer labels explicitly if JSON booleans are to be rejected. Reject malformed types/keys/shape consistently rather than leaking an indexing exception. In particular, an arbitrary direction string is currently treated as 'entry'.

Artifacts: audit.py, audit_results.json, independent_controls.py, independent_results.json, normal/receipt.json, optimized/receipt.json. Only session22-centers-review was written. No Git operations or remote writes were performed.
