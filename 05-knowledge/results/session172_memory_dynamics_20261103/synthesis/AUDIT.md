# Independent audit of S172 synthesis

**PASS.** The author applied and I verified two minor precision edits: the source set X is explicitly nonempty, and the entropy distinction is H(class)<=log2 q<=ceil(log2 q), with equality throughout only for uniform classes and power-of-two q. Neither changes any theorem or application.

## Closure proposition

The deterministic lower bound and attainment are valid for records supporting the declared finite-horizon responses, without an autonomous-update claim. The positive-state lower bound correctly uses disjoint nonempty supports and extremality of deterministic outputs; positivity and exact pure-source response assumptions are essential. Refinement is monotone. Equality of two successive relations makes each generator preserve the current relation: append that generator to each word of length <=k, then use equality at k+1. Induction yields all-future equivalence. Finite equal class counts imply equal partitions because a proper refinement of a finite partition increases its class count. Bounded nondecreasing finite class counts must have equal successive terms; unbounded counts preclude any finite exact autonomous record. The horizon-index difference from the observer sibling is explicitly and correctly stated.

## Sibling cross-checks

Read exactly [sibling observers](../observers/NOTE.md), [sibling collatz](../collatz/NOTE.md), [sibling carry](../carry/NOTE.md), [sibling spectrum](../spectrum/NOTE.md), and [sibling squares](../squares/NOTE.md).

- The affine full-branch observer hypotheses include unit slopes and branch integrality. The synthesis states these, preserves the b-adic domain, and correctly separates same-horizon failure from lower-precision update.
- The 5k-vertex tape, all-word parity bijection, carry recurrence, fixed-weight offset faithfulness and hidden conditional entropy are accurately summarized. The exact descent inequality and absence of a bounded horizon remain arithmetic claims separate from observation. The 5x+1 cycle and the 3x+1 rational point 19/101 agree with the source.
- The strong wrapper summary preserves M>=2, the exact H formula, full-profile identity, canonical articulation, code-distance/measurement-error distinction, and absence of an arithmetic-novelty claim. Its unique articulation language is correct because each initial side has at least five vertices.
- The carry summary accurately states triangular trace blindness and external lower-left probe recovery for the marked fixed-diagonal family. J is not represented as a legal Collatz operation. The guards and exact carries match the source.
- The square summary uses oriented counts throughout, the correct rank n(n+1)/2 for n>=3, fixed labelled read-head 0 versus 4, actual two-context 0 versus 48, and inherited arithmetic joining coefficients 3 and 10. It does not transfer the all-graph rank to the restricted arithmetic family.

## Agenda and links

The proposed cylinder tree expressly excludes 1 from strict descent obligations and requires genuine well-founded control rather than finite density. Two-feedback-vertex storage is explicitly conjectural. Arithmetic graph interfaces and legal carry probes are new proof obligations, not consequences claimed from shared terminology. Resource categories are kept distinct.

All seven Markdown relative sibling-link occurrences have valid target notes under the intended published sibling directory layout (`observers`, `collatz`, `carry`, `spectrum`, `squares`). They do not resolve directly from the current flat `/tmp/s172_*` work-directory arrangement; publication must preserve the stated sibling layout. No external URL appears in the synthesis itself.

## Independent finite controls

`check.py` exhausts all 5,898 deterministic systems on one through three states with two total transition maps and binary observations. It checks monotone class counts, equal-count/equal-relation equivalence, congruence, permanent stabilization and finite-state stabilization, plus the declared link mapping: 93,421 checks. Normal and Python -O receipts agree. These controls support transcription; the all-source and infinite-source claims rest on the proof above.

No forbidden path was read or scanned. No repository mutation was made.
