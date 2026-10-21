# Third independent review: terminal-score refinement

Date: 2026-10-21 UTC. Verdict: PASS. This separately reviews the strengthened section 8 giving degree at most 47; neither the initial review nor SECOND_REVIEW.md is retrospectively claimed to cover this refinement.

The setting is explicitly strong self-converse T with Delta<=6, exact marked-cycle cut index q>=5, and odd shore size m>=5. The repaired checker is unchanged at SHA256 44ae88a6d25a922348e668ec537d70f567bdb1b21120635c567fbd920ec4a992. Its finite scoped controls are not a census supporting this new analytic claim.

For a score class of one shore, the previously verified four-prefix argument provides four black orbits outside the marked-cycle cuts. For a score class of c=2,3,4 shores, choose a four-cycle S in one shore and an external terminal from a different shore of the same full score. Rotation makes the A-to-B row sums constant; 4k=2m is impossible for odd m. Thus a terminal with at least three exits or entries exists.

In the exit case every constructed word has the four cycle vertices followed immediately by that terminal, giving five consecutive vertices of the same full score. Because c<q, the marked cycle's cyclic shore-score sequence has no run of that score longer than c. This remains true across the boundary between successive q-blocks and for every linear cycle cut. Therefore none of these selected orbits is a cut orbit. Their reversal/converse mates have five consecutive vertices of the complementary score. Self-converseness preserves the score-class cardinality, so the complementary score also occupies exactly c shores and cannot supply a cut word either.

The entry dual is valid without any hidden score mismatch: work in the opposite tournament and reverse the marked Hamilton cycle. All full degrees become n-1-s, including both S and the chosen terminal, so their equal-score relation is preserved. The exact cut index, automorphism group, black-orbit count, and shore-score multiplicities are unchanged. Hence the same argument supplies four non-cut black orbits there and therefore in the original tournament.

Adding these four to the inherited cut-family lower bound q-gcd(2,q) gives at least eight black orbits. Thus every score class must occupy at least five shores. With q<=8 there can only be one score class, so the tournament is regular. Even q=6,8 gives even order and rules out regular tournaments; q=5,7 remains possible. The regular connectivity/sixteen-cycle argument excludes every regular tournament of order at least 49 in this deficit range, not merely q=1, so q=7,m=7 is ruled out.

The resulting maximal orders by q are 47,14,21,28,35,18,35,24 for q=1,...,8, respectively (q=1 bound is odd-order). Hence degree R6<=47. Under the separately inherited numerator-degree bound degree R6+12 and denominator recurrence, coefficients satisfy the safe recurrence from n>=60. This does not establish minimal degree, an earliest recurrence onset, or any new coefficient census.

As a supplemental hostile boundary control, score_run_control.py checks all binary cyclic shore-score words with q=5,...,8 and c=1,...,4, repeated five times, at every cut. Their maximum run never exceeds c. This finite word control is explicitly not an enumeration of tournaments. Actual total is recorded in score_run_results.json.

Only the review directory was written. No Git or remote writes were performed.
