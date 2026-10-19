# Independent review: the five-tree two-sided Gold cover

19 October 2026. Verdict: PASS. The theorem is justified for the specified ten-vertex full root core and arbitrary finite seven-thin exteriors. No claim about other cores or general seven-thin Gold follows.

## Reviewed result

For the induced core with predecessor masks

    (0,0,0,6,0,19,2,18,86,151)

and roots (0,3), one of the five displayed two-comparison trees has every path margin strictly greater than E/10000, where E is the ambient extension count. The same second comparison is used on both branches of each candidate tree. The choice of tree may depend on the ambient core law.

Reviewed author sources:

- TWO_SIDED_GOLD_ADAPTIVE_COVER.md
- verify_two_sided_gold_cover.py
- verification.json
- The predecessor BOUNDARY_CHAIN_COMPRESSION.md and ONE_SIDED_ARBITRARY_EXTERIOR_GOLD.md

The numerical discovery file cover_certificate.json was not used as a proof certificate.

## Independent finite reconstruction

The accompanying independent_gap_audit.py does not import author code, profile lists, numerical discovery arrays, or event counts.

1. It constructs profiles from nested order ideals (lower side) and order filters (upper side), enforcing root zeros, support cardinality, and each coordinate's degree slack. This independently yields exactly 46 lower profiles and 48 upper profiles. Coordinatewise shared budgets leave 1,764 compatible pairs.
2. For every pair it constructs the complete ambient relation as a set of ordered pairs and independently takes transitive closure. Closure adds no edges. All 1,764 models are acyclic, seven-thin, have precisely the required full root core, and retain critical-root predecessor/successor inclusions. Their orders range from 10 to 17.
3. It enumerates the 11,971 core linear extensions, grouping them into 3,546 boundary-word classes, and records all five joint answer-cell vectors in each class.
4. It obtains interleaving multiplicities from an integer transfer table for nonincreasing gap words. This is a different counting method from the author's direct ideal dynamic programming on each enlarged atom. All arithmetic is integral, with an explicit bound excluding signed-64-bit overflow.
5. It reconstructs every atom's total and every joint cell. The resulting complete table has SHA-256

       d4caae3e8354652af27e52ed8c33e213ff12e462168933719abeeb44c6f65547

   exactly matching the author's independent direct-DP table.
6. It verifies all ten universally safe path inequalities, all eight positive-integer certificate inequalities, and coverage of all 24 nonautomatic failure patterns. Every stated exact minimum agrees with the author's receipt. The smallest normalized certificate bound is 3205/27793233, strictly larger than 1/10000.

Normal and optimized independent runs produce identical receipts. The audit uses NumPy only for bounded integer matrix multiplication; the author's standalone verifier remains standard-library-only.

## Infinite-exterior argument

The full-core hypothesis supplies the required exterior decomposition without an additional unverified assumption. An exterior vertex is comparable to both incomparable roots, so it lies either below both or above both. A lower exterior vertex cannot be above a core vertex, and the upper statement is dual. The two exterior sides are mutually ordered through either root. Critical orientation is inherited because each outside vertex is simultaneously below both roots or above both roots.

For a fixed lower exterior linear extension, let alpha(q) be the number of exterior letters after q's last exterior predecessor. Once that exterior order is fixed, q's permissible insertion gaps are exactly the terminal alpha(q)+1 gaps. If h=max alpha, the initial exterior segment of length |A|-h must precede every core letter and can be removed uniquely. A chain of length h therefore preserves each individual core-word lift count. The upper side is dual. Since A<B and the roots separate the two interleaving ranges, the two multiplicities multiply.

Summing these exact per-core-word identities over exterior linear-extension pairs gives the claimed compression identity. Normalizing gives each occurrence weight e(C_alpha,beta)/e(P), not the uniform exterior-order weight. Repeated profiles can be aggregated. Every joint core comparison event is preserved by the same identity, so all eighteen homogeneous linear certificate inequalities transfer simultaneously to the actual ambient law.

The complete profile universe is justified: actual terminal gaps are bounded by actual exterior incomparability counts, are monotone in the required directions, vanish at both roots, and have support contained in the incomparable neighborhood of the final lower or first upper exterior vertex. Thus support size is at most seven, and the lower and upper profiles share each core vertex's slack budget. Conversely, the enumerated profiles give genuine thin chain models, as the independent relation checks confirm.

## Cover logic and strictness

The five-step contradiction in Section 4 is valid even if a tree has several bad paths. Assuming every tree has at least one path with margin at most epsilon E, the safe list and I1 force T1 to fail on FT or FF; I2-I3 force T2's bad path to be TT; I4-I5 force T3's bad path to be TF; I6 forces T1's bad path to be FT; and I7-I8 exclude both remaining T5 bad paths. This contradicts failure of every tree.

Each certificate is a positive combination of the hypothesized bad margins, with strict lower bound epsilon E times its coefficient sum. It therefore excludes simultaneous badness exactly as stated. Transfer occurs before elimination, so the proof does not assume that individual Gold atoms automatically imply a Gold mixture. It proves the stronger cover of the entire atom hull. The claimed strict margin E/10000 follows without a limit or boundary gap.

## Editorial observation

In Section 5, 'six distinct first-comparison numerators' should read 'six distinct comparison numerators'; only four of the six comparisons occur as first questions. This does not affect the proof or code. The author applied this correction before freezing checkpoint1.

No blocking mathematical or implementation issue was found. No author, repository, checkout, or Git files were edited.
