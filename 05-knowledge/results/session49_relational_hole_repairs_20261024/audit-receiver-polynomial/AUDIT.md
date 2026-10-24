# Final-revision review: primitive completion polynomial

This supplement preserves the already-frozen independent receiver-threshold audit. It reviews the producer's final mathematical revision, pins its final manifest, and verifies the additional completion-polynomial consequence. It does not replace or silently rewrite the earlier audit's pre-supplement source snapshots.

## Mathematical finding

Let N_(s,c) count the normalized interval-indecomposable trades of exact slot/clock cost(s,c), and put t=(c−2s)/2. In this receiver t is exactly the changed-triple indicator. Each stripped trade has equal endpoint grade s/2, changes s−3t of the16 independent base coordinates, and leaves16−s+3t base coordinates untouched. If t=0, the common triple may independently be UUU or VVV.

The exact primitive-pair grade generating function is therefore

    sum_(s,c) N_(s,c) z^(s/2) (1+z)^(16−s+3t) (1+z^3)^(1−t).

This counts unordered endpoint pairs once. Removing common V coordinates gives their unique normalized trade; conversely, common context choices preserve grade/charge equality and give a bijection of accepted endpoint intervals after subtracting the common constants. Primitivity is consequently preserved. No charge projection is used to invent a compatible assignment: the frozen native receiver has exactly these independently selectable base records and its one equality-linked triple.

The polynomial coefficients over grades0..19 are

    0,0,2,72,862,5245,18827,43929,72367,90616,
    90616,72367,43929,18827,5245,862,72,2,0,0.

Their sum is463,840 primitive endpoint pairs, in exact agreement with the independent native census.

## Why all worst strict endpoints are counted

Every nonzero primitive has support at least4. A same-grade endpoint pair in a19-slot binary receiver has even support, at most18. A nonprimitive pair recursively decomposes into at least two nonzero primitive events whose changed-slot supports are disjoint. Hence every constituent has support at most18−4=14, and the actual receiver cost formula bounds its clocks by30. Thus every nonprimitive pair admits a strict14-slot/30-clock route.

The primitive catalogue's only cost34 type is(16,34), with196 normalized trades. Each changes the triple plus13 base bits, has stripped grade8, and leaves exactly three common base choices. Therefore the complete strict34-clock endpoint polynomial is precisely

    196 z^8(1+z)^3.

Its coefficients196,588,588,196 at grades8,9,10,11 count all1,568 critical unordered endpoint pairs. There are no uncounted decomposable critical pairs. In contrast, exactly two ordinary endpoint pairs require34 clocks, one at grade8 and its grade11 complement, as independently checked from the complete pair histogram. This is a structural consequence of the independently checked finite catalogue, not an extrapolation or a new all-size supply theorem.

## Evidence and reproduction

check.py independently expands each primitive's untouched-coordinate product and compares it with the closed binomial formula. It compares all grade counts and aggregate cost counts against the frozen independent threshold audit, all grade/cost counts against the final producer certificate, and every critical count against the independent exact strict-pair histogram. The final producer's entire normalized primitive list and all-pair optimum digest must still agree with the original independent audit.

SOURCEPINS.json includes the final producer manifest, its exact mathematical source snapshots and successful final-revision portability log, and the unchanged independent audit manifest/result and portability receipt. Snapshots are data or text archives, never executed producer code. The checker verifies that reviewed source bytes match their final manifest entries. This supplement replays the completion argument; the main independent packet separately retains full native and dyadic authentication.

Run python -B reproduce.py for normal/-O equality, or python -B verify_manifest.py --relocate for a real manifest-only copy. Portability receipts are external to the frozen manifest, avoiding self-reference.
