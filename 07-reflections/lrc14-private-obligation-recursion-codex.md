# LRC14 private-obligation recursion

This session tried to reframe the last LRC(14) gap-side obstruction after the
fresh THM-526 arc-width lemma.  The key move was to stop asking whether
divisibility pressure is globally high, and instead ask whether a parked runner
is **privately responsible** for one of the q-grid obligations in THM-523.

The exact classifier found a useful split.  In a scout of `103` primitive
q-covering 13-sets, `94` were certified by THM-526 arc-width, and the `9`
arc-width residuals all had a parked runner with a private q-obligation.  No
LRC breaks appeared.  This is not a proof, but it says the residual has a much
smaller shape than "all covering sets": it is a small-m parked-runner packet
where deleting the parked speed uncovers a specific q-grid.

The other useful correction is scope.  The older `{1,...,11,13,84}` row with
`M=7/89` is the natural champion when the non-parked core remains
unit-residue-complete modulo 14.  But in the broader THM-523 q-covering sense,
the principal row

`{1,2,3,4,5,6,7,8,9,10,11,12,182}`

is closer to the floor:

`M=14/183`, with slack `13/2562`.

Its binding pair is `(1,182)`, so `D=183`, `j=14`, and `14j-D=13`.  It misses
unit residue `13` modulo 14, which is exactly why it was invisible to the
unit-residue-complete story.  That does not damage the old result; it separates
two proof surfaces that had been touching.

The tournament audit was deliberately modest.  Vertices were sampled covering
rows; the exact hardness tournament is transitive because it is ordered by
`M-1/14`.  The interesting comparison was against a blocking-height proxy
based on q-cover counts, parked private obligations, arc-width margin, and unit
coverage.  It agreed on about `65%` of pairwise orientations and flipped on
`35%`.  So "dominance grows with blocking height" is directionally true but not
sharp enough.  The missing coordinate is private ownership of an obligation at
the parked runner plus the forced THM-524 binding flank.

The proof target that now feels most alive:

> If `S=A union {14m}` is q-covering, not arc-width-certified, and `14m` has a
> private q-obligation, prove that the forced binding crossing has index
> `j >= D/14`.

For principal single-drop towers this is essentially the closed-form monotone
collapse already found: the infinite tower reduces to `k=1`.  For mixed cores,
the hoped-for recursion is: either no parked runner has private obligations and
we delete/recurse, or a private obligation exists and forces a small finite
flank law.  This is the most concrete "keyhole" I found in the last little bit:
not runners, not raw residues, not global dominance, but the private q-debt
that a parked runner alone pays.
