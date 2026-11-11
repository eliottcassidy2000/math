# Independent review of the decorated H9 fibre

**REVIEW PASS**, 11 November 2026. This review was performed by the arithmetic
worker, separately from the author of the H9 fibre proof. It is a mathematical
and independent-algorithm review, not an external referee report.

Reviewed files: `../fibres/H9_FIBRES.md` and `../fibres/h9_fibres.py`.

The complete-core proof is valid: after removing leading I and trailing O,
the first O and last I allow 2^(l-2) paths with the apex internal, in addition
to the two endpoint-apex paths. Thus H9 forces l<=4. Independent direct
permutation enumeration of all seven trimmed words of lengths2,3,4 gives
precisely OOOI,OIOI,OIII at H9. This independent oracle counts permutations
whose consecutive edges occur in the literal adjacency matrix; it does not
use the inherited append formula or subset dynamic program.

The A response was independently rederived from the OI block map
`(H,P)->(H+2P+2,4P+2)`. The three coefficients and constants are exactly
`(100,53),(220,77),(508,149)` in `(c*4^b+6b+d)/9`.

For C, repeated OOI blocks give
`P=6(8^a-1)/7` and `H=1+[24(8^a-1)+28a]/49`. For B, repeated OIO
blocks give `P=5(8^a-1)/7` and `H=1+[10(8^a-1)+28a]/49`.
Appending a substituted core contributes `alpha*P+h(core)-1`; these identities
give the stated B/C formulas, including their linear terms. The first O
removes dependence on the old last-bit coordinate at exactly the required
point. No independent branch labels have been mixed.

The fixed-translation proof was checked separately. At equal exponents the
minimum nonzero coefficient difference is120. At adjacent exponents the
minimum absolute normalized difference is27. At exponent gap at least2,
`100-508/16>27`. The lower-order term is bounded by `6M+96`, so a
30000 difference requires `27*4^M-6M-96<=270000`. The left side is increasing
and equals442230 at M=7. A separate direct evaluation of all21 formula values
with b=0,...,6 yields no pair at difference30000.

The scopes are correctly separated: three decorated cores correspond to two
ordinary strong classes; the ordinary generating function counts only this
one-apex H9 subfamily. The source `C3 -> C3` is outside it. The result does
not classify every H9 tournament, all ordinary tournament classes, or all
possible 30000 translations of rooted counts. No correction was needed.

## Extension-content follow-up

The subsequent `../fibres/extension_content.py` was independently challenged.
For both minimal two-triangle sources `C3 -> TT_b -> C3`, b=0,1, every
one-vertex extension was counted by direct permutation enumeration. All64
cuts at b=0 and all128 cuts at b=1 agree with the quadratic polynomial
computed by the worker's path-cover method. The gcd of the actual counts
is respectively1 and3. The minimum/maximum extension counts are9/99 and
9/207, respectively.

The all-order mechanism is sound: a coefficient with its named endpoints
not meeting both triangles admits a free rotation of an untouched triangle,
so is divisible by3. The remaining cross-triangle pair has the two ordered
cover counts `1,9` when b=0, but `0,9*2^b` when b>=1. A middle singleton
cannot belong to either forward cover at b>=1 because it lies after the
specified first endpoint and before the specified second endpoint. The
reverse cover has both full triangles, with each middle vertex assigned
freely to either of the two paths. The internal triangle coefficient is−6,
so the separated case has content exactly3 rather than9.

The two ordinary strong H9 cores have actual extensions H17 and H19,
respectively; ordinal exterior padding preserves these witnesses. For the
separated-triangle sources, extending each triangle to S4 gives H25 after
two births. A two-source/sink extension also retains H9, so the gcd of all
two-birth extension counts is1 by gcd(9,25)=1. Both witnesses are needed for
this last inference; a single H25 witness only excludes divisibility by3.

The full all-H9 claim uses the inherited complete H<=9 strong-core language
in `../../uniform_return_metagraph_20261111/metagraph/FAMILIES.md`, in addition
to the newly checked extension mechanism. This review does not reprove that
separate inherited structural classification.
