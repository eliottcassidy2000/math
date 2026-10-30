# Independent audit: PRIME bank ranking

**Status: PASS, no mathematical blocker found.** Reviewed frozen Session 113 `prime_ranking/PRIME_RANKING.md`, implementation, literal ring extractor, and checker. No Git writes. Approval applies to named binary O_k bank orbits, relative to the explicitly inherited ideal-module theorem, not full graph isomorphism or an open-problem solution.

## Proof review

- The even valuation of nonzero principal scalar displacements and odd valuation of nonzero clock displacements force trivial intersection; the stated orders multiply to the principal-unit group order. Residue representatives and scalar/clock factors therefore give the required coordinate bijection.
- The leading-residue orientation is preserved by both generators and by translations inside each applicable ideal. Thus every ideal cut contains one common fixed origin. There are no omitted affine signed offsets in the cut dictionary.
- Quotient depth r-d leaves scalar exponent modulo p^floor((r-d-1)/2), establishing the period-cap dictionary, independently for each named residue and valuation.
- Exact ORIGINAL period p^v is indispensable when v<=q: these are distinct labelled fixed words, not one packet-period category. The implementation correctly separates these counts. For v>q, packet period determines original period and the inherited necklace routine applies.
- Stabilizer and surviving-cut-mask updates are sufficient: all representatives of one exact-period class have the same stabilizer and same remaining cut compatibility. Conditional suffix orbit counting is legitimate; no independent scalar shifts are assumed.
- Terminal empty mask is exactly the inherited prime condition. Weighted cumulative counts and binary search preserve lexicographic representative order and give mutually inverse rank/unrank.
- The full 2^(2k-1) mask set is polynomial in N=p^k-1 for p>=5, not polynomial in k. Counts and threshold matrix entries have polynomial bit size. The complexity claim is appropriately parameterized.
- `bank_rank.py` is byte-identical to the checked Session 112 dependency.

## Independent evidence

`independent.py` constructs the common coordinate gauge in the actual p=5,k=2 ring, then recursively enumerates canonical blocks using explicit lists of group elements and their literal stabilizers. It does NOT use the new symbolic stabilizer recurrence, cut masks, period classes, or completion counts to enumerate/classify orbits. At each leaf it evaluates literal signed ideal-cut equations.

Results:
- All **692,224** named-group canonical orbits enumerated.
- Exactly **665,856** pass the literal prime predicate.
- Independent sorted-list indices match production rank and unrank at **46** boundary and seeded random ranks, including first/last and their neighbors.
- Six invalid rank inputs rejected.
- 36 length-25 thresholds check original periods 1 and 5 hidden inside packet period 1 at stride 25.

`thresholds.py` additionally exhausts **all 1,539** length-nine binary threshold/stride cases (strides 1,3,9). Expected values come from literal rotated words and original periods. This deliberately uses p=3 only as an abstract period-arithmetic control, not as a ramified-ring theorem.

The supplied checker passes under normal Python and `python -O` (2,101 explicit checks each). Added raw-ring families `(p,k,Delta)=(7,2,14),(11,2,22),(5,3,15)` pass representative, generator, every-cut partition, total, canonicality, and literal-prime rank checks. Larger supplied total-only controls also pass; these are not large rank/unrank tests.

Published independent receipts: `independent_results.json`, `threshold_results.json`, `additional_rings.json`. Producer receipts live beside the [ranking implementation](../prime_ranking/PRIME_RANKING.md). Redundant duplicate producer receipts are omitted.

Publication makes the two independent scripts location-independent and mechanically replaces assertions by explicit exception gates, preserving each condition. Run `python independent.py` and `python thresholds.py` here; both also support `python -O`.

## Boundaries

This audit does not reprove the inherited ideal-module characterization for all tournament modules. No claim that free action is equivalent to primeness is used. Exhaustive orbit enumeration is p=5,k=2 only; the general theorem rests on the reviewed proof and inherited dependencies. No universal full-normalizer quotient or major-conjecture closure is established.
