# Independent audit: exact continuous bounded-cone certificate

**PASS.** The coefficient certificate proves the claimed sign for the whole continuous simplex of nonnegative real lower weights at each order1 through7, for every h>=2. It is not an integer-weight census or an extrapolation to larger orders.

## Independent barycentric construction

The source expands in lower weights and homogenizes with a slack coordinate. The independent audit instead works directly with the simplex's extreme profiles. For L=n!-1, these are t^(n)+L t^(r) for each lower same-parity generator, together with t^(n) at the slack vertex. Every permitted profile is their barycentric combination. Since the cover contrast is cubic, its complete homogeneous coefficient at each three-index monomial is the corresponding trilinear polarization times its exact permutation multiplicity.

This construction independently reproduces every retained coefficient, including every zero, using standard-library exact rational arithmetic. Every coefficient is nonpositive. Nonnegative barycentric coordinates therefore give the asserted sign on the entire simplex. Checking only the diagonal vertex values would not suffice; this audit checks all mixed cubic coefficients as well.

The n=1 zero-budget case and the one-coordinate cases are handled directly. The exponent cutoff is valid because multiplication cannot lower any variable exponent. The high-h tail vanishes once the numerator's minimum total degree exceeds the available three profile degrees; the retained boundary is conservative and safe. Thus no higher cover coordinate is omitted.

Both normal and optimized audit runs are byte-identical. Fresh normal and optimized source replays also match the frozen receipts exactly. There is no floating-point approximation, optimization solver, or removable assertion.

## Scope and reproduction

The result is a bounded-order continuous-domain theorem for the stated ordered-cover contrast. It neither classifies tournament-realizable profiles nor establishes an all-order sign or general named reconstruction. The independently proved covered-ambient substitution implication can be applied separately when its actual-tournament hypotheses hold.

Run `python3 verify.py > normal.json`, `python3 -O verify.py > optimized.json`, and compare. The audit imports no producer modules. Frozen allowlist: `REVIEW.md`, `verify.py`, `normal.json`, `optimized.json`, `SHA256SUMS`. Exact source pins, every simplex vertex profile, and the complete coefficient table are retained in both audit receipts. Redundant producer-replay copies are excluded from the allowlist.
