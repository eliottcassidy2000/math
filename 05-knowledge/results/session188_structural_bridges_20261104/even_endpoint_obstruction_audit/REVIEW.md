# Independent review: even endpoint-parity obstruction

**PASS. No mathematical corrections.** This review binds the frozen source proof and exact executable pins in `normal.json`. It is a criterion obstruction, not a full-response collision.

## Mathematical scope

The physical quotient has rows `(238,4,240,118,195,146,34,74)`. Independent literal enumeration of all vertex words recovers exactly 249 Hamilton paths, start vector `(87,13,13,42,31,19,18,26)` and end vector `(13,65,29,22,13,29,52,26)`. The two parity masks are both 55. Checking all proper subsets of size at least two shows that the quotient is prime.

The smaller even-order census covers all labelled tournaments at orders 2, 4 and 6, without isomorphism reduction, and has no endpoint-balanced case. Thus order 8 is the first positive even order where this sufficient endpoint criterion can fail. No larger-order inference is drawn.

Section 2 correctly invokes the inherited complete coarse modulo-four law only for nontrivial children. That law says all coarse contrasts vanish modulo four here; it does not say any exact response vanishes. The bare physical query `(0,1)` gives 423 and 1313. The full bare-converse quartic contrast is zero by the inherited whole-converse symmetry, but cubic contrast is not zero. All coefficients, including these zero quartics, are retained.

The pair-inflation shared-cubic battery uses uniform masks on each physical pair. Its three Boolean variables have distinct roles `a_i,a_j,c_j`, even though the latter two refer to the same pair. The eight signed actual-query values give exactly -56,897,620, congruent to 4 modulo eight. The even-subset endpoint matrix is independently recovered from literal path words. This is consistent with, and does not weaken, the earlier eight-pair separation theorem.

## Independent computation

The auditor imports no producer helper.

- It enumerates every directed word on every nonempty subset of the eight names. This independently recovers the complete Hamilton-path list, start/end vectors and every even-subset endpoint contribution.
- The order-2/4/6 census enumerates Hamilton words first. For each word it enumerates all tournaments completing its forced adjacent edges, toggling endpoint parity masks. This differs from the producer's per-tournament subset recurrence and reproduces the complete histogram.
- A separate C++ predecessor recurrence counts all 2,818 physical augmented tournaments. The producer uses forward extension. Every retained integer response agrees, including all 1,369 bidegree-at-most-(2,2) bank coordinates, 32 full-mask controls and eight pair-inflation battery queries on both graphs.
- The audit independently reconstructs all finite-difference coefficients and all 32 full masks. The nonzero contrast counts by total degree are `(0,16,56,444,0)`.
- All pair interiors, physical names, masks, repeated-index roles and coefficient signs are independently checked.

Both normal Python and optimized Python produce identical receipts, with explicit checks active in both modes. C++ counts are bounded by 18!, which fits unsigned 64-bit integers. Compilation uses the installed compiler and installs nothing.

Frozen audit allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.
