# Session245: physical spectrum closure attempts

**SESSION COMPLETE, 10 November 2026. Independently audited. Chan–Pak remains OPEN.**

The goal is every positive odd Hamilton-path count except 7 and 21. This session
does not close that conjecture. It produces three exact physical tools and
identifies precisely why their apparent shortcuts do not prove coverage.

## New proved mechanisms

1. **A small pure-count switch beyond the word-family divisor.** Two strong
   seven-vertex tournaments have counts 69 and 71 with identical rooted
   endpoint, response and deletion counts (13,77,13). The deletion count 13 is
   coprime to 3. Full seven-coordinate interfaces are retained and independently
   recounted. [Proof and exact contextual laws](interfaces/PROOF.md).
2. **Physical response cancellation by reversal.** An opposite-facing module
   with a transitive separator cancels the response error that could not
   cancel in the earlier similarly oriented chain. Exact strong examples have
   gap 812 at order 19 and gap 14812 with deletion 47515 coprime to 3 at order 22.
   The latter is proved by constituent enumeration plus the composition law,
   not by a full order 22 subset recount. [Proof](interfaces/PROOF.md).
3. **Synchronized outer blocks and exact defects.** Every fixed O/Z/J middle
   supplies an explicit quadratic in 2^k. For binary middles the linear term
   cancels exactly for palindromes. Three half-word norm formulas and a complete
   target inverse retain actual graph recipes. Root-changing middles give,
   for example, the all-k family 436·4^k−37. Mersenne coverage and the original
   binary block law are explicitly inherited. [Proof/compiler](descent/PROOF.md).
4. **Simple constant-run spectrum.** Every individual letter in the inherited
   k-apex tape is a sum of commuting projections, with square-free annihilator
   ∏(x−j), j=0,...,k+1. Exact run interpolation has no polynomial prefactors;
   the independent-star family's leading coefficient is k!.
   [Proof](forest/PROOF.md).

## The remaining barriers are now exact

A mixed h13/h37 chain has coefficient gcd exactly 2 at every length. However,
its 2^k scalar outputs have minimum spacing 50·37^(k−2) for k>=2. This is a real
escape from the growing-divisor obstruction but still not an interval.

For a fixed apex graph and alphabet, at most r constant runs supply only
O((log N)^r) scalar values through N. Even r=o(log N) gives N^{o(1)} values.
This uses the inherited finished-word deletion kernel; it excludes only that
run budget, not the unrestricted fixed-apex language. [Proof](forest/RUN_BUDGET.md).

The all-target obligation remains a physical selector or well-founded
lifting that fixes the requested TOTAL count, rather than only a gap, modular
residue, coefficient gcd, endpoint count, or formal predecessor.

## Verification and provenance

[Independent audit](audit/AUDIT.md) passes the scoped proofs, literal small
closures, arbitrary apex/signature controls, and ordinary/optimized replays.
The aggregate additionally checks all stored receipts and failing sentinels.

Inherited sources are Session199's physical grammar, adjoint and binary
coefficients; Session221's fixed-apex transfer and SCC kernel (crediting
Sessions215/216); THM-4748's Mersenne fan and two-arc formula; and Session244's
seven-coordinate gluing. None is relabelled a new coverage theorem.

From this packet directory:

    python3 -I -B verify_all.py > REPLAY.json
    cmp AGGREGATE.json REPLAY.json
    sha256sum -c SHA256SUMS

All computations are finite checks of implementations or finite subclasses.
The universal statements are proved separately in the linked documents.

Canonical statements: [THM-5117](../../../01-canon/theorems/THM-5117-physical-cut-switches-reversal-cancellation-and-primitive-sparse-cubes.md), [THM-5118](../../../01-canon/theorems/THM-5118-synchronized-source-blocks-and-palindromic-defect-inverse.md), [THM-5119](../../../01-canon/theorems/THM-5119-commuting-gap-projections-and-fixed-apex-run-budget.md). Session245 has ended.
