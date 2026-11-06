# Uniform-incidence highest curvature

Status: **PROVED, independently audited in the stated scope.**

For a transitive m-bit spine and k>=1 core vertices which uniformly copy or complement those bits, the Hamilton-count Boolean polynomial has exact degree 2k whenever m>=2k. Its mixed difference on any consecutive 2k-bit block has exact 2-adic valuation 1+2v2(k!). The remaining core can be arbitrarily large but must have uniform incidences toward the spine.

Therefore two such contexts cannot satisfy H(child)=a H(parent)+c throughout the bit cube for any nonzero even integer a and bit-independent c, provided at least one sensitive count is positive and m is at least twice both sensitive counts. This includes length-dependent cores only at lengths satisfying that inequality.

The result does not cover positional or nonlinear incidences, changed spine edges, restricted bit languages, or the regime k>m/2. It does not exclude odd multipliers, controlled nonconstant corrections, or multiplier-one equal-profile substitutions. It proves no globally missing Hamilton count and no eventual spectrum threshold.

## Proof, dependencies and retained evidence

PROOF.md contains the full block-count bijection and valuation argument. Its only general dependency is the standard oddness of Hamilton-path counts of tournaments. No tournament bank, strong-order estimate, or conjectured coverage statement is used.

REVIEWED_PROOF.md retains the frozen input. STATUS_PROMOTION.diff records the sole status-line change. check.cpp and NORMAL.txt are unchanged producer controls on 38 contexts and 1,520 literal graphs.

audit/ preserves every entry of the independent audit manifest and the manifest itself. REVIEW.md gives the mathematical review. Independent controls cover 88 contexts, 4,256 literal graphs and 324 consecutive derivatives, with 768 graph counts also checked by separate literal DFS. Positive odd-multiplier controls and scattered-block/nonuniform counterexamples enforce the scope boundaries. Normal and optimized-Python/debug-C++ receipts match exactly.

## Reproduction

From this directory:

    g++ -O3 -std=c++17 check.cpp -o /tmp/uniform_curvature
    /tmp/uniform_curvature > replay.txt
    cmp replay.txt NORMAL.txt

From audit/:

    sha256sum -c SHA256SUMS
    g++ -O2 -std=c++17 count.cpp -o /tmp/uniform_curvature_independent
    python3 verify.py /tmp/uniform_curvature_independent > replay.json
    cmp replay.json NORMAL.json

See audit/REVIEW.md for the separate debug/optimized replay and detailed controls. The finite controls test the all-order bijection, rather than replacing it with extrapolation.
