# A divisor-sensitive order bound for actual rooted profiles

**PROVED as a cited corollary, independently audited.** This applies Busch's published strong-tournament minimum to the physical edge extension of the root component. It is an order bound, not a profile-existence theorem.

For a realized positive-start pair (B,s) with B>s, delete later singleton strong components. The remaining tail count q divides gcd(B,s). If S(t) is the maximum strong order allowed by Busch's minimum and G(t) the maximum reduced order allowed by the inherited component minimum, the same induced representative has order at most

    max over q dividing gcd(B,s), S(B/q)>=4:
        S(B/q)-1+G(q).

This is at most G(B)-1 and can be combined with the elementary bound s+1+floor(log2 B). The separate B=s boundary has bound 1+G(B).

At the staggered base X=137, all eighteen candidate pairs have a complete order cap of at most 10, compared with the earlier elementary maximum 33. This does not assert that any requested pair exists, or that the cap is attained. F(8,137) and the eventual threshold 23 remain OPEN.

Read PROOF.md, the exact dependency provenance in dependencies/PROVENANCE.json, and audit/REVIEW.md. The reviewed original is preserved, and STATUS_PROMOTION.diff changes only the status line.

Reproduce:

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify > replay.json
    cmp replay.json NORMAL.json

The producer and independent audit check all 33,867 labeled root-zero tournaments through order 6. The independent audit also verifies minimum-product partitions and 259,649 numerical pairs using exact integers. These controls supplement the cited all-order proof.
