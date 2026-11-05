# Signed parity-class coverage and physical fiber repairs

Status: PROVED and independently audited. The proof snapshot retains its pre-audit qualifier; `audit/AUDIT.md` records completed review of its exact SHA256.

- Every reduced interior rational slope with odd denominator covers its full signed parity class on actual strong unit-interface sources.
- Explicit populated rational fibers permit genuine target-preserving C2 J-swaps, including common target 7649 on two directly counted 15-vertex graphs.
- The odd-numerator 5/9 fiber has an explicit 34-vertex pair with common target 2972894005, certified by source words, the independent g recurrence, and the physical graft theorem.
- A prescribed positive target still fixes an unproved arithmetic index. Universal odd Hamilton-count coverage remains OPEN.

Reproduce from this directory:

    python3 verify.py > /tmp/signed-repair-normal.json
    python3 -O verify.py > /tmp/signed-repair-optimized.json
    cmp /tmp/signed-repair-normal.json /tmp/signed-repair-optimized.json
    python3 audit/check_signed.py
    python3 -O audit/check_signed.py
    python3 audit/signed_literal_graphs.py
    python3 -O audit/signed_literal_graphs.py

The author verifier checks 106 slopes through denominator 31 and 1111 compressed recipes. Independent arithmetic checks cover 297 slopes through denominator 53; an independently implemented subset DP rechecks both 15-vertex C2 graphs. All tests use explicit exceptions. The 34-vertex certificate is not claimed to have a direct Hamilton-path census.

Frozen dependency proofs, exact recipes, author receipts, and independent audit checks are included. No named external conjecture is assumed, and no claim about Diophantine integer-ring definability is used to infer existence.

The separately audited [short block transport addendum](../signed-block-transport/README.md) gives a second physical proof and an exact suffix normal form for the 5/9 response fiber.
