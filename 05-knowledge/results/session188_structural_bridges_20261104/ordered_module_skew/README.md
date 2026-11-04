# Ordered-module cover skew and cyclic prime quotients

PROOF.md gives a support-preserving cover injection with a quantitative positive
skew for every ordered module A->B of size at least three, in any exterior.
Together with the marked-pair quartic identity this proves quotient-converse
separation for all five-vertex quotients with nontrivial children and at least
one pair child. It also gives an unbounded cyclic strong-prime quotient family
with arbitrary other children, under the stated opposite-child size condition.

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

The checker uses the adjacent marked_pair_quotient/verify.py primitives, whose
hash is pinned. It checks 311 ordered-module configurations, 451 literal cover
injections including 24 covers with other module visits, all 5,120 five-vertex
quotient/root classifications, fifteen ten-name prime-five batteries, cyclic
prime structure through order15, and actual quotient-order7/9 batteries.

Publication files: PROOF.md, README.md, verify.py, receipt.json,
receipt_optimized.json, SHA256SUMS. Independent review passed in the adjacent ordered_module_skew_audit packet.
