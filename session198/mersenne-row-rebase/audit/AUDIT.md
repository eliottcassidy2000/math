# Independent audit: balanced witnesses for 3*2^b-1

Verdict: PASS.

- PROOF.md SHA256: 9973a3bf57d3317b416023d16ba7334e51ae16f2d9dc1f782f628bab440218a7
- verify.py SHA256: 15cf5f40b6baa5b88b820d8423171351ddd3df307cf1e620343dbc99b6f1fed7
- NORMAL.json and OPTIMIZED.json SHA256: 184f13329f2c09b2f01643bf1d3ae15214bb1f51e7382f7a56ed93d043f29761

Both displayed physical word families are correct at all orders. After the
first O^i J block the profile is (2x+1,1,2x+1,1). The stated OO and OZO
prefixes are exact. Iterating Z by y=2^j gives
D_y=D-E+(y-1)U-(y-2)d and U_y=yU-(y-1)d, with E_y=d;
these formulas follow by induction and supply both displayed source profiles.
The final J leaves D and U consecutive, proving the claimed primitivity.

Before tying the two run lengths together, the balanced readouts are exactly
48xy-18x+18y-1 and96xy-30x+30y-1. Setting x=y cancels the linear terms and gives
every even exponent b>=6 and odd exponent b>=7. The separate b=5 word yields
(19,0,38,5), hence95. There is no gap between the three exponent cases.

The chronological JZ return is the physical M2 extension and has profile
(D+2U,D,D+2U,D). Its unique rooted Hamilton path closes via the new final
vertex, so it is strong C1; no unproved balanced-source converse is required.
The exact word lengths give b+4 vertices in every case, including nine atb5.
The residual-row identification is a direct arithmetic equality. It does not
assert coverage of the other rows or cofinite odd coverage.

Isolated normal and optimized replays reproduce both supplied receipts, and
all four review-manifest hashes match. Additional independent checks certify
the two bivariate degree-(1,1) polynomial identities on an unisolvent four-point
grid, and reconstruct all twelve source/final graphs for b5 through10. A
separate forward subset recurrence, rather than the packet's backward
recurrence, verifies complete source profiles, root-start counts, balanced
profiles and exactly one Hamilton cycle. Final counts are95,191,383,767,1535,
3071, at orders9 through14. See ../mersenne-independent.py and
../MERSENNE_INDEPENDENT.json.

No correction to the frozen mathematical statements was needed.
