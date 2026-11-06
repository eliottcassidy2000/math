# Session 208: native tree-capacity parity

Independent elementary proof: canonical tree weighted degree is the original degree plus twice internal transit. When all graph degrees are odd, the odd-capacity edges form the tree's unique V-join, and minimizing their number over spanning trees equals the minimum V-join size.

For the actual labelled order-four tournament path-complement metagraph, the sharp minimum is40, not the general64/2=32 bound. The ordinary-count strata H=1 and H=3 form an independent40-vertex set; an explicit40-edge odd forest extends to the supplied63-edge spanning tree. The maximum matching size is24.

- PROOF.md: general identities, sharp finite result and scope.
- check.py and check.out: exhaustive64-state construction, cut-based verification and explicit certificate.
- independent_check.py and independent.out: independently generated fixed-path graph and route-based load verification.
- REVIEW.md: mathematical audit.

Run `python3 check.py` and `python3 independent_check.py` with assertions enabled. Both use only Python's standard library. Their outputs should match the recorded files. `sha256sum -c SHA256SUMS` checks the frozen payloads.

This work is inspired by the fundamental-cut capacity identities in Hajiaghayi, *Logarithmic Embeddings into Probabilistic Spanning Trees for Distances and Capacities*, arXiv:2611.05189v1, https://arxiv.org/abs/2611.05189. Its logarithmic theorem is not assumed or independently certified. Native labelled degree equals ordinary H; quotient, loop and capacity conventions must not be substituted silently. No literature-priority claim is made.
