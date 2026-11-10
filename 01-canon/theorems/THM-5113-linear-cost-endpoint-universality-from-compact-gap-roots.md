---
id: THM-5113
title: "Every positive endpoint count has a strong logarithmic-order realizer with total count at most three times the endpoint and controlled response"
status: "PROVED, independently audited. Endpoint universality is inherited; quantitative full-state control is new here. Chan–Pak remains OPEN."
session: Session244-even-gap-anchoring-20261110
source: 05-knowledge/results/session244_even_gap_anchoring_20261110/selector/ROUND3.md
audit: 05-knowledge/results/session244_even_gap_anchoring_20261110/audit/REPORT.md
---

# Efficient endpoint realization with its full physical sidecar

For every N>=1, a strong tournament T_N has a distinguished vertex v with
exactly N Hamilton paths ending there, trivial root stabilizer, and

    |T_N| <= ceil(log2 N)+3,
    2N+1 <= H(T_N) <= 3N,
    B_v(T_N) <= 4N and B_v(T_N) <= 4H(T_N)/3.

Endpoint existence and trivial-stabilizer realization are inherited from
THM-5001-every-positive-endpoint-quotient-has-a-small-strong-realizer.md.
The present construction quantitatively controls all four rooted coordinates.

Let g(1)=g(2)=1 and use the inherited binary response selector

    g(2m)=g(m),  g(2m-1)=g(m)+m-(m mod2), m>=2.

For N>=2 build the canonical source with state(g(N),N mod2,N,1), append
OZ, and mark the O newborn. N=1 uses the directed triangle. With epsilon=N mod2,
the resulting full state (H,endpoint,response,deletion) is exactly

    (g(N)+2N, N, 2g(N)+2N-epsilon, g(N)+N-epsilon).

The marked vertex has a unique outneighbor, the final Z vertex, which itself
has the old root as its unique outneighbor. Fixing the marked vertex therefore
fixes these two vertices and then the remaining transitive chain, proving
trivial stabilizer. The last Z makes the graph strong.

The lower bound is attained exactly when N is a power of two. The upper bound
is attained at N=1 and N=2^k+1, k>=1, within this constructor. No global optimum
is claimed away from the inherited lower-bound equality slice.

The full H image is the compact-gap endpoint image already analyzed by
THM-5114. Thus even efficient complete endpoint control does not select every
full count. Keeping this root and appending t>=1 Z moves gives exactly

    H_t=(2^t+1)(g(N)+N)-epsilon.

For even N this is composite. For odd N, H_t+1 is divisible by the odd
number 2^t+1. Hence this specific continuation never reaches a Mersenne-prime
target above3. Another root or operation is required; no global Mersenne
exclusion is asserted.

[Proof and standalone literal checks](../../05-knowledge/results/session244_even_gap_anchoring_20261110/selector/ROUND3.md).
