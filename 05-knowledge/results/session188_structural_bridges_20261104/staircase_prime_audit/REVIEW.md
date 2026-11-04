# Independent staircase-prime review

Status: PASS. No mathematical correction required. Producer proof d840a72a4d8a2ff2be45cdb167a07e7f76ab29da6aa7ecb7ccd56cf1bf7f505f and checker 71e2ee4d2ba7723538d82cd0eb30e7c8263fd5bcfb125eb1763b703fb8bd1a49 are pinned in both independent receipts.

## Analytical review

The cyclic first-block reversal equals the stated two-chain staircase convention. The displayed cycle is directed, including its final edge. For primeness, two A-names force a B-name and two B-names force an A-name; a mixed pair forces B_0 and B_r. The latter endpoints force all A-names, which force every intermediate B-name. This works at r=1 without an empty-case inference.

Every deletion is correct with physical labels retained. In particular, deleting B_j exposes the directed pair A_j -> A_(j-1), not its reverse. Deleting B_0 leaves A_0 a sink, and deleting B_r leaves A_(r-1) a source. A module of the deleted quotient lifts to the union of its complete child blocks. The marked-pair exact coefficient and ordered-module injection therefore apply without any restriction on other child orientations.

The positivity hypothesis is the physical size of the exposed module, not quotient primeness by itself. For pair deletions its two child orders must sum to at least three. At endpoint deletions and r>=2 the full exterior already has at least four names. At r=1 both exterior children singleton give a zero quartic. The analogous two-singleton module in a larger quotient also gives zero: the pure terminal/initial-run swap is then a bijection, since the two prescribed endpoints already use both module names. These are zeros of the chosen coefficient, not full-response collisions.

The sixteen-query finite difference is on four disjoint groups of physical bits, and is applied to actual Hamilton counts. Neither cubic alias terms nor quartic terms are discarded. There is no replacement of an arbitrary child by a pair, no arbitrary strong outer-ambient lift, and no conclusion about all original-module families.

## Independent executable checks

`audit.py` imports no producer or response checker. It constructs U by reversing the first r names of a cyclic tournament, reconstructs physical blow-ups independently, and uses predecessor subset dynamic programming for path counts. It checks:

- Strong spanning cycles and exhaustive module exclusion through quotient order 17.
- All 80 displayed deletion structures through order 17, with exact ordered sides.
- All 24 producer batteries, including 18 positive and 6 zero-boundary batteries.
- Every one of 384 actual two-port query values on each graph, independently recomputed.
- Each finite-difference quartic against a separately computed ordered-cover skew.
- Exact normal/optimized producer replay and frozen source hashes.

Both independent modes have byte-identical receipts. The finite checks support the implementation; the all-order result comes from the analytical proof above and the previously reviewed marked-pair and ordered-module lemmas.
