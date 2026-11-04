# Independent audit: known six-pair reconstruction

PASS. No mathematical correction was required. The proof is specific to a known partition into six two-name modules; the reconstruction count is not an exhaustive bank claim or an unrestricted twelve-name result.

## Proof review

The four-name endpoint expansion is correct. Summing reversed candidate Hamilton paths cancels cubic terms; the remaining terms reduce to the displayed constant, linear and ordered quadratic sum. Across all four-name subsets of a five-name tournament, constants and quadratic terms occur an even number of times and each anchor incidence occurs three times. This proves the degree-parity identity.

For six quotient names, deleting j therefore gives h_ij=d_i+A_ij. Its signature cannot vanish: rowwise constancy would make every vertex a source or sink, impossible for six names. Equality of two signatures makes adjacency differences constant in each row; tournament antisymmetry makes those six row constants identical. The two possible quotients are exactly whole converses. The explicit reconstruction formula recovers both candidates, and a nonzero shared-cubic contrast distinguishes them using the same coarse query battery.

The 94 low-mask queries and their reflections are disjoint, giving exactly 188. Every required eight-term finite difference is included. Coarse queries ignore internal pair directions exactly, by run reversal. After recovering the quotient, the six isolated pair queries distinguish H=1 from H=5 in their four-name cores. All other pairs remain modules. The counts 194 adaptive and 380 fixed queries are valid nonoptimal upper bounds.

The ambient reflection step is correct: only the six core attachment bits are reflected. The two actual fourteen-name cores, including both ports, are compared inside the same rooted exterior K. Even-core transfer uses one query-independent odd multiplier; subtraction and finite differences therefore preserve each residue 4h_ij. The pair decoders similarly retain their gap four in every fixed ambient context. No reflection of the exterior attachments, changing exterior, or unknown pair partition is assumed.

The producer's integer bounds are safe: every exact physical count is on at most sixteen names, and 16! fits in uint64. All assertion-style checks remain active under Python optimization.

## Independent finite checks

The audit imports no producer code. It enumerates literal paths to prove all 1,024 five-name endpoint instances, checks every entry of the retained 32,768-signature bank and its 16,384 exact converse classes, and independently constructs the 188-query battery.

A separately written predecessor subset DP recomputes all retained physical data, additional quotient samples and a complete 188-query strong ambient signature. It checks all 64 pair-bit states under six decoders. The exact number of distinct physical tournaments is in the paired receipts. Both normal and optimized audit modes agree and bind the source and full-data hashes. These controls validate implementation; the two-million-state injectivity conclusion follows from the structural proof.

Exact publication allowlist: REVIEW.md, audit.py, count.cpp, normal.json, optimized.json, SHA256SUMS. No compiled binary or producer import is required. The audit expects the source's retained literal.json and h_vectors.bin at their bound paths.
