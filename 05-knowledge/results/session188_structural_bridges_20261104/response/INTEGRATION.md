# Incoming parity theorem closes the single-module branch

The fresh main checkpoint [`76b3b1a77d124100f03026420f6b821533cdec15`](https://github.com/eliottcassidy2000/math/commit/76b3b1a77d124100f03026420f6b821533cdec15) contains [the cover/current parity certificate, E1-E4](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md) and its [same-session review](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/REVIEW.md). Both were read directly at that exact commit during this session.

## Exact implication

Put D_uv=(C_uv-C_vu)/2 modulo 2 and h_uv=[a_u a_v c_v](F_T-F_(T^op))/2 modulo 2. Incoming E1 proves, for every n>=3 and every anchor v,

    sum_(u!=v) D_uv + sum_(v->u) h_uv = 1 modulo 2.

Its second proof partitions the endpoint contributions of even-support paths in T-v. Prepending v produces odd paths starting at v. Complementing the last endpoint direction separates odd paths ending at v from one contribution H(T[A]) for each nonempty even subset A of T-v. Each H is odd and there are 2^(n-2)-1 such subsets, an odd number. This proves the compatibility obstruction at all n>=3, not merely at twelve. At n=2 the constant is zero, and incoming E3 supplies a separate literal query.

Consequently every nontrivial tournament differs from its same-named whole converse by 2 modulo 4 at some actual two-port query. Incoming E4 gives the ports the replaced root's exterior incidences and applies the inherited modulo-four substitution transfer. The odd H(Q) factor preserves residue 2. Thus **no single nontrivial module converse is hidden in any fixed ambient quotient**.

This supersedes the earlier open status of the endpoint-balanced strong prime twelve-core branch. It is an inherited all-order closure, not a conclusion of the new quartic census. The displayed transfer in E4 is a congruence modulo four; it does not assert exact integer cancellation of a strong quotient.

## What our packet still contributes

The quartic normal form, deletion-capacity identities, sharp arbitrary-comparison threshold at five old vertices, and exponential transitive-exterior quartics remain valid. They describe the exact integer top-degree response and show why an internal module reversal cannot be reconstructed by discarding quartics. They are complementary to the incoming lower-degree modular separator; they are not needed to close the now-resolved single-module question.

The finite order-five quartic classification is a separate extension. Its four-member exceptions are only quartic-tensor collisions and are explicitly separated by actual full-response queries. They do not contradict E4.

Unrestricted reconstruction from the full named response, exact collisions between arbitrary unrelated tournaments, and simultaneous changes of several modules are not settled by E1-E4 or this packet. Any next search must preserve those distinctions rather than relabel the resolved single-module residual as open.
