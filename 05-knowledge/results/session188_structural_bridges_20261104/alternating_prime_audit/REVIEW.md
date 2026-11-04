# Independent audit: alternating-hub quotient theorem

PASS. No mathematical correction was required. The audit independently checked the proof, source receipts, structural cases and every retained physical query. The producer files were not edited.

The primeness proof is valid at every r>=1: two chain vertices force their entire interval and hence the hub, while the hub together with any chain vertex forces both endpoints and then all chain vertices. In the latter step the stated incidences are correct even at an endpoint, where that endpoint is already present.

Deletion of the hub or either chain endpoint exposes a global cut. Deletion of an interior chain name exposes its two chain neighbors; these have equal parity, ensuring uniform hub adjacency. The marked-pair identity and ordered-module inequality apply to their inflated unions without assumptions on other child interiors.

For r>=2 the global-cut exterior already has at least four names with merely nonempty children. The interior deletion needs the two neighboring children to have total size at least three. The main all-nontrivial-child hypothesis implies this, but is stronger than necessary. At r=1 the two-singleton exterior is exactly a zero boundary of this quartic, not a complete-response collision. The proof makes no outer-ambient or arbitrary-child-without-a-pair claim.

Independent checks enumerate all possible proper modules through quotient order 17 and check all 80 deletion structures. They independently recover all 17 retained grouped quartic coefficients from literal integer subset-DP queries: 14 positive and 3 zero. The 272 attachment queries per graph retain every finite-difference term. Both audit modes agree, and fresh normal/optimized producer replays exactly match the frozen receipt.

Files: REVIEW.md, audit.py, normal.json, optimized.json. The checker imports only independently written primitives from ../ordered_module_skew_audit/audit.py, which must remain available. Source and dependency hashes are bound in each receipt. No finite result is used to infer the all-order theorem.
