# Independent audit: path-locked cycle moves

Verdict: PASS, with no unresolved mathematical blocker. Audited producer REPORT.md and check.py in /tmp/session56-path-locked-cycle-moves. Producer's initial bare-assert weakness was reported and has been corrected to explicit exception checks before the frozen version.

## Mathematical audit

- Balanced disagreement arcs decompose into directed edge-disjoint cycles; reversing those supports reaches the second labelled orientation exactly.
- A chord of a directed r-cycle yields two smaller directed reversals in the order determined by the current chord direction. Each recursive call restores all off-cycle edges. Their lengths sum to r+2, so the induction gives at most r-2 induced reversals. Applying this to a disagreement decomposition gives M-2q.
- The sharpness orientation of an induced cycle directs all outside-to-cycle edges inward and orders outside vertices acyclically. Thus no directed cycle leaves or enters the outside part; the sole directed simple cycle is the chosen cycle. Every balanced disagreement is empty or that whole cycle. This proves a two-state fibre and the exact universal threshold, not merely a counterexample in an induced subgraph.
- For a forest lock, an induced k-cycle in the complement would force k(k-3)/2 forest edges on k vertices. For k>=5 this is at least k, contradicting the forest bound k-1. This holds for every forest, regardless of degree.
- Fixed full labelled scores and locked directions are equivalent to fixed labelled residual scores. Every move avoids the locks. Degree-determined internal properties survive; Hamilton-path counts and their residues need not.
- The n=5 sharp pair, full scores (3,2,1,2,2), residual scores (2,1,0,1,2), and two-state argument are correct. The two HP counts are 11 and 13, hence unequal modulo 8. The fixed-path residual graphs for n<=4 are forests, proving minimality.
- Hodge identities, full complete-graph dimension, and connected residual-graph dimension (n-1)(n-4)/2 for n>=4 are correct. Cube-feasibility and native-operation limitations are explicitly retained.

## Independent finite verification

The independent audit.py finds cycles as connected induced 2-regular subgraphs, rather than producer permutation enumeration, and determines components by BFS rather than union-find. It reproduces all producer fibre counts, complete size histograms, induced-cycle counts, and move-edge counts for n=1..7. It independently checks every graph and orientation through n=5, including exact two-state sharpness fibres for every cyclic graph. It independently recomputes the sharp HP pair and absence of residual directed triangles.

The independent replay.py loads only function definitions and imports from the producer checker, without running its top-level writes. It validates the actual producer's 27,620 constructed fixed-path routes via explicit exception checks, and validates 35,857 generic directed-cycle simulations across all graphs through n=5, including every orientation of edges outside the chosen directed cycle. Replay checks intermediate directedness, inducedness, exact final reversal, and the r-2 move bound. This complements the independent BFS proof and checks the constructive code rather than only the existence theorem.

Both auditors pass normally and under python -O with byte-identical reports. No producer files were edited by the auditor; no remote writes were made.

## Literature check

The primary source https://arxiv.org/html/2108.12605v1 was checked directly. Theorem 3.1 proves the balanced-disagreement mechanism and the proof of Theorem 3.4 explicitly splits chorded directed cycles. The manuscript correctly credits these established mechanisms and carefully specializes the source's automorphism-based score-list setup to fixed labelled scores. No novelty or complexity inference is asserted.

Artifacts: audit.py, audit-result.json, audit-result-O.json; replay.py, replay-result.json, replay-result-O.json; corresponding normal and optimized stdout logs. These are independent audit artifacts, not replacements for the producer certificate.
