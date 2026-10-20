# Walsh-support bijections: session 10 packet

**Current result:** for each fixed odd prime p, an explicit modular Walsh matching bijects tournament types and Royle-sign-even graph types while preserving the entire S_n-conjugacy class of a Sylow p-subgroup of their automorphism groups. It therefore preserves every p-group symmetry predicate simultaneously. This proves the prime-power portion of HYP-9340(a). The mixed-prime odd-group case remains OPEN.

## Reading order and scope precedence

1. `PRIME_VALUATION_THEOREM.md`: frozen first checkpoint; full weighted-matrix theorem, self-contained integral Walsh mechanism, deterministic O(N^3)-field-operation matching, and the first joint-prime obstruction at n=5.
2. `SYLOW_ACTION_ADDENDUM.md`: **supersedes the weaker limitation in the first checkpoint** that only subgroup order is preserved. The actual theorem preserves the Sylow subgroup's permutation action, and closes all odd p-group H symmetry counts. Its proof does not alter or invalidate the initial valuation proof.
3. `COKERNEL_ADDENDUM.md`: explicit odd-primary finite abelian group isomorphism; equivalent primewise Smith data, despite failure to match displayed cyclic summands by order.
4. `QUOTIENT_COROLLARY.md`: extension to every invariant binary reversal subspace U; switching classes versus sign-even Euler graphs are a special case. This is proved, not separately numerically enumerated here.
5. `CONVERSE_PARITY_ADDENDUM.md`: a further choice of matching preserves Sylow action and sends converse pairs to graph pairs of opposite edge parity; signed trace identity refined by Sylow action.
6. `SUPPORT_COUPLING_ADDENDUM.md`: canonical doubly stochastic coupling, total support, strict Hall expansion, arbitrary single row/column deletion, and an exact probability interpretation.

None is a relabelling-equivariant map on labelled objects, a polynomial-in-n algorithm, or a simultaneous-all-primes symmetry-preserving bijection. The word “even” always means orientation-sign-even unless Eulerian is separately stated. The earlier SCC–Lyndon rank codec remains a different construction with its own decomposition invariant; these maps use nonzero signed-copy coefficients instead of equal-rank pairing.

## Reproduction

The complete class universe is all 0<=n<=6. No inherited enumeration scripts or nauty are used. Python 3 and SymPy are required.

    python3 -B walsh_support.py --max-n 6 --out walsh_support.json
    python3 -B walsh_converse.py
    python3 -B walsh_sylow.py

The JSON files contain every matrix, matching, and the selected group-action witnesses. The `.out` files contain short exact run reports. `walsh_sylow.py` checks the signed localization congruence and actual relabelling witnesses, not merely automorphism orders. The `-O` reruns of the two addendum scripts match their standard output exactly; assertions critical to correctness use an explicit `need` function, so optimization does not disable them. The root also independently rebuilt and checked the initial packet.

`CHECKPOINT_SHA256.txt` pins the original four-file first checkpoint. `FINAL_SHA256.txt` pins the final research packet, excluding its own hash and redundant optimized transcripts. No repository edits, commits, or Git writes were performed by this worker.

## Proof dependencies and primary-source boundaries

The Boolean Walsh basis and signed invariant-orbit construction are inherited from THM-4721 and THM-4531 and are reproduced in the initial proof note. Royle et al.'s 2023 paper is the cited enumeration theorem and source of the natural-unlabelled-bijection question. The doubly stochastic vocabulary is classical; Sinkhorn–Knopp's primary paper is linked in the support note. The new research steps in this packet are the concrete modular matching, p-group signed localization, action-preserving consequence, explicit quotient extension, and converse/parity refinement. No global priority claim is made for these applications or for the underlying general matrix and p-orbit mechanisms.

## Inheritance, concept board, and stopping boundary

The closest inherited mechanism is the integral signed Walsh basis. The canonical hostile is the empty graph: it has the full symmetric automorphism group and cannot receive an arbitrary labelled equivariant tournament map. The corrected near miss is mistaking linear basis equivalence for a class set map; the present cofactor recursion explicitly supplies that set map. The least-used sidecar is odd-primary integrality of the inverse matrix, which forces the valuation filtration.

The live concepts were: signed-copy support; integral lattice pairing; p-group orbit localization; converse eigenspaces; naturality versus canonical class selection. Their connection is concrete: a p-unit signed coefficient forces a fixed graph copy for a tournament's p-Sylow action. The primary remaining question is mixed-prime subgroup symmetry; the first obstruction to a single all-primes map is already n=5. Improving the atlas-dependent construction to a local graph operation remains open.
