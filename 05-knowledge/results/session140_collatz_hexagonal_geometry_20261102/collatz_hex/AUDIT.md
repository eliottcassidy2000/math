# Independent review and exact reproduction

2026-11-02, Session140.

The parallel Hadwiger--Nelson lane (`recover_hadwiger_nelson_bridges`) independently reviewed `PROOF.md` and reported PASS:

- Reopened the primary Kühn PDF and checked Theorem 2(iii), printed page 4. The disconnected-graph proviso is load-bearing: a VAP-free drawing alone would not suffice. The proof retains the stronger pseudoforest/no-infinite-nested-cycles condition.
- Checked the proper interval-level drawing using finite BFS levels and separate component strips.
- Checked that the brick-wall vertical matching contracts to square-grid neighbours `(x+-1,y+-1)` in the displayed coordinates.
- Checked the subcubic branch-set median construction, including repeated attachment vertices and zero-length arms, and preservation of directed one-outdegree after routing.
- Independently checked the five-step affine predecessor formula, unit survival, disjoint same-depth fibres, and the sharp-clock hostile rooted at 7.

The verifier ran successfully under ordinary Python and `python3 -O`; the JSON outputs are byte-identical. Explicit exception checks remain active in optimized mode. The script independently checks both forward composition and exhaustive inverse-fibre construction in the stated finite universes, along with the parallel lane's two odd-predecessor constructions and the known-cycle hexagon routing.

Scope: the review checks the new argument and correct application of the cited classical theorem, not a new proof of Kühn's theorem or of the inherited 3-adic separation lemma. The latter is used only in the graph-rigidity discussion, with the explicit 3-unit restriction and forward-extension proof. None of the four main theorems depends on inherited rigidity. No Collatz convergence conclusion is made.
