# Independent review of the constant-region circuit packet

**Verdict: PASS for the incidence identity, equality classification, and stated numerical deletion bounds. The minor endpoint scope clarification found during review is now resolved.**

This review read `circuits/PROOF.md` and `circuits/verify.py`, then independently derived the edge count and tested a different structural universe. It did not import producer functions into the independent checker.

## Mathematical review

For a closed constant region C, every boundary gate has precisely one distinct C parent. With internal indegrees 0, 1, and 2 counted by r, f, and |C|-r-f, respectively, the internal edge count is 2|C|-2r-f. Liveness and exclusion of the output from C give outdegree at least one for each C vertex. Therefore total outgoing edges are |C|+E, and subtracting the internal edges proves |C|+|B|=2r+f+E. Repeated pins are correctly collapsed to distinct-parent edges.

Equality requires f=E=0. Every C vertex then has exactly one outgoing gate edge, and every nonsource C vertex has two distinct C parents. A finite acyclic functional graph of this form has full binary merge trees as its weak components. Each component has one exit, and closure forbids two component roots from sharing a boundary gate. The proof characterizes the certified set, not necessarily every gate subsequently removed by semantic simplification.

The permissive-output deletion argument is correct in the full B2 basis: replace each C gate by its certified value, replace each B gate by a constant or signed predecessor, and absorb internal signs into surviving truth tables. All substitutions use ancestors, so no cycle is introduced. Under strict output accounting there is at most one terminal replacement gate; the D-1 bound is safe, and the nonconstant two-essential-input criterion correctly removes the terminal charge in the fixed-coordinate model.

## Endpoint scope issue and verified resolution

The hypotheses require output not in C, but do not require the restricted output to be semantically nonconstant when C is only a certified subset. For example, take A=a AND b, O=A AND x, fix a=0, and choose C={A}. Closure holds and O lies in B, yet the restricted output is the constant 0. If terminal constants are charged in the chosen strict model, the explanation should allow a terminal **constant-or-NOT gate**, rather than mentioning only a new terminal NOT. Alternatively, explicitly assume a nonconstant restricted output for that paragraph. The bound D-1 itself remains valid. This is a scope/wording clarification, not a counterexample to the theorem.

The author has now added an explicit semantic-nonconstancy hypothesis to Section 4 while keeping the incidence identity under its original weaker hypotheses. I reread that repair and verified it resolves the issue. No numerical bound or checker changed.

The independent checker includes this constant-output boundary control. The producer's semantic census deliberately filters constant outputs, so it does not test this particular case.

## Replays and independent finite checks

- The producer checker completed in normal and optimized Python. Outputs are byte-identical and match `circuits/receipt.json`.
- `independent.py` exhausts all live DAG presentations through four gates with two external vertices and **zero, one, or two input pins**, including repeated pins. This extends beyond the producer's exactly-two-pin topology universe.
- Presentation counts by gate number are 6, 60, 900, and 18,900. The independent check covers 4,776 closed constant regions, including 2,952 equality cases.
- Both normal and optimized independent outputs are byte-identical.
- All checks raise explicit exceptions; optimized Python does not remove them.

Artifacts: `normal.json`, `optimized.json`, `independent.py`, `independent.json`, and `independent.optimized.json` in this review directory.
