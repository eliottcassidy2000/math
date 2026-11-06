# Independent audit of finite strategy5 certificate and finite-menu criterion

Status: PASS in their explicitly finite/conditional scopes. Neither is a threshold proof.

## Literal finite strategy certificate

Reviewed `STRATEGY5.tsv`, `STRATEGY5_SEEDS.txt`, `verify_strategy.cpp`, and `STRATEGY5_LITERAL.json` in `session204-interval-induction`.

The literal tournament convention is consistent: apex 0 beats apex 1, chain vertices are transitively ordered by word position, and each letter is the bit mask of apices beaten by that chain vertex. Every unordered vertex pair has exactly one directed edge. A map digit 4 means erasure, not a fifth incidence type. Rewriting rebuilds the entire chain using the image word in the same position order. End contexts and the one indexed interior insertion are then applied to that image word. These are actual new tournaments, without a claim that the original graph remains induced.

The verifier reconstructs every rewrite literally, directly counts both endpoints using Hamilton-path subset DP, checks the exact 2H-1 or 2H+1 law, and forms a directed graph with a separate child set for each sign. Its readiness predicate has the correct quantifiers: for each sign there exists an eligible child that is recursively ready at the remaining depth. It does not demand a single universal rewrite, and does not confuse separate existential trajectories with a common five-level both-branch strategy. The selected child may depend on the current word and remaining depth, which is valid for the stated finite strategy.

An independent rerun exactly reproduces the author receipt: 744 edges, 755 literal graphs, maximum order 12, twelve seeds with counts 23,25,...,45, and depth-five readiness for every seed. All odd targets 23 through 1471 occur among the recounted certificate vertices. An additional separately written reverse-path subset recurrence independently recounted all 755 graphs and checked all listed count values; PASS. Its source is `independent.cpp` alongside this review.

The supplied C++ verifier uses signed 64-bit counts; the actual certified order bound 12 makes this safe since even 12! is below the signed limit. Its generic n<22 guard alone would not guarantee overflow safety for arbitrary replacement certificates, so no broader such guarantee is inferred.

This certificate proves only the displayed finite strategy readiness and finite coverage. It gives no all-order closure, no invariant infinite witness language, and no new scalar range beyond what the project already knew. The mechanism tested is genuine target-dependent physical rebuilding with erasures and insertions.

## Finite-menu bridge and tensor-product criterion

Reviewed `FINITE_MENU_BRIDGE.md`. Its sufficient conditions are correct and explicitly unsatisfied.

For a fixed morphic/context rewrite p phi(w) q, chronological matrix products give the displayed 57-dimensional difference series: two 28-coordinate forest states and one constant coordinate. Erasure correctly uses the identity matrix. Exact reachable spans indexed by guard-DFA states prove the identity on every accepted word, not merely finite samples. Guard coverage and morphic-image containment are regular-language checks.

Saturating nonnegative integer forest coordinates at 23 commutes with all transfers and the clipped H readout, so the language H(w)>=23 is regular. For each sign and a fixed finite menu, the pointwise product of the difference series vanishes exactly when at least one menu branch succeeds. Tensor-product representations implement this pointwise product, and guard-restricted reachable-span zero testing is finite. No regular selector is required. Products of sums of squared defects also express the stated finite conjunction/existential conditions because all evaluations are real integers.

If both sign menus passed these all-word tests and literal seeds existed for 23..45, induction on the unique odd parent would realize every odd target at least 23. Successful outputs stay eligible because 2H-1>=45. This logical implication is valid; the note correctly says no actual menu has passed the required test.

A freely chosen interior position is not one fixed context/morphism rewrite. The note explicitly excludes silently treating the finite strategy's position choices as a fixed finite menu. A separate guarded or transducer encoding would be needed before applying the tensor criterion to such choices.

No required mathematical corrections. Keep finite strategy certification and the unsatisfied all-word criterion clearly separated.
