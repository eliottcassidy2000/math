# Independent review: Session 197 two-port interface and H49 supplement

2026-11-05. **PASS**, after the minor empty-product order clarification and two explicit certificate-coverage gates described below. This review concerns the exact stated construction language, not global spectrum coverage.

## Main mathematical audit

- The 27 signatures retain the endpoint identity of each component through its unique start label. They contain exactly the information needed to choose the new vertex's predecessor and successor components and to exclude joining both ends of one component into a cycle.
- In the inverse restriction, deleting the newborn vertex can create a new component start only at an old port, because the newborn loses to every old non-port. Restoring a forgotten old port can only reduce the number of non-port starts. The old at-most-one condition therefore follows from the counted output condition. This closes the potentially delicate forgetting step.
- The raw matrices are nonnegative. Coordinate capping commutes with every matrix, including operations that later discard a large coordinate's contribution. The H-readout and induced-path injection justify the scalar pruning. The negative reduced 17-state matrices are correctly not used for saturation.
- The 17-dimensional upper certificate is checked on the full 27-coordinate observation matrix, not merely a finite collection of reachable inputs. Seventeen physical prefix/suffix words, each of length at most three, give the independently counted nonsingular continuation matrix with determinant -2^21. The empty suffix is present.
- The H67 and H173 words produce the displayed actual graphs. Their counts are checked by independent adjacency subset-DP, literal permutation enumeration, and full cover-state propagation.
- The stationary full-cover clone proof correctly bounds the total number of clone runs by three: each extra run consumes an old port not already consumed as a component start. This applies even when some components contain no clone.
- The ordered-join reset really replaces both exposed old ports by a fresh transitive seed and sends every old-to-future arc forward. The 23-smooth factor compiler handles all odd 7-exponent residues, with precisely 7 and 21 excluded from the stated domain.

One minor statement was corrected: exact order=sum of seed orders≤12t is stated for N>1,t≥1. The empty product N=1 uses the two-vertex starting seed.

## Reproduction performed independently

- Full `verify.py --bound 21` normally and under `-O`; outputs agree byte for byte.
- 820 input profiles and 9840 physical birth checks passed; all selected witnesses, 16 smooth seeds, the rank certificate and four ordered-product controls passed.
- `emit_automaton.py --check` confirmed both serialized automata agree with the direct combinatorial generator.
- Recompiled `census.cpp` in a separate temporary executable and reproduced the complete ceiling 173 closure: 485591 states, maximum shortest word 19, exactly 85 odd scalar values and missing 7,21.

The main proof reviewed has SHA256 `0883968b2441bc349412a4727d672939abef90b394051fc5cbbdbbe9220e51d9`. The transition source has SHA256 `733155121ac7a5936a38a0ec1f3a602759f868c1e4e547fd4b6f89bced25f9cb`; the rank certificate has SHA256 `b73807dedcb4c82157c86bdbc1ca67f3ffe12fc1930e96ba62ae72256a3a6096`.

## H49 supplement

The exact named S196 global-induction source was read at commit `cf3b44b9a442c3cff0f28cab924e73c4c00bb00c`. It does state the physical H49 input and its old arbitrary-gadget descendant support `{49,85,147}` through 151.

The new two-port supplement correctly permits forgetting root 5 later. It starts from independently counted physical profiles for all six choices of second port; there is no assumption that this initial tournament is itself reachable from the main seed.

A separately compiled C++ replay produced all six entire sorted state sets, not merely their cardinalities. They match the Python certificates. Python normal and optimized verifier runs agree. The case sizes/supports are:

- p0: 335 states, `{49,85,147}`.
- p1,p3,p4: 347 states each, `{49,85,147}`.
- p2: 486 states, `{49,85,115,147}`.
- p6: 485 states, `{49,85,91,147}`.

All have maximum shortest capped depth 10. The displayed two-seam H115 graph is independently counted by both graph DP and permutations. The exclusion of 149 and 151 follows from closure and reachability, not a depth cutoff.

The certificate verifier was hardened to require exactly the six distinct second-port keys and to require the positive-witness keys equal the computed support in each case. The original current data already satisfied both properties; these gates prevent silent weakening after future edits.

The supplement proof has SHA256 `759fd3ce52820d0ea8b5d060e90ce250071d8eb4cedc606f75d8fb5004883586`, and its full state certificate has SHA256 `5b9ab1324db8dd6262a80c35f15a982ff1bfc7b519d8e57de0a547be1a19da4d`.

No remaining mathematical blocker was found. All repository reads were exact assigned/named paths. No worktree was created or removed, and no unrelated checkout was changed.
