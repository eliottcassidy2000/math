# Independent review: arbitrary gaps in the 1122 crossing skeleton

**PASS.** The four-parameter Hamilton-path formula and its exact eventual pure-ternary cancellation criterion are sound for the stated physical 1122 skeleton. The proof handles every zero-size boundary explicitly. It does not classify arbitrary crossing words or arbitrary unit clocks.

## Frozen source binding

Directory: `session197-prime-atoms/transitive-crossword-gaps/`.

- PROOF.md: `3027ec2684d7be7d1a497b1602b701c5f96c5f19ce7e550673c568fcdca415a3`
- verify.cpp: `9097be97e0fb8391c7c3b0291c7f0870084c259eacb54753ff95dd2225e4d6d4`
- receipt-normal.json: `82faea8249aa01f8465767dd2dfff4f22da354567be4158c3ccfc720993dc3b5`

The original files remain unchanged. Copied optimized and unoptimized library-debug builds reproduce the frozen receipt exactly.

## All-order proof audit

For each gap module, collapse all but a single representative before identifying its crossing graph. The earlier chain segment contains precisely the indicated anchors and any intervening zero modules. Later chain vertices and root-module vertices lose to the entire earlier segment. The outside vertex u has no outgoing chain edge. Hence the i and j crossing graphs are nonempty stars centered at a, and the k and root crossing graphs have cover {a,b} and two disjoint edges. Their ranks are uniformly one, one, two, two, including arbitrary zero sizes of the other modules.

The run-collapse/Stirling argument therefore puts each positive-coordinate restriction in the stated exponential span. Adding a delta function at zero supplies the independent deletion value. Successive univariate interpolation on the full 3-by-4-by-3-by-4 rectangle is legitimate: after extending one coordinate, every remaining coordinate retains a uniform bound independent of that extended value. It does not require a guessed total degree or finite stabilization.

The proposed formula belongs separately to the same four function spaces. In particular U_t=2^t-1+delta_t, so its zero boundary and the delta_j term fit the announced spaces exactly. Thus agreement on all 144 grid graphs proves the full formula on the nonnegative parameter lattice.

The cancellation criterion is exact. Both summands inside -D/2 are nonnegative. The first vanishes exactly at k=0 or 1, using the displayed recurrence for 3^k-2^(k+1)+1; the second vanishes exactly at i=0 or 1. All multiplicative factors are strictly positive. A is always positive. Since H(m+1)-3H(m)=-D*2^m-2C, an eventual constant affine offset is possible exactly when D=0. At those parameters U_i=U_k=1, so the offset is necessarily -1 at j=0 and +1 at j>=1. This proves the eventual claim as well as the exact claim.

The full run polynomial follows from the independently certified m=0 value, not from extrapolating a positive-clone formula. Leading zero padding and trailing root-preclone padding are genuine graph identifications and do not alter the cancellation criterion. No extension to a different nonzero skeleton, a type-3 letter, or a different outside orientation follows.

## Independent computation

A new checker builds the graph using the proof label order a,u,b. It evaluates every one of the 144 grid graphs by inclusion-exclusion over induced vertex sets of unrestricted walks of length n-1. This is independent of both source methods, endpoint-subset Hamilton DP and literal Hamilton-word DFS. Intermediate walk arithmetic uses signed 128-bit integers; final values are checked before conversion. All 144 counts agree, every positive-root graph is strongly connected, and the finite cancellation controls agree.

Optimized and unoptimized library-debug independent runs produce identical receipts.

- Independent checker: `c45caa362ecd116676b1b3ea1f9e25aedef175692f5f1923db0d5072d8b20222`
- Independent receipt: `5dd18295bac4219cb23d9b38f720d5ae53eab20b70f96bdf90d17e676f1aaab0`

These controls verify a complete tensor certificate whose sufficiency is proved above. No mathematical correction or unresolved blocker remains.
