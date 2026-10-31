# Independent symmetry-corollary proof review

**Verdict: PASS.** Reviewed 2026-10-31 by the concurrent
`audit_binary_bank_census` agent, separately from the cut-response producer.
This note records a proof review, not an additional exhaustive computation.

Given the previously independently audited canonical cut law with positive
symmetric weights, the following claims were checked:

1. For `n>=2`, cardinality invariance of `F(S)` makes all pair mixed second
   differences equal, hence all pair weights equal to a positive `w`.
   Equal singleton values then force equal outdegrees, so the tournament is
   regular. Conversely, regularity and uniform weights give outgoing cut
   size `s(n-s)/2`, hence `F(S)=H+w*s*(n-s)`.
2. For an arc-transitive tournament of order at least three, choose a vertex
   that is both an arc head and tail. Arc transitivity makes its vertex
   orbit contain every head and tail, hence all vertices. It is therefore
   vertex-transitive and regular. Two-path-cover equivariance makes all
   pair weights equal.
3. The order-two exception is necessary: the single arc is arc-transitive
   but not regular, and its two singleton cuts have different responses.
4. If ordered two-path covers allow one empty path, their empty-path cases
   contribute `2H`. The remaining `n(n-1)` ordered endpoint entries sum to
   `n(n-1)w`, so `E=2H+n(n-1)w`.

No issue was found. The exact Paley 7 and Paley 11 numerical controls have
their own producer normal/optimized/relocation receipts in
`../receipts/verify_examples.*.json`; the proof reviewer did not claim a
second independent enumeration of those examples.
