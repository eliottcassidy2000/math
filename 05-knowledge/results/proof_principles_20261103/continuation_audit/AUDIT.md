# Independent review of contextual one-probe invisibility

**Verdict: C1--C4 pass independent mathematical review and exact controls.**
The reviewer is the parent agent, separate from the continuation author.
The checker here imports no producer code and counts directed Hamilton
words by literal permutations, instead of the producer's endpoint subset DP.

For C1, a Hamilton word has a unique exterior skeleton and ordered cover
of the module. An arbitrary extra vertex touches zero, one, or two module
runs. Converse preserves unmarked covers, exchanges fixed-slot starts and
ends, and transposes the two-slot endpoint matrix. Under endpoint balance
its skew difference has zero row and column sums; the cut expression is
zero because its linear and symmetric-quadratic parts vanish separately.
The two adjacent runs are necessarily distinct, and empty exterior ends
are included. This covers orphan paths without relying on deletion being
a path-preserving operation. Vertex-transitivity proves the required
balance separately for each run count, not just at degree level.

For C4, mixed differences of the equal response polynomial recover the
same positive capacities. Equality of the weighted degrees then makes the
changed arcs a balanced integral circulation. Cycle decomposition keeps
all changed arcs inside old SCCs. Applying the same argument in reverse
forces the same ordered SCC partition. This is a necessary condition:
the capacities depend on the tournament, so a balanced flip measured in
the old weights need not keep the new weights or H fixed.

The independent finite universe is all 29 choices of a labelled quotient
of order1 through3 and its distinguished hole, substituting C3 or its
converse; all840 one-vertex attachment comparisons agree. Every one-step
response collision through order4 is also checked (nine unordered pairs).
This is a bounded control, not the basis of the general proofs.

Three hostile computations use direct permutations as well:

- The strong pair has H15, common first response31, and second responses
  81 and87. All32 first responses agree.
- The regular seven-vertex example has H171 and unequal endpoint counts
  24/25 at vertex0; its next responses are403/407.
- The old-weight-balanced triangle has three capacities4 but flipping it
  changes H11 to15 and changes the weight matrix.

The script reports427 explicit checks. Its source-stability check remains
active in optimized Python. Normal and optimized receipts agree. The
producer independently checks the larger universe described in its own
receipt; these universes are not conflated. Universal two-port faithfulness
and a continuation-closed compressed state remain open here.
