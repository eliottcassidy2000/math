# Independent audit: canonical weighted-cut response

**Verdict: main theorem, strong pair-deletion bound, SCC equality characterization and exact algorithms pass.** This is a structural identity, not an unrestricted attainment theorem or a novelty-priority claim.

## Independent computation

`independent_check.py` enumerates literal permutations and possible break locations to obtain H, endpoint counts and every W entry. It constructs each child tournament explicitly and counts its Hamiltonian paths using an independent adjacency DP. It imports only the implementation under audit for comparisons and optimization calls; it does not import the producer verifier or endpoint tables.

For all1,100 labeled parent tournaments of orders0–5, the independent checker verifies:

*1,100 literal endpoint/two-path-cover tables;
*33,867 actual child cut counts;
*5,405 endpoint-divergence identities;
*10,650 positive integral pair weights and pair-deletion lower bounds;
*83,506 mixed Boolean differences;
*1,100 exact **sets** of SCC-suffix minimizers, not just their cardinalities;
*6,446 constrained minimum values and returned physical cuts;
*1,100 minimum-positive-extension optimizations, including None at orders0/1;
*1,100 gain-content comparisons against gcd of all actual child-DP gains;
*5,461 arbitrary endpoint-comparison words for free-orientation block insertion.

The normal and optimized receipts are in `normal.json` and `optimized.json`. The producer checker is also replayed independently in `producer.normal.json` and `producer.optimized.json`; its extra random orders6–9 and malformed-input controls supplement, rather than replace, the universal proof.

## Proof checks

### Coefficient balance and canonical orientation

For an arc u→v, W_uv equals the number of Hamiltonian paths containing that arc, by concatenating and cutting. At a vertex v, total outgoing use minus incoming use is s_v−t_v. Hence

`L_v=s_v-t_v+sum_u W_uv=sum_(v→u)(W_vu+W_uv)`.

This exactly converts the deletion polynomial into a weighted outgoing cut in the original orientation. No prescribed-degree orientation theorem or feasibility assumption is needed. The implementation's `degrees[v]` and all-cuts recurrence agree with this identity: adding v to an existing cut changes its capacity by weighted outdegree(v) minus the sum of **all** pair weights to vertices already inside.

### Integrality

The two distinguished endpoints belong to different paths in the definition of W_uv. Reversing their mutual arc therefore preserves both W_uv and W_vu. The reversal changes H by W_vu−W_uv, which is even by Rédei oddness. Consequently their sum is even. The cited exposition was independently read: Irving–Omar, [Theorem5](https://arxiv.org/html/2412.10572v1), explicitly states and proves tournament Hamiltonian-path oddness. No unverified parity strengthening is needed.

### Strong positivity

For each Hamiltonian path P of T−{u,v}, insert the consecutive block uv or vu with its internal orientation free. The number of valid insertions is exactly the stated endpoint-comparison sum. Its telescoping decomposition has constant2 and nonnegative remainder. The resulting ordered two-path covers are distinct: concatenate across the marked gap, delete the two marked vertices, and recover P, position and order. Thus `W_uv+W_vu >= 2H(T−{u,v})`. It is important that this is an injection into all covers, not a claim that deleting the marked endpoints from every cover yields a Hamiltonian path.

### Small-order conventions

For n=0, H(empty)=1 and adding the sole vertex gives1; there are no pair weights. This must be discharged separately before the displayed deletion expansion, whose sums are empty. For n=1, every extension has H=1 and the usual expansion works. For n=2, W_uv=W_vu=1 and b_uv=1; for the parent0→1, the four cut responses in bitmask order are1,3,1,1. The implementation handles all these cases correctly.

### Equality and submodularity

All pair weights are strictly positive, so equality H_new=H occurs exactly when the cut has no outgoing original arc. Such a cut cannot split a strong component; in the totally ordered tournament condensation it is exactly a suffix union. The empty tournament has zero SCCs and one cut, agreeing with k+1. Strict submodularity on incomparable subsets follows from the negative pair coefficients; nesting gives equality rather than strictness.

### Algorithms and bit costs

Endpoint DPs and the cover convolution cost O(n²2^n) arithmetic operations and O(n2^n+n²) integer storage. The all-values recurrence costs O(n2^n) further operations. Each W_uv is at most(n−1)!: concatenate a cover into a permutation with the ordered consecutive block uv. Thus b_uv≤(n−1)!; all counts and capacities have O(n log n) bits. Arithmetic-operation bounds are not unit-cost bit-complexity claims.

For disjoint forced memberships, the extra capacity `1+sum b` exceeds every feasible original cut. Source-reachable vertices in the final residual network provide the correct constrained cut. The implementation checks feasibility and exact flow/capacity equality. Edmonds–Karp yields a polynomial number of arithmetic operations after compilation (O(n⁵) is a sufficient dense-network bound), with polynomial-bit capacities. An integral-value termination argument alone would only show termination, not this polynomial bound. Computing the input weights remains exponential; no polynomial-time Hamiltonian counting result follows.

## Follow-on small-gain and packing audit

The producer has explicitly discharged n=0 before the deletion polynomial and added the bit-size clarification. These editorial issues are resolved.

For positive gain2k, every crossing arc has `H(T−{u,v})≤b_uv≤k`, because all its weights contribute positively to a sum of k. For gains2 or4, positivity and oddness force all these pair-deletion counts to be1. A tournament has H=1 exactly when transitive: any directed triangle has three paths, and subsequent vertex additions cannot lower H by the already-proved response identity. Thus the parent has a feedback vertex set of size at most two. This is necessary, not sufficient, and does not concern zero-gain SCC suffixes.

A fixed packing of ν disjoint directed triangles contributes3^ν to the nonnegative odd-cycle collection formula. Deleting any two vertices leaves at least max(ν−2,0) triangles of that fixed packing. Therefore every positive extension gain is at least `2·3^max(ν−2,0)`. The cited Irving–Omar HTML was independently checked: Section3.4, Corollary19 expresses paths in the complement using literal odd cycles in T. Reversing paths identifies this count with H(T); complementary loops cannot occur in Hamiltonian vertex words. This resolves the complement convention without conflating different cycles with equal vertex supports.

`independent_consequences.py` additionally checks32 independently seeded random parents of orders6–9 and the ordinal sums of three/four triangles. It computes maximum triangle packing by an independent mask recursion, pair-deletion H values by actual induced adjacency DP, all weighted-cut gains, and selected actual child counts. Its receipts are `consequences.normal.json` and `consequences.optimized.json`; the initial line replays the full n≤5 independent audit. Producer consequence checks are also separately replayed.

The new `minimum_positive()` method is correct: every positive cut crosses at least one original arc, and forcing the endpoints of any arc to opposite sides guarantees a positive cut. Taking the minimum across these constrained minima therefore gives the true least positive extension count. For orders0/1 there are no arcs and the method returns None. Our independent all-parent audit compares this method to exhaustive cut values.

The exact gain-content identity also passes: writing G(S)=F(S)−H, all G(S) are multiples of twice the pair-weight gcd, while `G({u,v})−G({u})−G({v})=−2b_uv` proves the reverse divisibility. Both gcds are0 at orders0/1. The implementation is independently compared against the gcd of all actual child-DP gains, not merely against another calculation on the weights. This congruence constraint is a statement about extensions of a fixed parent; it is not automatically inherited by later descendants whose weights change.

No assertion of universal target-cut existence, closed recursive compiler domains, new globally forbidden H-values, or complete tournament spectrum is made.
