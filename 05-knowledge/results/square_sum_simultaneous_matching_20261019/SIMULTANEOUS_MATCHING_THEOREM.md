# An infinite family of simultaneous prescribed square-sum edges

## Theorem

Put Q=49 and T_k=(107·49^k−1)/2. For every integer k≥3 and every integer N≥T_k, there is a square-sum Hamiltonian cycle on [1,N] containing the consecutive seed482,3362,6242 and, simultaneously, every edge

    {a, 4·49^k−a}

for integers a satisfying

    (49^k+1)/2 ≤ a ≤ (3·49^k−1)/2,
    |a−49^k| > (49^(k−1)−1)/2.

These are 48·49^(k−1) pairwise vertex-disjoint prescribed edges. Their common sum is (2·7^k)². For any fixed finite collection of such levels, the same constructed cycle contains their union once N is at least the largest corresponding threshold. All these edges coexist with the prescribed seed.

This is a corollary of the frozen constructive transition theorem, proved below by retention of internal edges. No additional glue search or new base census is required. It is not a theorem about arbitrary prescribed edges or arbitrary connected prescribed paths, and no novelty claim is made.

## Definitions

Let L_k be the set of allowed lower endpoints a displayed above, and let E_k={{a,4Q^k−a}:a∈L_k}. The same definition is useful for k=1 and k=2. In particular,

    L_1 = {25,…,73} minus {49}.

The thresholds satisfy T_0=53 and T_k=49T_(k−1)+24. They guarantee one additional level of the certified construction. The eventual-cycle theorem uses B nice pairs, protected Aplus/Aminus nice pairs, and the final split lift. Its definitions and complete finite certificate tables are dependencies of this corollary.

## 1. The internal-edge lifting lemma

Suppose both paths of a nice input pair contain an undirected edge {a,b}, and none of the allowed cuts removes that edge. The common positional parity of a in the input pair fixes a sign σ∈{−1,1}, independently of which input member is selected at any offset c. Since a and b are adjacent, their signs are opposite.

For every offset −24≤c≤24, the selected lifted block therefore retains the edge

    {49a+σc, 49b−σc}.

A whole block or segment can be reversed without changing this undirected edge. Thus all 49 edges are present after gluing, irrespective of the word's orientations or its other connector edges. As c ranges across the 49 offsets, the first endpoint runs through the complete integer interval [49a−24,49a+24].

If a+b=4Q^k, all these new edge sums are4Q^(k+1). Distinct old lower endpoints produce disjoint consecutive 49-integer intervals. Therefore

    ⋃_(a∈L_k) [49a−24,49a+24] = L_(k+1).

To check the endpoints, write h_k=(Q^k−1)/2. The full allowed outer interval has endpoints (Q^k+1)/2 and(3Q^k−1)/2. Its lifted endpoints are those of the next level. The excluded central interval [Q^k−h_(k−1),Q^k+h_(k−1)] lifts to [Q^(k+1)−h_k,Q^(k+1)+h_k], since49h_(k−1)+24=h_k. This proves the set identity for all k, not just checked sample values.

## 2. Retention in B nice pairs

Claim: for every k≥1 and every n≥T_k, both paths of the constructed B_n contain every edge of E_k.

For k=1, every B input has prefix1,3. The only allowed B split is the prefix split at offset c=0 in the two residue48 states. Hence every nonzero-offset block retains its lifted1–3 edge:

    {49−c,147+c},  −24≤c≤24, c≠0.

These are exactly E_1. This argument applies for every n≥T_1=2621, because its input index is at least53.

For the induction step, write n=49m+d with24≤d≤72. If n≥T_(k+1), then m≥T_k. Both B_m paths contain E_k by induction. Every permitted split cuts only the input edge1–3, while all endpoints in E_k are at least25. Therefore no E_k edge is cut. Apply the internal-edge lifting lemma. It gives E_(k+1) in both output paths, independently of their selected inputs, block orientations, and prefix-split states. This proves the claim.

## 3. Retention in the protected A pairs

For k≥2 and n≥T_k, the B-to-A input index is at leastT_(k−1). Both input paths therefore contain E_(k−1). In either A family the only permitted cut is again the input prefix edge1–3, at one of offsets0 or−2. Since k−1≥1, all E_(k−1) endpoints are at least25. None is that prefix edge, so every copy is retained. The lifting lemma gives E_k in both A output paths.

The protected singleton edge2–7 introduced by the A glue does not remove any retained internal block edge. The known A prefix and endpoint conditions remain valid because the same verified transition words are used.

## 4. Retention in the final cycle

Fix k≥3 and N≥T_k. Its final input index n=floor((N−24)/49) is at leastT_(k−1). The selected protected A pair contains all E_(k−1) edges by the previous section.

The final lift cuts only the input edge2–7. Here k−1≥2, so the smallest endpoint of E_(k−1) is at least1201. Thus none of these edges is cut. Apply the internal-edge lifting lemma once more: the cycle contains every E_k edge. The same final certificates separately preserve the seed482–3362–6242. No new gluing requirements have been added.

For multiple levels, each retention argument applies to the same recursive construction, not a separately chosen cycle. At a common threshold all the levelwise conclusions therefore hold simultaneously.

## 5. Matching property and size

For each level, all lower endpoints lie below1.5Q^k while the corresponding upper endpoints exceed2.5Q^k. The map a↦4Q^k−a is injective, so the edges are pairwise vertex-disjoint. The outer lower-endpoint interval has Q^k integers, and the omitted central interval has Q^(k−1), leaving48Q^(k−1) edges.

Different levels also have disjoint vertex ranges: the highest endpoint at level k is less than3.5Q^k, while the lowest endpoint at level k+1 exceeds24.5Q^k. At levels k≥3 all endpoints exceed the seed vertices, so coexistence causes no degree conflict.

## 6. Finite checks and a large streaming control

`verify_matching_retention.py` checks the entire frozen294-state typed table for precisely the stated split exceptions and permissible atom types, checks all392 final states use the stipulated2–7 split model, verifies the interval identities at the first three levels, and records the general recurrence identities used in the proof. The universal transition and final-table replay remain the source of those tables' Hamiltonian validity.

A separate streaming control constructs a cycle of length6469494, without saving a second huge cycle. It deliberately uses the repaired B2694/2695 transition at parity0 residue48, then Aplus at index132030, and the final lift. It checks all115248 prescribed level3 matching edges of sum470596=686², both seed edges, exact interval coverage and every cyclic edge. It passes12938991 checks, in both normal and optimized modes. The level3 sufficient threshold is6294221.

The streaming control is evidence about one large instance. The all-k, all-N result follows from the retention induction above and the already complete universal construction, not from extrapolation of that instance.

## Reproduction and dependencies

Place this directory next to `../seeded_eventual_cycle_theorem_20261019/`, `../seeded_final_split_table_20261019/`, and `../square_sum_B_base_interval_20261019/` under the repository results directory. Run `python verify_matching_retention.py` and `python -O verify_matching_retention.py`. The optional large control is `python check_large_matching_control.py`; it streams rather than stores the large cycle. Both normal and optimized receipts are included. The retention checker hashes all consumed transition tables. The release manifest additionally pins the construction and original dependency manifests.
