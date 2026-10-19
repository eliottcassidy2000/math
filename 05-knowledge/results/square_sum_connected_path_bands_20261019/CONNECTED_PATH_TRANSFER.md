# Connected prescribed paths on explicit geometric bands

## General transfer lemma

Let Q=49 and h_k=(Q^k−1)/2. Fix a certified B base pair B_b, 53≤b≤2620. Suppose the same undirected path P=(p_0,…,p_(r−1)), with r≥2 distinct vertices, occurs consecutively in both members, allowing the whole path to be reversed. Suppose P does not use the edge {1,3}.

For every k≥2 and every integer

    Q^k b + h_k ≤ N ≤ Q^k(b+1) + h_k − 1,

the certified seeded-cycle construction contains simultaneously all Q^k paths

    P_(k,t) = (Q^k p_j + (−1)^j t)_(j=0,…,r−1),
    −h_k ≤ t ≤ h_k.

These paths are pairwise vertex-disjoint. Every such cycle also contains the consecutive seed 482,3362,6242. This is a conditional transfer lemma whose hypothesis can be checked on one finite base pair; it is not a theorem for arbitrary prescribed paths.

### Exact bands and recursion depth

A construction step descends by f(n)=floor((n−24)/49). Induction gives

    f^k(N)=floor((N−h_k)/Q^k),

because h_(k+1)=49h_k+24. Hence the displayed band is exactly f^k(N)=b. Reading the construction upward, the first k−2 steps are B-to-B, the next is B-to-A, and the final one is A-to-cycle. Since b≥53, every B index after the first lift is at least2621; since b≤2620, the recursion stops at exactly the certified base b. The least N is at least Q^k·53+h_k=T_k, so the eventual construction applies.

### Edge retention, signs, and orientations

At the first lift the only possible cut is input edge1–3, which P avoids. At later B and A steps the only possible cuts are again1–3. Every once-lifted vertex is at least25, so these later cuts cannot touch a lifted copy. At the final step the input has already been lifted at least once; its vertices are at least25, so the final2–7 cut cannot touch a copy either. Reversal of a whole block or retained segment reverses traversal only, preserving each undirected copy of P.

Each pair shares the positional parity of every common vertex. Thus the selected input member at an offset changes neither the labels nor the alternating parity pattern of any retained path copy. For a copy with parameter t after l lifts, another lift sends its jth vertex to

    Q^(l+1) p_j + (−1)^j (49t + σ_t c),

where σ_t∈{−1,1} is fixed for that copy by its common starting-vertex parity, and −24≤c≤24. The sign may depend on t, but cannot affect the set of generated parameters: for each t it fills [49t−24,49t+24]. These disjoint consecutive intervals, over all t∈[−h_l,h_l], partition [−h_(l+1),h_(l+1)]. At the first lift the common initial parity merely replaces t by−t. This proves that all centered parameters occur, independently of orientations and parity signs, in both pair outputs at every nonfinal stage and in the final cycle.

### Disjointness and square sums

For any fixed p_j, the possible labels form exactly the interval [Q^k p_j−h_k,Q^k p_j+h_k], of length Q^k. Distinct integer p_j have disjoint such intervals. Within each interval the map t↦Q^k p_j+(−1)^j t is injective. Thus no vertex appears in two different copies, or twice in one copy. Each adjacent sum is Q^k(p_j+p_(j+1)), a square since Q^k=7^(2k). Positivity and containment in [1,N] also follow from the certified expansion; directly, the largest possible label is at most Q^k(b+1)+h_k, but actual P vertices belong to the shorter base path and hence are at most b, giving the sharper bound Q^k b+h_k≤N.

## Certified twelve-vertex instance

The frozen base pair B54/55 has the common prefix

    1,3,46,35,29,52,48,33,16,20,44,37,12.

Consequently take

    P=(3,46,35,29,52,48,33,16,20,44,37,12).

Its adjacent sums are

    49,81,64,81,100,81,49,36,64,81,49.

It avoids1–3. For every k≥2 and every N in

    [54·49^k+(49^k−1)/2, 55·49^k+(49^k−1)/2−1],

the cycle therefore contains 49^k simultaneous disjoint prescribed twelve-vertex paths P_(k,t), as well as the original seed. For k=2 this is every N from130854 through133254, with2401 prescribed paths and28812 prescribed vertices in each cycle.

The bands are pairwise disjoint and have large gaps; this statement does not say that this path family occurs for all sufficiently large N. For each k the prescribed labels are fixed throughout that band's N range. Increasing k changes the prescribed paths.

## Verification scope

The accompanying base-and-expansion checker reads the frozen B54 pair, verifies its common prefix and parity, then streams one cycle at the lower endpoint130854. It checks exact interval coverage, every cyclic square edge, the original seed, and every edge of all2401 twelve-vertex paths. The all-N and all-k conclusions follow from the transfer proof and complete frozen transition tables, not from that single expansion.

## Reproduction

Place this folder beside `../seeded_eventual_cycle_theorem_20261019/`, `../seeded_final_split_table_20261019/`, and `../square_sum_B_base_interval_20261019/`. Run `python verify_connected_path_transfer.py` and `python -O verify_connected_path_transfer.py`. Both receipts are included; no cycle file is generated. The sibling eventual theorem supplies the complete transition proofs. This release adds the transfer induction and finite base-path check. No historical novelty claim is made.
