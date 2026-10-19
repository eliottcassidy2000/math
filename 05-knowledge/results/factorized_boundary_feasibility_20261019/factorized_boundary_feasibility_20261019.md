# A rank-one feasible boundary array for the second (3,3,1) obstruction

Packet 146, 2026-10-19. Root-reviewed exact relaxation counterexample.

## Exact result

For predecessor masks (0,28,0,0,12,12,13,125,125) and order (2,3,0,4,5,6,1,7,8), there is an exact nonnegative factorized boundary array whose probability for every forward incomparable pair strictly exceeds 2/3. Consequently, retaining rank-one factorization alone cannot exclude this core/order. This is a counterexample to the proposed stronger relaxation, not a counterexample to the poset theorem: actual exterior realizability and seven-thin exterior degree constraints are not asserted.

## Array

The only nonzero prefix increments are U_empty=1, U_(2)=1, U_(2,3,4,5)=6. The only nonzero suffix increments are V_empty=1, V_(8)=1, V_(5,6,7,8)=6. All unlisted increments are zero. Every generator coefficient is U_u V_v, including unsupported combinations: the product 36 of the two length-four words contributes zero because these two words both contain vertex 5 and cannot occur as disjoint prefix and suffix of an extension. Unsupported products are not set to zero as array entries.

For each core extension w, let h be the first position of roots 0 and 1 and k the position immediately after the later root. Give w weight

    (sum over i<=h of U_(w[:i])) (sum over j>=k of V_(w[j:])).

All weights are nonnegative, and each extension has weight at least 1. The empty normalization is strictly positive on each side.

There are 408 core extensions and 133 supported prefix-suffix generators. Eight products have nonzero supported contributions. Their total weighted mass is 948. For the forward pairs

    (0,1),(2,0),(3,0),(0,4),(0,5),(5,1),(6,1),
    (1,7),(1,8),(2,3),(4,5),(4,6),(5,6),(7,8),

the exact scaled margins 3 N_forward - 2 mass are

    705,300,84,3,246,246,3,84,300,12,3,228,3,12.

In particular, the minimum forward probability is 211/316 = 2/3 + 1/948.

## Consequences and next obstacle

No sound polynomial-multiplier certificate can exclude all nonnegative factorized arrays for this core/order, regardless of multiplier degree. Any exclusion must use additional necessary conditions for actual exteriors. A particularly concrete feature is that the only increments containing common vertex 5 are the two length-four words, each with value 6. The shared external incomparability allowance at that vertex is only 3, so deriving bounds that jointly constrain these two increments is a natural next step. This packet itself makes no such derivation or impossibility claim.

## Independent replay

Run replay.py in this directory, or invoke it by absolute path. Python normal and optimized (-O) execution produce identical receipts. The checker uses only Python's standard library and exact integers/Fractions. It enumerates all 9! permutations, tests predecessor masks directly, computes each boundary sum and each forward pair mass, reconstructs the support count, and checks strict positivity with explicit exception guards of every margin. It imports neither the scan engine nor an optimizer. It writes replay_report.json and weighted_extensions.json.

Optimizer artifacts are excluded from this packet. The numerical candidate had many tiny increments; direct rounding gave the much simpler exact integer witness above. Exact independent replay, rather than the optimizer, establishes the result.

The verifier also checks the exact displayed pair and margin lists and retains the unsupported product 36 while checking its extension count is zero. The witness is independently excluded by the additional actual-exterior ancestor-cap inequalities in packet 147; that does not change this packet's rank-one-relaxation conclusion.
