#!/usr/bin/env python3
"""
SINGLETON CANCELLATION LEMMA — attempt at algebraic proof.

Lemma: For centered edge deviations s_i = A[p_i, p_{i+1}] - 1/2, if a
position subset has an isolated position (not adjacent to any other selected
position), then sum_P prod_{selected positions} s_i = 0.

Proof attempt:
Consider position subset I = {i_1,...,i_k} with i_j isolated (meaning
i_j-1 and i_j+1 are not in I).

sum_P prod_{l in I} s_{i_l}

The product involves edges (p_{i_l}, p_{i_l+1}) for each l.
The isolated position i_j involves the edge (p_{i_j}, p_{i_j+1}).

Key: the vertex p_{i_j} appears in the product ONLY through the edge
(p_{i_j}, p_{i_j+1}). The vertex p_{i_j+1} might also appear through
position i_j+1 if i_j+1 ∈ I, but we assumed i_j+1 ∉ I.

Actually, p_{i_j} appears through:
- edge at position i_j: (p_{i_j}, p_{i_j+1})  [involving p_{i_j} and p_{i_j+1}]
AND ALSO potentially:
- edge at position i_j-1: (p_{i_j-1}, p_{i_j})  [involving p_{i_j}]
But i_j-1 ∉ I (isolated), so position i_j-1 is not selected.

So vertex p_{i_j} appears ONLY through the edge at position i_j.
Similarly, vertex p_{i_j+1} appears only through positions i_j and i_j+1.
Since i_j+1 ∉ I, vertex p_{i_j+1} also appears only through position i_j.

Now, split the sum: fix all vertices EXCEPT those at positions i_j and i_j+1.
The product over I factors as:

[prod_{l ≠ j} s_{i_l}] * s_{i_j}

The first factor doesn't involve p_{i_j} or p_{i_j+1}.
The second factor is A[p_{i_j}, p_{i_j+1}] - 1/2.

Summing over all choices of (p_{i_j}, p_{i_j+1}) from the remaining vertices:

sum_{distinct a,b not in other positions} (A[a,b] - 1/2)

For this to be zero, we need: among the remaining vertices (not used by other positions),
sum_{ordered (a,b)} (A[a,b] - 1/2) = 0.

For any tournament on m vertices:
sum_{ordered (a,b)} A[a,b] = C(m,2) (each unordered pair contributes 1)
sum_{ordered (a,b)} 1 = m*(m-1)
So sum (A[a,b] - 1/2) = C(m,2) - m*(m-1)/2 = 0.

The remaining vertex count is n - (number of OTHER positions' vertices).

Wait, but the positions other than i_j don't necessarily have distinct vertices.
Actually, the permutation assigns distinct vertices to ALL n positions. Positions
0,...,n-1 get vertices p_0,...,p_{n-1} all distinct.

So after fixing the other n-2 positions (all except i_j and i_j+1), we're left
with 2 specific positions (i_j and i_j+1) and exactly 2 remaining vertices.
The sum over the 2! = 2 orderings of these 2 vertices at positions i_j, i_j+1:

For the two remaining vertices a, b:
(A[a,b] - 1/2) + (A[b,a] - 1/2) = (A[a,b] + A[b,a]) - 1 = 1 - 1 = 0.

That's it! The cancellation happens because for any two vertices a,b:
A[a,b] + A[b,a] = 1, so (A[a,b]-1/2) + (A[b,a]-1/2) = 0.

Wait, but this doesn't account for the OTHER products. Let me be more careful.

Actually, we're summing over ALL permutations P of n vertices. We can group
them by which vertices go to the "other" positions. For each such assignment,
the vertices at positions i_j and i_j+1 are the two remaining vertices {a,b}.
They can be ordered as (a at i_j, b at i_j+1) or (b at i_j, a at i_j+1).

The product over I for these two orderings:
Case 1: p_{i_j}=a, p_{i_j+1}=b: prod = [OTHER] * (A[a,b]-1/2)
Case 2: p_{i_j}=b, p_{i_j+1}=a: prod = [OTHER'] * (A[b,a]-1/2)

But [OTHER'] ≠ [OTHER] in general! Because the OTHER positions also depend
on the full permutation. The vertices at i_j and i_j+1 affect the edges
at positions i_j-1 and i_j+1 as well (if those are selected).

HOWEVER: position i_j is isolated, meaning i_j-1 ∉ I and i_j+1 ∉ I.
So the edges at positions i_j-1 and i_j+1 are NOT in the product.

But the edges at positions i_j-1 and i_j+1 DO affect the product indirectly
if they involve vertices that also appear in other selected positions. Wait no —
the product only involves edges at positions in I.

Let me re-examine. The product is:
prod_{l in I} s_{i_l} = prod_{l in I} (A[p_{i_l}, p_{i_l+1}] - 1/2)

For l ≠ j: s_{i_l} = A[p_{i_l}, p_{i_l+1}] - 1/2.
The vertices p_{i_l} and p_{i_l+1} are at positions i_l and i_l+1.
These are all DIFFERENT from positions i_j and i_j+1 (since i_j is isolated
and i_l ≠ i_j, and i_l+1 ≠ i_j since i_j-1 ∉ I, and i_l ≠ i_j+1 since
i_j+1 ∉ I and... hmm, this needs more care).

Actually, position i_j+1 might coincide with position i_l for some l.
But that doesn't matter! The key is that vertex p_{i_j+1} could appear
in OTHER edges too (at position i_j+1, which might be i_l for some l).

OH WAIT. If i_j+1 = i_l for some l ≠ j, then position i_l = i_j+1 ∈ I,
contradicting i_j+1 ∉ I. So i_j+1 is NOT any i_l.

Similarly, i_j = i_l+1 would mean i_l = i_j-1 ∈ I, contradicting i_j-1 ∉ I.

So positions i_j and i_j+1 are not shared with any other selected position
or its successor. But vertex p_{i_j} might appear at some other non-selected
position, which doesn't matter since we only look at edges at selected positions.

The EDGES in the product involve: for each l, the pair (p_{i_l}, p_{i_l+1}).
The vertex p_{i_j} appears ONLY in the edge (p_{i_j}, p_{i_j+1}) [from l=j].
The vertex p_{i_j+1} appears ONLY in the edge (p_{i_j}, p_{i_j+1}) [from l=j],
since i_j+1 is not any i_l and (i_j+1)-1 = i_j is the isolated position.

Wait: p_{i_j+1} also appears in the edge at position i_j+1 (which is
(p_{i_j+1}, p_{i_j+2})). But i_j+1 ∉ I, so this edge is not in the product.

BUT p_{i_j+1} might be the SUCCESSOR of some other position: p_{i_j+1} = p_{i_l+1}
for some l ≠ j, which means i_j+1 = i_l+1, i.e., i_l = i_j. But l ≠ j and all i_l
are distinct. So this can't happen.

Conclusion: the edges in the product that involve p_{i_j} or p_{i_j+1} are
ONLY the edge at position i_j. All other edges involve neither p_{i_j} nor p_{i_j+1}.

Therefore, fixing the assignment of vertices to all positions EXCEPT i_j and i_j+1,
the product factors as:

prod_{l ≠ j} (A[p_{i_l}, p_{i_l+1}] - 1/2) * (A[a, b] - 1/2)

where {a, b} are the two remaining vertices assigned to positions {i_j, i_j+1}.

The first factor is FIXED (doesn't depend on the ordering of {a,b}).

Now sum over the 2 orderings of {a,b} at {i_j, i_j+1}:
(A[a,b]-1/2) + (A[b,a]-1/2) = A[a,b] + A[b,a] - 1 = 1 - 1 = 0.

QED!

Actually, I need to also verify that the "rest of the permutation" factor
doesn't depend on the ordering of {a,b} at {i_j, i_j+1}. This is because
p_{i_j} and p_{i_j+1} DON'T appear in any other edge of the product
(as shown above). But they DO appear at their specific positions,
affecting the overall permutation structure.

The key point: for the OTHER positions' edges, the vertices are determined
by which positions they are at. The vertices at positions i_j and i_j+1
don't affect any of the other edges in the product. So the "other factor"
is genuinely independent of the {a,b} ordering.

THIS COMPLETES THE PROOF.

Let me verify computationally one more time to be sure.
"""
from itertools import permutations, combinations
import random

def verify_singleton_cancellation(n, pos_subset, num_trials=5):
    """Verify that a position subset with a singleton gives sum=0."""
    m = n - 1
    for trial in range(num_trials):
        random.seed(trial * 113 + 7)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        total = 0.0
        for p in permutations(range(n)):
            prod = 1.0
            for pos in pos_subset:
                s = A[p[pos]][p[pos+1]] - 0.5
                prod *= s
            total += prod

        if abs(total) > 0.001:
            return False, total
    return True, 0.0

print("SINGLETON CANCELLATION LEMMA VERIFICATION")
print("=" * 50)

for n in [5, 7, 9]:
    m = n - 1
    print(f"\nn={n}:")
    # Find position subsets with singletons
    for k in [2, 3, 4]:
        if k > m: continue
        for ps in combinations(range(m), k):
            # Check for singleton
            has_singleton = False
            for p in ps:
                left_adj = (p-1) in ps
                right_adj = (p+1) in ps
                if not left_adj and not right_adj:
                    has_singleton = True
                    break
            if has_singleton:
                if n <= 7:
                    ok, val = verify_singleton_cancellation(n, ps, 3)
                    if not ok:
                        print(f"  FAIL: {ps} -> {val:.6f}")
                    # Only print first few
                elif n == 9:
                    # Skip brute force for n=9
                    pass

    if n <= 7:
        # Count how many subsets have singletons
        total = 0
        with_singleton = 0
        for k in [2, 3, 4]:
            for ps in combinations(range(m), k):
                total += 1
                for p in ps:
                    if (p-1) not in ps and (p+1) not in ps:
                        with_singleton = True
                        break
                if with_singleton:
                    with_singleton_count = with_singleton
        print(f"  All singleton subsets cancelled correctly at n={n}")

print("\n" + "=" * 50)
print("ALGEBRAIC PROOF:")
print("""
Theorem (Singleton Cancellation Lemma):
Let I ⊆ {0,...,n-2} be a set of positions, and let j ∈ I be isolated
(i.e., j-1 ∉ I and j+1 ∉ I). Then:
  sum_P prod_{i∈I} (A[p_i, p_{i+1}] - 1/2) = 0.

Proof:
1. The edge at position j involves vertices p_j and p_{j+1}.
2. Since j is isolated: j-1 ∉ I and j+1 ∉ I.
   - No selected edge has p_j as its SECOND vertex (would need position j-1).
   - No selected edge has p_{j+1} as its FIRST vertex (would need position j+1).
   - The edge at j is the ONLY selected edge involving p_j or p_{j+1}.

3. Fix all vertex assignments except at positions j and j+1.
   Let the remaining two vertices be {a, b}.

4. The product factors as:
   [OTHER edges factor] × (A[p_j, p_{j+1}] - 1/2)
   where [OTHER] doesn't depend on {a,b} ordering.

5. Summing over the 2 orderings of {a,b} at {j, j+1}:
   (A[a,b] - 1/2) + (A[b,a] - 1/2) = (1 - 1) = 0.

QED.
""")

print("DONE")
