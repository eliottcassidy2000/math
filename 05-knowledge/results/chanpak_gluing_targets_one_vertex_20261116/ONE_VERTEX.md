# One shared vertex never helps; a corrected overlap bound

Notation:
- **Witnesses.** A and B are induced subtournaments of T, with W = A ∩ B.
- **Counts.** p = H(A), q = H(B), r = H(W), a = H(A − W), b = H(B − W), and H(∅) = 1.
- **Packs.** The pack expansion (OCF, THM-002) writes H as a sum over sets of vertex-disjoint directed
  odd cycles, a set of k cycles weighing 2^k.
- **Mixed weight M.** Packs inside A or inside B weigh p + q − r. Write M ≥ 0 for the weight of all other
  packs of T[A ∪ B] (THM-4114: every mixed Boolean difference is a nonnegative pack sum). Then
  H(A ∪ B) = p + q − r + M.

## Theorem 1 (PROVED, all orders)

Let p, q be odd primes, not both 3. If A and B share at most one vertex, every tournament containing
them as induced subtournaments has H ≥ p + q + 7.

This sharpens the bound p + q + 3 of
[CONDITIONAL.md](../conditional_determinant_vitali_20261116/CONDITIONAL.md). Consequently **none** of
its three rows (p + q − 1, p + q + 1, p + q + 3) can be implemented with induced witnesses meeting in at
most one vertex. The single exception is 3 + 3 + 3 = 9: two triangles through one vertex, H = 9.

*Proof.* By monotonicity (THM-1370 (2)) it suffices to bound H(A ∪ B). Each witness has a prime count,
so it has exactly one nontrivial strong component: K_A and K_B, with H(K_A) = p and H(K_B) = q.

- **Case 1: v is missing from a core.** Then K_A ∩ K_B = ∅. Packs of K_A and packs of K_B unite
  injectively, so H ≥ pq. Since (p − 1)(q − 1) ≥ 8 for odd primes not both 3, pq ≥ p + q + 7.
- **Case 2: v lies in both cores.** A strong tournament has a directed triangle through every vertex
  (Moon), so pick triangles t_A ∋ v in K_A and t_B ∋ v in K_B. Their five-vertex union is strong, so
  its count is at least 9. Its only odd cycles inside one side are t_A and t_B. Hence it contains two
  mixed odd cycles X and Y, which contribute 4 to M (the bowtie of CONDITIONAL.md).
  - **2a: K_A has an odd cycle C avoiding v.** C misses B, so {C, t_B} is a mixed pack of weight 4,
    distinct from {X} and {Y}. Then M ≥ 8. The case where K_B has such a cycle is the same.
  - **2b: every odd cycle of both cores passes through v.** Then K_A − v is a transitive tape
    a_1 → … → a_k and K_B − v is a tape b_1 → … → b_l. Strongness forces v → a_1, a_k → v, v → b_1 and
    b_l → v. We exhibit four distinct mixed odd cycles, so M ≥ 8.
    - **Y1, from the arc between a_k and b_1.** If b_1 → a_k, take the triangle v b_1 a_k. If
      a_k → b_1, take the 5-cycle v → a_1 → a_k → b_1 → b_l → v.
    - **Y2, from the arc between b_l and a_1.** If a_1 → b_l, take the triangle v a_1 b_l. If
      b_l → a_1, take the 5-cycle v → b_1 → b_l → a_1 → a_k → v.
    - Y1 ≠ Y2: either their vertex sets differ, or the successor of v differs.
    - **If p ≥ 5, then k ≥ 3: a third cycle Y3 through a middle vertex a_j (1 < j < k).**
      - If v → a_j: take the triangle v a_j b_l when a_j → b_l, otherwise
        v → b_1 → b_l → a_j → a_k → v.
      - If a_j → v: take the triangle v b_1 a_j when b_1 → a_j, otherwise
        v → a_1 → a_j → b_1 → b_l → v.
      - Y3 contains a_j; Y1 and Y2 do not.
    - **If also q ≥ 5: a fourth cycle Y4 through a middle b_i (1 < i < l).**
      - If v → b_i: take the triangle v b_i a_k when b_i → a_k, otherwise
        v → a_1 → a_k → b_i → b_l → v.
      - If b_i → v: take the triangle v a_1 b_i when a_1 → b_i, otherwise
        v → b_1 → b_i → a_1 → a_k → v.
      - Y1–Y3 use only b_1 and b_l, so Y4 is distinct.
    - **If q = 3:** K_B is the triangle v → b_1 → b_2 → v. Pass to the converse if necessary so that
      v → a_j for some middle j; the converse preserves counts, odd cycles, both tape structures and the
      triangle. Then there is a fourth cycle Y5:
      - if a_k → b_1: v → a_j → a_k → b_1 → b_2 → v;
      - if b_1 → a_k and a_j → b_1: v → a_1 → a_j → b_1 → b_2 → v;
      - if b_1 → a_k, b_1 → a_j and a_1 → b_1: v → a_1 → b_1 → a_j → a_k → v;
      - if b_1 → a_k, b_1 → a_j and b_1 → a_1: v → b_1 → a_1 → a_j → a_k → v.

      Each Y5 differs from Y1, Y2 and Y3 by its vertex set or by v's successor.
    - The case p = 3 ≤ q is symmetric.

So M ≥ 8 and H ≥ p + q − 1 + 8. ∎

**Census (FINITE-EXACT).** `jglue.c` checks this over every cover with |W| = 1 at orders 6–9:
- for prime pairs not both 3, the minimum of H − p − q is exactly 7, attained only by (3, 5);
- for (3, 3) the minimum is 3.

`onevx.c` lists the least one-vertex excess for each pair (`onevx.out`):
- (3, q) and (5, q) attain H = pq within the census range;
- (11, 11) needs H ≥ 135 at orders ≤ 9; its disjoint-core product 121 needs order 10;
- (13, 13) needs H ≥ 187.

## HYP-9862: one-vertex supermultiplicativity (OPEN)

**Conjecture.** For induced A and B with |A ∩ B| ≤ 1, H(T[A ∪ B]) ≥ H(A)·H(B).

**Evidence (`supmult.c`, `supmult.out`).** Exhaustive over every tournament of orders 4–9 and every
cover with exactly one shared vertex. At order 9 this is 218,925,648 instances; there are no
violations, and 139,658 equalities.

**Sharp at two shared vertices.** The conjecture fails as soon as two vertices are shared. Two triangles
on an edge give H = 5 < 9, and at order 9 there are 695,729 two-vertex failures among 434,403,648
instances.
- Disjointness (|A ∩ B| = 0) is the classical product bound.
- The conjecture says the first shared vertex buys nothing.

**Consequences.**
- It would give Theorem 1 in the strong form H ≥ pq.
- Every gluing that keeps prime witnesses induced would then need at least two shared vertices to
  undercut the product. PRIME_GLUE.md
  ([fresh_unit_prime_gluing_20261116](../fresh_unit_prime_gluing_20261116/PRIME_GLUE.md)) proves that
  (J) needs four shared vertices for primes ≥ 11.

**A natural attack.** Splice two cycles through v. If C_A = v → α … α′ → v and C_B = v → β … β′ → v,
then one of α′ → β or β′ → α closes an odd cycle of length |C_A| + |C_B| − 1. When neither holds, the
two triangles v α β′ and v β α′ appear. So every conflicting pair (C_A, C_B) yields two mixed odd
cycles. The missing step is to make this injective; the triangle splices remember only (α, β′).

## Correction: the quantitative overlap bound of CONDITIONAL.md (MISTAKE-775)

**What CONDITIONAL.md argued.** Class 1 pairs a nonempty pack of A − W with a pack of B that uses
B − W; class 2 is the mirror. It stated that the two classes "intersect exactly in the product of
nonempty exterior packs", and printed the bound

    H ≥ p + q − r + (a−1)(q−r) + (b−1)(p−r) − (a−1)(b−1).

**What the intersection actually is.** A pack lies in both classes iff it is P ∪ P′ ∪ R, where:
- P is a nonempty pack of A − W;
- P′ is a nonempty pack of B − W;
- R is **any pack of W**.

So the intersection weighs r(a − 1)(b − 1), and the correct bound is

    H ≥ p + q − r + (a−1)(q−r) + (b−1)(p−r) − r(a−1)(b−1).

The two agree when W is transitive (r = 1).

**Counterexample to the printed bound (order 9, brute-force permutation count).**
- T = `111011111110111111111111111111111101`: `gentourng` upper triangle, row by row, 1 meaning i → j.
- A = {0..5} and B = {1, 3, 5, 6, 7, 8}. W = {1, 3, 5}, A − W and B − W are directed triangles.
- p = 33, q = 9, r = a = b = 3, and H(T) = 99.
- The printed bound gives 107; the corrected one gives 99, which is tight.

**Census.**
- The printed bound fails 135 times among the 1,787,030,880 covers at order 9. It fails nowhere at
  orders ≤ 8: the first failure needs |W|, |A − W|, |B − W| ≥ 3.
- The corrected bound fails nowhere at orders ≤ 9 (`cover_census_summary.out`).

**What stands.** CONDITIONAL.md's corollary r ≥ (2q − 1)/3 (when H ≤ p + q + 1 and a ≥ 3) uses class 1
only and is unaffected.
