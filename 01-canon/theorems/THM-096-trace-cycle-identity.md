# THM-096: Trace-Cycle Identity for Tournaments

**Status:** PROVED
**Proved by:** kind-pasteur-2026-03-07-S39b
**Verified computationally:** n=3-6 exhaustive (all tournaments), n=7 sampled (2000), 0 failures

---

## Statement

For any tournament T on n vertices with adjacency matrix A:

**tr(A^k) = k * c_k(T)  for k = 3, 5**

where c_k(T) is the number of directed k-cycles in T.

Equivalently:
- c_3(T) = tr(A^3) / 3
- c_5(T) = tr(A^5) / 5

This does NOT hold for k >= 7 (compound walks exist starting at n=6).

---

## Proof

In a tournament, A[i][j] * A[j][i] = 0 for all i != j (no bidirectional edges). Therefore any closed walk has length >= 3.

**Claim:** Every closed walk of length k <= 5 in a tournament is a simple directed cycle (visits k distinct vertices).

**Proof of claim:** Suppose a closed walk w_0, w_1, ..., w_k (w_k = w_0) of length k visits some vertex v at positions i and j (0 <= i < j < k). This splits the walk into two closed sub-walks:
- Walk 1: w_i, w_{i+1}, ..., w_j (length j - i)
- Walk 2: w_j, w_{j+1}, ..., w_k, w_0, ..., w_i (length k - (j - i))

Both are closed walks in the tournament, hence both have length >= 3. So k = (j-i) + (k-j+i) >= 3 + 3 = 6.

For k = 3 or 5: no vertex repetition is possible (would require k >= 6). So the walk visits exactly k distinct vertices and is a simple directed cycle.

**Counting:** Each directed k-cycle contributes exactly k to tr(A^k) (one contribution per starting vertex in the cycle). Since all contributions come from simple cycles:

tr(A^k) = k * c_k(T) for k = 3, 5. QED.

---

## Failure at k >= 7

For k = 7: a closed walk can consist of a directed 3-cycle and directed 4-cycle sharing a vertex v. The walk traverses the 3-cycle, returns to v, then traverses the 4-cycle. This visits 6 distinct vertices but has length 7.

**Example (n=7, bits=4):** c_3=2, c_4=1, c_7=0 but tr(A^7)=14. The excess 14 = 2 compound walk arrangements * 7 starting positions.

For k = 6: a walk can consist of two directed 3-cycles sharing a vertex (3+3=6 edges, 5 vertices), confirmed by tr(A^6) != 6*c_6 at n=5.

---

## Computational Consequence

**c_5(T) = tr(A^5) / 5** gives an O(n^3) algorithm for counting directed 5-cycles (via matrix multiplication), replacing the O(n^5) subset-DP enumeration.

This provides a 2.8x speedup at n=7 and grows to O(n^2) speedup at larger n.

Combined with Moon's formula c_3 = C(n,3) - sum_v C(s_v, 2), both c_3 and c_5 are now computable in polynomial time without explicit cycle enumeration.

---

## Corollary: Moon's Formula via Trace

Moon's classical formula c_3(T) = C(n,3) - sum_v C(s_v, 2) is equivalent to tr(A^3)/3 = C(n,3) - sum_v C(s_v, 2). This connects the combinatorial (score sequence) and algebraic (adjacency matrix) descriptions of 3-cycle counts.

**Verified:** Exhaustive at n=4,5,6 — all tournaments satisfy both formulas simultaneously.

---

## Scripts

- `04-computation/trace_cycle_theorem.py`: Exhaustive verification and compound walk analysis
- `04-computation/savchenko_verify.py`: Related Savchenko formula verification
- `04-computation/c5_formula_regular.py`: c_5 analysis for regular tournaments
