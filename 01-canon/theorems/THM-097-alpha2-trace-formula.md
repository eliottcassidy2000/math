# THM-097: Alpha_2 Trace Formula

**Status:** PROVED
**Proved by:** kind-pasteur-2026-03-07-S39b
**Verified computationally:** n=5,6 exhaustive (all tournaments), n=7 sampled (500), 0 failures

---

## Statement

For any tournament T on n vertices with adjacency matrix A:

**alpha_2(Omega_3(T)) = C(c_3, 2) - sum_v C(t_3(v), 2) + s_2**

where:
- c_3 = number of directed 3-cycles (via Moon's formula: c_3 = C(n,3) - sum_v C(s_v, 2))
- t_3(v) = (A^3)[v][v] = number of directed 3-cycles through vertex v
- s_2 = sum over directed edges (a->b) of C((A^2)[b][a], 2) = number of 3-cycle pairs sharing exactly 2 vertices
- alpha_2 = number of vertex-disjoint pairs of directed 3-cycles

All quantities are computable in O(n^3) via matrix multiplication.

For n <= 7: alpha_2(Omega_3) = alpha_2(Omega) (full conflict graph), since no 5-cycle can be vertex-disjoint from any 3-cycle (would need 5+3=8 > 7 vertices).

---

## Proof

Partition pairs of 3-cycles by overlap count:
- s_0 = pairs sharing 0 vertices (vertex-disjoint) = alpha_2
- s_1 = pairs sharing exactly 1 vertex
- s_2 = pairs sharing exactly 2 vertices

Total pairs: C(c_3, 2) = s_0 + s_1 + s_2.

**Claim 1:** sum_v C(t_3(v), 2) = s_1 + 2*s_2.

*Proof:* Each pair sharing exactly k vertices (k=1 or 2) contributes C(k,1) = k to the sum (counted once per shared vertex). So sum = 1*s_1 + 2*s_2.

Note: sharing 3 vertices is impossible in tournaments (unique directed cycle per vertex triple).

**Claim 2:** s_2 = sum_{directed edges a->b} C((A^2)[b][a], 2).

*Proof:* Two 3-cycles sharing exactly 2 vertices {a,b} share a directed edge. If the shared edge is a->b (T[a][b]=1), both cycles contain the path a->b, and the third vertex c satisfies T[b][c]=T[c][a]=1. The number of such c is (A^2)[b][a]. So the number of 2-overlap pairs through edge a->b is C((A^2)[b][a], 2).

Note: T[b][a]=1 would mean the shared edge is b->a, giving the same formula with arguments swapped. Since T[a][b] + T[b][a] = 1, each undirected edge contributes once.

**Conclusion:** alpha_2 = C(c_3, 2) - (s_1 + s_2) = C(c_3, 2) - (sum_v C(t_3(v), 2) - s_2) = C(c_3, 2) - sum_v C(t_3(v), 2) + s_2. QED.

---

## Computational Consequence

Combined with THM-096 (c_3, c_5 via trace), this gives an O(n^3) formula for H(T) at n <= 7:

H(T) = 1 + 2*(c_3 + c_5 + c_7) + 4*alpha_2

where c_3 = Moon's formula (O(n^2)), c_5 = tr(A^5)/5 (O(n^3)), c_7 is computable at n=7 in O(2^7) = O(128), and alpha_2 via the formula above (O(n^3)).

At n >= 8: need additional terms for alpha_3 and for 5-cycle contributions to alpha_2.

---

## Extension to n >= 8

At n = 8: vertex-disjoint (5-cycle, 3-cycle) pairs can exist (5+3=8 vertices). Computing these disjoint pairs requires knowing both the 5-cycles and the 3-cycles and checking disjointness — currently O(n^5) for 5-cycle enumeration. With THM-096 giving c_5 = tr(A^5)/5, we know the COUNT but not the vertex sets. Extending this to alpha_2 at n >= 8 requires either:
1. A formula for the number of disjoint (3-cycle, 5-cycle) pairs in terms of matrix entries, or
2. Fall back to cycle enumeration for 5-cycles.

---

## Scripts

- `04-computation/alpha2_formula.py`: Implementation and verification
- `04-computation/trace_ocf_bridge.py`: Derivation and testing
