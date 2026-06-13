# THM-082: Deletion-Contraction for Hamiltonian Path Counts

**Status:** PROVED
**Proved by:** kind-pasteur-2026-03-07-S35
**Verified computationally:** n=4 exhaustive (384 edge tests), n=5 exhaustive (10240 edge tests), 0 failures
**Scope:** All digraphs (not just tournaments)

---

## Statement

Let D be a digraph on n vertices and e = (u → v) a directed edge of D. Define:

- **D \ e** (deletion): the digraph obtained by removing edge e (same vertex set, one fewer edge).
- **D / e** (contraction): the digraph on n−1 vertices obtained by merging u and v into a single vertex w, where:
  - For x ≠ u,v: edge (x, w) exists iff (x, u) exists in D
  - For x ≠ u,v: edge (w, x) exists iff (v, x) exists in D
  - Edges among other vertices: unchanged
  - (Convention: w inherits **IN-edges from the tail** u and **OUT-edges from the head** v)

Then:

**H(D) = H(D \ e) + H(D / e)**

where H(·) denotes the number of Hamiltonian paths.

---

## Proof

Partition Ham(D) = {Hamiltonian paths of D} into two disjoint sets:

1. **Paths not using e:** These are exactly Ham(D \ e), since removing e from D doesn't affect any path that avoids e. Contributes H(D \ e).

2. **Paths using e:** We construct a bijection φ : {P ∈ Ham(D) : e ∈ P} → Ham(D / e).

### Bijection φ

Let P = (p₀, p₁, ..., pₙ₋₁) be a Hamiltonian path in D using edge e = (u → v). Then u = pₖ and v = pₖ₊₁ for some unique position k (0 ≤ k ≤ n−2).

Define φ(P) = P' = (p'₀, ..., p'ₙ₋₂) where:
- p'ᵢ = pᵢ for i < k  (vertices before u)
- p'ₖ = w            (merged vertex replaces u→v)
- p'ᵢ = pᵢ₊₁ for i > k  (vertices after v, shifted)

**φ(P) is a valid Hamiltonian path in D/e:**
- **Visits all n−1 vertices:** P visits all n vertices of D. In P', u and v are replaced by w, and all other vertices appear exactly once. ✓
- **Edge (p'ₖ₋₁, w) exists in D/e** (when k ≥ 1): p'ₖ₋₁ = pₖ₋₁, and (pₖ₋₁, u) is an edge of D (from path P). By contraction definition, (pₖ₋₁, w) exists in D/e. ✓
- **Edge (w, p'ₖ₊₁) exists in D/e** (when k ≤ n−3): p'ₖ₊₁ = pₖ₊₂, and (v, pₖ₊₂) is an edge of D (from path P). By contraction definition, (w, pₖ₊₂) exists in D/e. ✓
- **All other edges** (pᵢ, pᵢ₊₁) for i ≠ k−1, k, k+1 involve vertices other than u,v, so they exist unchanged in D/e. ✓
- **Boundary cases:** w at start (k=0) or end (k=n−2) omit the corresponding in/out check. ✓

### Inverse bijection φ⁻¹

Given P' = (p'₀, ..., p'ₙ₋₂) ∈ Ham(D/e), let k be the position where p'ₖ = w. Define:

φ⁻¹(P') = P = (p'₀, ..., p'ₖ₋₁, u, v, p'ₖ₊₁, ..., p'ₙ₋₂)

**This is valid in D:**
- Edge (p'ₖ₋₁, u): exists because (p'ₖ₋₁, w) ∈ D/e implies (p'ₖ₋₁, u) ∈ D. ✓
- Edge (u, v) = e: exists in D by hypothesis. ✓
- Edge (v, p'ₖ₊₁): exists because (w, p'ₖ₊₁) ∈ D/e implies (v, p'ₖ₊₁) ∈ D. ✓

### Conclusion

φ is a bijection, so |{P ∈ Ham(D) : e ∈ P}| = H(D/e). Together with H(D\e) = |{P ∈ Ham(D) : e ∉ P}|:

**H(D) = H(D \ e) + H(D / e)** □

---

## Remarks

1. **Convention matters:** The reverse contraction (w inherits IN from v, OUT from u) does NOT satisfy the identity. This is because the reverse convention would require v→u to be the contracted edge, not u→v.

2. **Generality:** The proof works for ANY digraph, not just tournaments. No completeness or tournament constraint is needed.

3. **Connection to Mitrovic (arXiv:2504.20968):** This is the commutative specialization of Mitrovic's noncommutative deletion-contraction W_X = W_{X\e} − W_{X/e}↑. The sign difference (+ vs −) comes from the "↑" operation (promotion/superization) in the noncommutative setting. When specialized to counting (all variables = 1), promotion becomes identity and the sign becomes +.

4. **Iteration:** Repeated application decomposes H(D) as a sum over all 2^m subsets of a chosen edge ordering, each evaluated on the corresponding minor. This is analogous to the Tutte polynomial expansion but for Hamiltonian paths.

5. **Not a matroid invariant:** Unlike the Tutte polynomial, H(D) depends on the ORDER of deletion-contraction. Different edge orderings give different intermediate digraphs (the contraction can create non-tournament digraphs even from tournaments). The final answer H(D) is of course order-independent.

6. **Connection to opus flip formula (S45):** For tournaments T with e = (u→v), flipping e gives T' with e' = (v→u). Then:
   - H(T) = H(T\e) + H(T/e)
   - H(T') = H(T'\e') + H(T'/e')
   - T\e and T'\e' have the same edge set (both remove the u-v edge)
   - But T\e ≠ T'\e' as digraphs if there are other edges between u's and v's neighborhoods... Actually T\e = T'\e' = T with edge u↔v removed (identical digraphs)
   - So H(T) - H(T') = H(T/e) - H(T'/e')
   - This reduces the arc-flip Hamiltonian path difference to a contraction-level identity on n−1 vertices.
