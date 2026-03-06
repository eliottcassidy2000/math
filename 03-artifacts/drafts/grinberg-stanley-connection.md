# Connection to Grinberg-Stanley (2023/2024)

**Discovered by:** kind-pasteur-2026-03-05-S12 (comprehensive audit)
**Status:** CONFIRMED — OCF is equivalent to their Corollary 20

---

## The Discovery

During a comprehensive audit of the project's mathematical results, a web search
revealed that the OCF formula H(T) = I(Omega(T), 2) is equivalent to a result
already proved in the literature:

**Corollary 20 of arXiv:2412.10572** (Grinberg & Stanley, December 2024):

> For a tournament D on [n]:
> ham(D̄) = Σ_{σ ∈ S(D), all cycles of σ have odd length} 2^{ψ(σ)}

where:
- S(D) = set of permutations whose nontrivial cycles are directed cycles of D
- ψ(σ) = number of nontrivial cycles of σ
- D̄ = complement digraph (for tournaments, D̄ = D^op, the converse)

---

## Why This Is Exactly OCF

### Step 1: D̄ = D^op for tournaments

For a tournament D, each pair {u,v} has exactly one arc. The complement D̄ has
(u,v) iff (u,v) ∉ D, which (since exactly one of (u,v),(v,u) is in D) means
(u,v) ∈ D̄ iff (v,u) ∈ D. So D̄ = D^op.

### Step 2: ham(D^op) = ham(D)

Path reversal: (v₁, v₂, ..., vₙ) is a Hamiltonian path in D iff
(vₙ, ..., v₂, v₁) is a Hamiltonian path in D^op. This is a bijection.

### Step 3: The RHS equals I(Omega(D), 2)

Each permutation σ ∈ S(D) with all nontrivial cycles odd corresponds to a
collection of vertex-disjoint odd directed cycles of D (the nontrivial cycles),
plus fixed points. The map σ → {nontrivial cycles of σ} is a bijection between:
- {σ ∈ S(D) : all cycles odd} ↔ {collections of pairwise v.d. odd directed cycles of D}

The weight 2^{ψ(σ)} = 2^{|collection|}.

Independent sets in Omega(D) are exactly collections of pairwise vertex-disjoint
odd directed cycles (since Omega(D) has vertices = odd directed cycles, edges when
two cycles share a vertex). So:

RHS = Σ_k α_k(Omega(D)) · 2^k = I(Omega(D), 2)

### Step 4: Conclusion

ham(D) = ham(D̄) = I(Omega(D), 2) = H(D).

This is OCF.

---

## The Original Paper

The result originates in:

1. **arXiv:2307.05569** (Grinberg & Stanley, "The Rédei-Berge symmetric function
   of a directed graph", July 2023). This paper introduces the Rédei-Berge
   symmetric function RB(D) and derives tournament identities from it.

2. **arXiv:2412.10572** (Grinberg & Stanley, "Revisiting The Rédei-Berge
   Symmetric Functions via Matrix Algebra", December 2024). This follow-up
   re-derives the results using matrix algebra and states Corollary 20 explicitly.

The proof technique is entirely different from any approach attempted in this
project: it uses symmetric functions, matrix algebra, and the Rédei-Berge
framework, rather than combinatorial bijections, arc-flip invariance, or
per-path analysis.

---

## What This Means for the Project

### Results that are now PROVED (were previously open or conditional)

| Result | Previous Status | New Status |
|--------|----------------|------------|
| OCF: H(T) = I(Omega(T), 2) | Conditional on Claim A | PROVED (Grinberg-Stanley) |
| Claim A: H(T)-H(T-v) = 2Σμ(C) | OPEN (verified n≤8) | PROVED (OCF + Claim B) |
| PROP-001: E(T) invariant | OPEN | PROVED (E(T) = 0 identically) |
| OPEN-Q-002 | OPEN | RESOLVED |
| OPEN-Q-009 | Partial progress | RESOLVED |

### Results that remain independently valuable

| Result | Why it's still valuable |
|--------|----------------------|
| THM-016/017 (even-odd split) | Proved for all n by our own induction; independent of G-S |
| THM-018 (coefficient identity) | Provides structural insight into WHY OCF holds per-vertex |
| THM-015 (polynomial identity) | Direct computational verification, independent check |
| Tribonacci theorem (tiling) | New result not in G-S; characterizes specific tournament family |
| Paley computations | H(T_11)=95095, H(T_19)=1,172,695,746,915 — new data |
| Omega(T) perfectness | Not in G-S; potential new theorem |
| Tiling class structure | Entirely new framework not related to G-S |

### What the project adds beyond Grinberg-Stanley

1. **The independence polynomial / conflict graph formulation.** G-S use permutation
   cycle language; our Omega(T) conflict graph formulation is more graph-theoretic
   and connects to the hard-core lattice gas model at fugacity 2.

2. **The mu-weighted vertex deletion recurrence.** Claim A gives a recursive
   structure that G-S don't explore. The mu weights mu(C) = H(T[V\V(C)]) provide
   a way to compute H(T) recursively by vertex deletion.

3. **The even-odd split and signed adjacency identity.** THM-016/017 are new
   identities about tournament Hamiltonian paths (proved for all n) that don't
   appear in G-S.

4. **Computational data.** Exhaustive verification tables, Paley tournament
   values, tiling class structure, Tribonacci connection — all new.

---

## Lesson for Future Work

The project independently discovered and partially proved a result that was
established in the literature by Grinberg and Stanley (2023-2024). The symmetric
function approach was entirely different from our combinatorial/computational
methods. This underscores the importance of literature searches early in a
research project.

However, the project's formulation (conflict graph + independence polynomial)
and the structural theorems (THM-016, THM-018, even-odd split) provide genuine
new insight that is not present in the Grinberg-Stanley papers.
