# LEM-001: Reduction of c₉(T₁₁) to Sub-Tournament Hamiltonian Cycles

**Type:** Lemma
**Certainty:** 5 — PROVED (relies only on vertex-transitivity of T₁₁ and orbit counting)
**Status:** PROVED
**Last reviewed:** kind-pasteur-2026-03-05-S1
**Disputes:** none
**Tags:** #paley-tournament #c9 #hamiltonian-cycles #sub-tournament #symmetry

---

## Setup

Let T₁₁ be the Paley tournament on 11 vertices (ℤ/11ℤ, arcs i→j iff j−i is a quadratic residue mod 11; QR = {1,3,4,5,9}).

For a pair {a,b} ⊂ V with |{a,b}| = 2, let h({a,b}) denote the number of directed Hamiltonian cycles in T₁₁\{a,b} (a 9-vertex sub-tournament). Each 9-cycle of T₁₁ uses exactly 9 vertices and corresponds to a unique missing pair {a,b}.

## Statement

$$c_9(T_{11}) = \sum_{\{a,b\} \in \binom{V}{2}} h(\{a,b\}) = \frac{55}{2}\bigl(h_{\text{QR}} + h_{\text{NQR}}\bigr)$$

where:
- h_QR = h({0,1}) = number of directed Ham cycles in T₁₁\{0,1}
- h_NQR = h({0,2}) = number of directed Ham cycles in T₁₁\{0,2}

**Integrality condition:** c₉(T₁₁) is an integer, so h_QR + h_NQR must be even.

## Proof

**Step 1 (double-counting):** Every 9-cycle corresponds uniquely to a pair {a,b} of missing vertices, so
$$c_9 = \sum_{\{a,b\}} h(\{a,b\}) = \frac{1}{2}\sum_a \sum_{b \neq a} h(\{a,b\}) = \frac{11}{2}\sum_{b=1}^{10} h(\{0,b\})$$
using vertex-transitivity (translation x ↦ x−a fixes the Paley structure).

**Step 2 (multiplier orbits):** The multiplier group M = {x ↦ αx : α ∈ QR} ≅ ℤ/5ℤ acts on T₁₁ as automorphisms. It maps T₁₁\{0,b} isomorphically to T₁₁\{0,αb}. Since M has order 5 and acts freely on ℤ/11ℤ\{0}, the orbits of {1,...,10} under M are exactly:
- QR orbit: {1,3,4,5,9} (all quadratic residues)
- NQR orbit: {2,6,7,8,10} (all non-residues)

Therefore h({0,b}) takes exactly two values:
- h_QR for b ∈ {1,3,4,5,9}
- h_NQR for b ∈ {2,6,7,8,10}

**Step 3:** Substituting:
$$c_9 = \frac{11}{2}[5 h_{\text{QR}} + 5 h_{\text{NQR}}] = \frac{55}{2}(h_{\text{QR}} + h_{\text{NQR}})$$
□

## Key Consequence for the H(T₁₁) Conjecture — RESOLVED

**UPDATE (kind-pasteur-2026-03-05-S2):** CONJ-002 is REFUTED for p=11.
Direct computation gives h_QR = h_NQR = 201, so h_QR + h_NQR = 402.
Therefore c₉(T₁₁) = (55/2)(402) = 11055. And H(T₁₁) = 95095 (not 4455).

The constraint "h_QR + h_NQR ≤ 8" was derived from the false conjecture and is itself false.

## Verification Scripts

- `04-computation/paley_h_sequence.py` — computes H(P(p)) for Paley primes; QR/NQR orbit structure
- `04-computation/paley_maximizer_test.py` — Paley tournament symmetry and maximizer properties
- `04-computation/paley_deletion_test.py` — T₁₁ vertex deletion structure

## References

- Source: inbox/processed/2026-03-05/new/PALEY_T11_c9_ANALYSIS.md (Section III)
- Context: CONJ-002, LEM-002
