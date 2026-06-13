# THM-482 — One doubling step thermalizes to d₂ₙ⁺: the universal forgetting theorem

**Status:** PROVED (claudebox-2026-06-11-S3) + VERIFIED exact (8→16, 16→32 code equality;
dims and structure at 64 implied by the proof). Resolves HYP-2409(1), strengthening it.
**Provenance:** dispatch "leverage the RM/self-dual understanding". **Companions:** THM-480
(the ladder data), THM-481 (the contrasting Paley/QR regime), THM-447 (the doubling).

## Statement

Let H be ANY skew-type Hadamard matrix of order n ≡ 0 (mod 4), n ≥ 8, whose tournament-gauge
binary rows all have even weight (equivalently, every row sum ≡ n mod 4; holds e.g. whenever
C(H) is doubly even, in particular for every level of the THM-447 tower and every Paley
border with q ≡ 7 mod 8). Let H' = [[H, H], [−Hᵀ, Hᵀ]] be its skew-Sylvester double. Then

   **C(H') ≅ d₂ₙ⁺,  regardless of H.**

In particular: (a) the THM-447/448 tower satisfies C(H_{2^k}) = d_{2^k}⁺ for ALL k ≥ 4
(HYP-2409(1), now proved — base C(H₈) = ê₈ has even rows, and d⁺ codes have even rows, so
the hypothesis propagates); (b) the doubling FORGETS the input code entirely — e.g.
double(border(Paley₂₃)) has code d₄₈⁺, not the extremal eQR(48) of THM-481.

## Proof

Write b_i = bin(row_i(H)), γ_j = bin(col_j(H)) ∈ 𝔽₂ⁿ. Rows of H': the top block rows are
(r_i, r_i) ↦ (b_i, b_i); the bottom rows are (−c_j, c_j) ↦ (γ̄_j, γ_j) = (γ_j, γ_j) ⊕ (𝟙, 0).

(1) **Gauge identity.** Skew-type means Hᵀ = 2I − H, so c_j = 2e_j − r_j as ±1 vectors:
binarizing, **γ_j = 𝟙 ⊕ b_j ⊕ e_j** (complement off the diagonal; the diagonal +1 survives).

(2) **Span shape.** (γ̄_j, γ_j) ⊕ (γ̄_k, γ_k) = PD(γ_j ⊕ γ_k), where PD(x) = (x, x) is the
pair-doubling embedding (pairs = positions {i, i+n}). Hence
C(H') = PD(W) + ⟨(γ̄_1, γ_1)⟩, with W = span({b_i} ∪ {γ_j ⊕ γ_k}).

(3) **W is the full even-weight code.** By (1), γ_j ⊕ γ_k = b_j ⊕ b_k ⊕ e_j ⊕ e_k, so W =
C(H) + ⟨e_j ⊕ e_k : j < k⟩ = C(H) + E_n = E_n (dim n−1), using the even-row hypothesis
C(H) ⊆ E_n.

(4) **PD(E_n) = d₂ₙ.** PD(x) is constant on pairs with 11-pairs exactly at supp(x); even
weight ⟺ evenly many 11-pairs — the defining property of d₂ₙ.

(5) **The glue is pairwise-odd.** On each pair {i, i+n} the vector (γ̄_1, γ_1) reads
(1 ⊕ γ_1[i], γ_1[i]) — odd on every pair. The self-dual extensions of d₂ₙ are exactly
i₂ⁿ (pair-constant glue) and d₂ₙ⁺ (pairwise-odd glue), and both pairwise-odd cosets give
codes equivalent to d₂ₙ⁺. Hence C(H') = d₂ₙ ∪ (d₂ₙ + glue) ≅ d₂ₙ⁺. Weight check: the glue
has weight n ≡ 0 (mod 4), so C(H') is doubly even and the even-row hypothesis propagates. ∎

## Verification (script gleason_qr_dplus_cbx3.py)

Exact code equality C(H₁₆) = PD(E₈) + glue and C(H₃₂) = PD(E₁₆) + glue as SETS of words;
rows+cols spans confirm step (3)'s ingredients; d⁺ identification at 16/32 from THM-480.
Negative control: the even-row hypothesis is necessary — dropping it would give
W = 𝔽₂ⁿ and a non-self-orthogonal candidate, which cannot occur for skew-Hadamard
(self-orthogonality of C(H') is automatic from MMᵀ = 2nI with 8 | 2n).

## Remarks

1. **Thermalization.** One doubling step is a complete memory wipe at the code level: the
   output is the crystalline d⁺ regardless of input. THM-451's "tower splits from Sylvester
   at 16" is the FIRST application; this theorem says the split lands at d⁺ and stays. In
   the S720 temperature language: skew-Sylvester doubling is the cooling map's fixed point
   arrival in ONE step; the Paley borders (THM-481) are the arithmetic/hot family the
   doubling instantly destroys.
2. **Contrast with the ±1 level.** At the matrix level the double REMEMBERS H (different H
   give inequivalent H' in general — e.g. THM-451's tower vs Sylvester); the binary code is
   a coarser shadow that forgets in one step. The code functor collapses the doubling orbit;
   the Pfaffian/H-invariants do not. Where exactly the memory survives (SNF? 4-rank? the
   ℤ₄-code?) is the natural next question — the ℤ₄-linear lift (Nordstrom–Robinson
   territory at length 16!) is flagged in t-0119.
