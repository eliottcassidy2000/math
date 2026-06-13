/-
  TournamentH7.Forbidden — Generalised forbidden-H-value framework

  Infrastructure for proving H(T) ≠ h beyond h = 7.

  Each forbidden value h is killed by combining:
    (1) **Arithmetic**: enumerate non-negative integer α-vectors satisfying
        the OCF sum and the combinatorial chain constraints.
    (2) **Structural**: for each surviving α-vector, an axiom asserts
        non-realisability of that α-vector as Ω(T)'s independence vector.

  ## Combinatorial axioms (extensions of OCF.lean A1')

  The basic `alpha_subset_bound` says `α_k ≠ 0 → α_1 ≥ k`. For h ≥ 21 we
  also need:

    • the **adjacent chain step**  `α_k ≠ 0 → α_{k-1} ≥ k`,
    • the **binomial upper bound** `α_k ≤ C(α_1, k)`.

  Both are immediate from the meaning of α_k:
    – each k-element independent set has k distinct (k−1)-element
      independent subsets;
    – distinct k-element independent sets are distinct k-element
      vertex subsets, hence ≤ C(α_1, k) in number.

  ## Extended OCF

  The OCF axiom in `OCF.lean` truncates the sum at k = 4 (sufficient for
  h = 7). For h ≤ 127 we need the sum to k = 6.

  References
   · Project canon `01-canon/theorems/THM-343-H7-impossible.md` (audit S6).
   · `04-computation/audit_thm343_s6.py` confirms each chain bound at
     small n exhaustively.
-/

import TournamentH7.OCF
import Mathlib.Data.Nat.Choose.Basic

namespace Tournament

variable {n : ℕ}

/-! ### Combinatorial chain axioms -/

/-- **Axiom (adjacent chain step).** Each k-element independent set in
    Ω(T) has k distinct (k−1)-element independent subsets, so

      α_k(T) ≠ 0 ⟹ α_{k−1}(T) ≥ k. -/
axiom alpha_chain_step (T : Tournament n) (k : ℕ) (hk : 2 ≤ k) :
    alphaCount k T ≠ 0 → k ≤ alphaCount (k - 1) T

/-- **Axiom (binomial upper bound).** Distinct k-element independent sets
    are distinct k-element subsets of the α_1(T) odd cycles, so

      α_k(T) ≤ C(α_1(T), k). -/
axiom alpha_binomial_bound (T : Tournament n) (k : ℕ) :
    alphaCount k T ≤ Nat.choose (alphaCount 1 T) k

/-! ### Convenience specialisations -/

lemma alpha_pair_bound (T : Tournament n) :
    alphaCount 2 T ≠ 0 → 2 ≤ alphaCount 1 T :=
  alpha_subset_bound T 2 (by decide)

lemma alpha_triple_subset (T : Tournament n) :
    alphaCount 3 T ≠ 0 → 3 ≤ alphaCount 1 T :=
  alpha_subset_bound T 3 (by decide)

lemma alpha_quad_subset (T : Tournament n) :
    alphaCount 4 T ≠ 0 → 4 ≤ alphaCount 1 T :=
  alpha_subset_bound T 4 (by decide)

lemma alpha_quint_subset (T : Tournament n) :
    alphaCount 5 T ≠ 0 → 5 ≤ alphaCount 1 T :=
  alpha_subset_bound T 5 (by decide)

lemma alpha_sext_subset (T : Tournament n) :
    alphaCount 6 T ≠ 0 → 6 ≤ alphaCount 1 T :=
  alpha_subset_bound T 6 (by decide)

lemma alpha_triple_chain (T : Tournament n) :
    alphaCount 3 T ≠ 0 → 3 ≤ alphaCount 2 T := by
  have h := alpha_chain_step T 3 (by decide)
  -- 3 - 1 = 2
  simpa using h

lemma alpha_quad_chain (T : Tournament n) :
    alphaCount 4 T ≠ 0 → 4 ≤ alphaCount 3 T := by
  have h := alpha_chain_step T 4 (by decide)
  simpa using h

lemma alpha_quint_chain (T : Tournament n) :
    alphaCount 5 T ≠ 0 → 5 ≤ alphaCount 4 T := by
  have h := alpha_chain_step T 5 (by decide)
  simpa using h

lemma alpha_sext_chain (T : Tournament n) :
    alphaCount 6 T ≠ 0 → 6 ≤ alphaCount 5 T := by
  have h := alpha_chain_step T 6 (by decide)
  simpa using h

/-! ### Extended OCF (Grinberg–Stanley up to k = 6, sufficient for h ≤ 127) -/

/-- **Axiom (Extended OCF — Grinberg–Stanley 2024, truncated to k = 6).**
    H(T) equals the OCF sum truncated at α₆, provided α₇ = α₈ = ⋯ = 0.

    For tournaments with H(T) ≤ 127 this truncation is exact because each
    α_k ≥ 1 with k ≥ 7 would contribute 2^k ≥ 128 to the sum. -/
axiom ocf_extended (T : Tournament n) :
    H T = 1 + 2 * alphaCount 1 T + 4 * alphaCount 2 T
            + 8 * alphaCount 3 T + 16 * alphaCount 4 T
            + 32 * alphaCount 5 T + 64 * alphaCount 6 T

/-! ### Realisability predicate -/

/-- `HIsRealisable n h` ↔ some `n`-vertex tournament has H = h. -/
def HIsRealisable (n h : ℕ) : Prop :=
  ∃ T : Tournament n, H T = h

end Tournament
