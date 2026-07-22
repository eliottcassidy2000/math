import Mathlib

/-!
# The additive orbit-sum contradiction (THM-2101 core, bypassing THM-1550)

The additive Duistermaat–van der Kallen proof (THM-2101) closes with a purely algebraic identity that
avoids the small-root **product**, Hensel factorization, logarithms and Wiener–Hopf (THM-1550).  In
the paper proof, a transitive monodromy group `G` acts on the full root set `Ω`; a residue weight `w`
has full-root sum zero (the Lagrange identity, since the residue at infinity vanishes when the
positive endpoint is present), while the vanishing-constant-term assumption forces the translated
sums over a small-root subset `S` to all equal `1`.  Summing over `G` and using transitivity gives
`|G| = |Stab|·(∑_Ω w) = 0`, impossible in characteristic zero.

This module proves that final contradiction abstractly.  It is the additive analogue of the
multiplicative `thm2067_contradiction` (THM-2067), and — unlike that route — it needs **no** product,
Hensel lift, or Wiener–Hopf input.  What remains external to it is the analytic/`t`-adic construction
of the two hypotheses `hA` (translated small-root sums are `1`) and `hB` (full-root sum is `0`).
-/

open scoped BigOperators

namespace GMC2AdditiveOrbitSum

variable {K : Type*} [Field K] [CharZero K]
variable {G : Type*} [Group G] [Fintype G]
variable {Ω : Type*} [Fintype Ω] [MulAction G Ω]

/-- Summing a weight over `Ω` is invariant under the permutation `σ • ·`. -/
theorem sum_smul_eq (σ : G) (w : Ω → K) : ∑ a : Ω, w (σ • a) = ∑ b : Ω, w b :=
  Equiv.sum_comp (Equiv.ofBijective (fun a => σ • a) (MulAction.bijective σ)) w

/-- The full orbit-sum `∑_σ w (σ • α)` is independent of the base point `α` (transitivity). -/
theorem orbitSum_const (htrans : MulAction.IsPretransitive G Ω) (w : Ω → K) (α α' : Ω) :
    (∑ σ : G, w (σ • α)) = ∑ σ : G, w (σ • α') := by
  obtain ⟨σ₀, hσ₀⟩ := htrans.exists_smul_eq α α'
  subst hσ₀
  simp_rw [← mul_smul]
  exact (Fintype.sum_bijective (· * σ₀) (Group.mulRight_bijective σ₀)
    (fun σ => w ((σ * σ₀) • α)) (fun σ => w (σ • α)) (fun σ => rfl)).symm

/-- With the full-root sum zero, every orbit-sum vanishes. -/
theorem orbitSum_zero [Nonempty Ω] (htrans : MulAction.IsPretransitive G Ω) (w : Ω → K)
    (hB : ∑ β : Ω, w β = 0) (α : Ω) : (∑ σ : G, w (σ • α)) = 0 := by
  have key : (Fintype.card Ω) • (∑ σ : G, w (σ • α)) = 0 := by
    calc (Fintype.card Ω) • (∑ σ : G, w (σ • α))
        = ∑ _α' : Ω, (∑ σ : G, w (σ • α)) := by rw [Finset.sum_const, Finset.card_univ]
      _ = ∑ α' : Ω, (∑ σ : G, w (σ • α')) :=
            Finset.sum_congr rfl (fun α' _ => orbitSum_const htrans w α α')
      _ = ∑ σ : G, (∑ α' : Ω, w (σ • α')) := Finset.sum_comm
      _ = 0 := Finset.sum_eq_zero (fun σ _ => by rw [sum_smul_eq σ w]; exact hB)
  rw [nsmul_eq_mul] at key
  exact (mul_eq_zero.mp key).resolve_left (Nat.cast_ne_zero.mpr Fintype.card_ne_zero)

/-- **The additive orbit-sum contradiction** — the algebraic core of the additive DvdK proof
(THM-2101), bypassing THM-1550 / Wiener–Hopf.  If a transitive finite group action carries a weight
`w` with full sum zero (`hB`), then no subset `S` can have all its `G`-translated `w`-sums equal to
`1` (`hA`). -/
theorem additive_orbit_contradiction [Nonempty Ω] (htrans : MulAction.IsPretransitive G Ω)
    (w : Ω → K) (S : Finset Ω)
    (hA : ∀ σ : G, ∑ α ∈ S, w (σ • α) = 1)
    (hB : ∑ β : Ω, w β = 0) : False := by
  have hcard : (Fintype.card G : K) = 0 := by
    calc (Fintype.card G : K)
        = ∑ _σ : G, (1 : K) := by rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one]
      _ = ∑ σ : G, ∑ α ∈ S, w (σ • α) := Finset.sum_congr rfl (fun σ _ => (hA σ).symm)
      _ = ∑ α ∈ S, ∑ σ : G, w (σ • α) := Finset.sum_comm
      _ = 0 := Finset.sum_eq_zero (fun α _ => orbitSum_zero htrans w hB α)
  exact (Nat.cast_ne_zero.mpr Fintype.card_ne_zero) hcard

end GMC2AdditiveOrbitSum

#print axioms GMC2AdditiveOrbitSum.additive_orbit_contradiction
