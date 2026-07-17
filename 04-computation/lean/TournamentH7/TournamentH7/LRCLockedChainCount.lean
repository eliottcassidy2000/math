/-
  TournamentH7.LRCLockedChainCount — THE LOCKED-CHAIN JOINT COUNT
  (death-star-2026-07-17-S47, HYP-7233; the liveCount-floor lane).

  THM-960's rung lock collapses JOINT band failures on exact-ratio pairs and
  chains to a SINGLE NARROW BAND on the bottom runner, making the joint counts
  of the deviation ledger EXACTLY computable:

  * `locked_pair_fail_iff` — for ratio `M ≤ 13` (witness form): both runners
    fail ⟺ the bottom runner fails the `14·M`-narrow band.  The (⇐) direction
    is monotonicity; the (⇒) direction is the rung lock.
  * `locked_chain_fail_iff3` — three members `u, M₂u, M₃u` (`M₂ ≤ M₃ ≤ 13`):
    all fail ⟺ the bottom fails the `14·M₃`-narrow band.  Only the bottom-top
    lock is needed; middle members ride along by monotonicity.
  * `card_mod_filter_eq` — THE MOD TRANSPORT (general): for `gcd(u, q) = 1`,
    `p ↦ (u·p) mod q` is a bijection of `(0, q)` onto itself, so filtering by
    ANY predicate of the residue has the same count as filtering the residues
    directly.  (Factored, predicate-agnostic form of the THM-942A bijection.)
  * `narrowFailN_count` — the residue-space count: for `1 ≤ M`,
    `#{r ∈ (0,q) : 14Mr < q ∨ 14M(q−r) < q} = 2·⌊(q−1)/(14M)⌋`.
  * `locked_pair_count` — the two composed: at coprime moduli the locked
    joint-failure count is EXACTLY `2·⌊(q−1)/(14M)⌋`.

  THE DEVIATION LAW this yields (recon `lockedchain_recon_deathstar_S47.out`):
  `D_pair(M) = 2⌊(q−1)/(14M)⌋ − (q−1)/49` — positive for `M < 7`, ≈0 at
  `M = 7` (the equilibrium ratio!), negative for `8 ≤ M ≤ 13`.  The pair rung
  of the B5 ledger is UNIFORM on locked strata — no per-stratum tables.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCRungLock

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **Locked pair, witness form**: for exact ratio `1 ≤ M ≤ 13`, the pair
fails jointly iff the bottom runner fails the `14M`-narrow band. -/
theorem locked_pair_fail_iff (u M : ℤ) (q p : ℕ)
    (hM1 : 1 ≤ M) (hM13 : M ≤ 13) :
    ((∃ w₁ : ℤ, 14 * |u * (p : ℤ) - w₁ * q| < q) ∧
     (∃ w₂ : ℤ, 14 * |(M * u) * (p : ℤ) - w₂ * q| < q))
      ↔ (∃ w : ℤ, 14 * M * |u * (p : ℤ) - w * q| < q) := by
  constructor
  · rintro ⟨⟨w₁, h₁⟩, ⟨w₂, h₂⟩⟩
    have hlock : w₂ = M * w₁ :=
      rung_lock u (M * u) w₁ w₂ M q p hM1 hM13 rfl h₁ h₂
    refine ⟨w₁, ?_⟩
    have hrw : (M * u) * (p : ℤ) - w₂ * q = M * (u * p - w₁ * q) := by
      rw [hlock]; ring
    rw [hrw, abs_mul, abs_of_pos (by linarith : (0:ℤ) < M)] at h₂
    calc 14 * M * |u * (p : ℤ) - w₁ * q|
        = 14 * (M * |u * (p : ℤ) - w₁ * q|) := by ring
      _ < q := h₂
  · rintro ⟨w, hw⟩
    have habs : (0 : ℤ) ≤ |u * (p : ℤ) - w * q| := abs_nonneg _
    constructor
    · refine ⟨w, ?_⟩
      nlinarith [hw, habs, hM1]
    · refine ⟨M * w, ?_⟩
      have hrw : (M * u) * (p : ℤ) - (M * w) * q = M * (u * p - w * q) := by ring
      rw [hrw, abs_mul, abs_of_pos (by linarith : (0:ℤ) < M)]
      calc 14 * (M * |u * (p : ℤ) - w * q|)
          = 14 * M * |u * (p : ℤ) - w * q| := by ring
        _ < q := hw

/-- **Locked chain of three**: only the bottom-top lock is needed; the middle
member rides along by monotonicity. -/
theorem locked_chain_fail_iff3 (u M₂ M₃ : ℤ) (q p : ℕ)
    (hM₂1 : 1 ≤ M₂) (hM₂₃ : M₂ ≤ M₃) (hM₃13 : M₃ ≤ 13) :
    ((∃ w₁ : ℤ, 14 * |u * (p : ℤ) - w₁ * q| < q) ∧
     (∃ w₂ : ℤ, 14 * |(M₂ * u) * (p : ℤ) - w₂ * q| < q) ∧
     (∃ w₃ : ℤ, 14 * |(M₃ * u) * (p : ℤ) - w₃ * q| < q))
      ↔ (∃ w : ℤ, 14 * M₃ * |u * (p : ℤ) - w * q| < q) := by
  have hM₃1 : 1 ≤ M₃ := le_trans hM₂1 hM₂₃
  constructor
  · rintro ⟨h₁, _, h₃⟩
    exact (locked_pair_fail_iff u M₃ q p hM₃1 hM₃13).mp ⟨h₁, h₃⟩
  · rintro ⟨w, hw⟩
    have habs : (0 : ℤ) ≤ |u * (p : ℤ) - w * q| := abs_nonneg _
    have hnarrow₂ : ∃ w' : ℤ, 14 * M₂ * |u * (p : ℤ) - w' * q| < q := by
      refine ⟨w, ?_⟩
      nlinarith [hw, habs, hM₂₃]
    obtain ⟨hbot, htop⟩ := (locked_pair_fail_iff u M₃ q p hM₃1 hM₃13).mpr ⟨w, hw⟩
    obtain ⟨_, hmid⟩ :=
      (locked_pair_fail_iff u M₂ q p hM₂1 (le_trans hM₂₃ hM₃13)).mpr hnarrow₂
    exact ⟨hbot, hmid, htop⟩

/-- **THE MOD TRANSPORT**: at coprime `(u, q)`, filtering `(0, q)` by any
predicate of the residue `(u·p) mod q` counts the same as filtering the
residues directly.  (Predicate-agnostic form of the THM-942A unit bijection.) -/
theorem card_mod_filter_eq (u : ℤ) (q : ℕ) (hq : 0 < q)
    (hgcd : Int.gcd u (q : ℤ) = 1) (P : ℕ → Prop) [DecidablePred P] :
    ((Finset.Ioo 0 q).filter fun p : ℕ => P ((u * (p : ℤ) % (q : ℤ)).toNat)).card
      = ((Finset.Ioo 0 q).filter fun r => P r).card := by
  set f : ℕ → ℕ := fun p => ((u * (p : ℤ)) % (q : ℤ)).toNat with hf
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hmod_lt : ∀ p : ℕ, ((u * (p : ℤ)) % (q : ℤ)) < (q : ℤ) :=
    fun p => Int.emod_lt_of_pos _ hqZ
  have hmod_nonneg : ∀ p : ℕ, (0 : ℤ) ≤ ((u * (p : ℤ)) % (q : ℤ)) :=
    fun p => Int.emod_nonneg _ (by omega)
  have hfval : ∀ p : ℕ, ((f p : ℕ) : ℤ) = (u * (p : ℤ)) % (q : ℤ) := by
    intro p
    rw [hf]
    exact Int.toNat_of_nonneg (hmod_nonneg p)
  obtain ⟨a, b, hab⟩ := Int.isCoprime_iff_gcd_eq_one.mpr hgcd
  have hva : u * a ≡ 1 [ZMOD (q : ℤ)] :=
    Int.modEq_iff_dvd.mpr ⟨b, by linarith⟩
  apply Finset.card_bij (fun p _ => f p)
  · -- maps into the residue filter
    intro p hp
    rw [Finset.mem_filter] at hp
    obtain ⟨hpIoo, hpP⟩ := hp
    rw [Finset.mem_Ioo] at hpIoo
    rw [Finset.mem_filter, Finset.mem_Ioo]
    refine ⟨⟨?_, ?_⟩, hpP⟩
    · rcases Nat.eq_zero_or_pos (f p) with h0 | hpos
      · exfalso
        have hzero : (u * (p : ℤ)) % (q : ℤ) = 0 := by
          have hv := hfval p
          rw [h0] at hv
          simpa using hv.symm
        have hdvd : (q : ℤ) ∣ u * (p : ℤ) := Int.dvd_of_emod_eq_zero hzero
        have hp_eq : (p : ℤ) = a * (u * (p : ℤ)) + (b * (p : ℤ)) * (q : ℤ) := by
          linear_combination (-(p : ℤ)) * hab
        have hqp : (q : ℤ) ∣ (p : ℤ) := by
          rw [hp_eq]
          exact dvd_add (Dvd.dvd.mul_left hdvd a) (Dvd.intro_left _ rfl)
        have hppos : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hpIoo.1
        have hple := Int.le_of_dvd hppos hqp
        have hplt : (p : ℤ) < (q : ℤ) := by exact_mod_cast hpIoo.2
        omega
      · exact hpos
    · have h1 : ((f p : ℕ) : ℤ) < (q : ℤ) := by rw [hfval]; exact hmod_lt p
      exact_mod_cast h1
  · -- injective
    intro p₁ hp₁ p₂ hp₂ hEq
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp₁ hp₂
    have hmodeq : u * (p₁ : ℤ) ≡ u * (p₂ : ℤ) [ZMOD (q : ℤ)] := by
      show (u * (p₁ : ℤ)) % (q : ℤ) = (u * (p₂ : ℤ)) % (q : ℤ)
      have hc : ((f p₁ : ℕ) : ℤ) = ((f p₂ : ℕ) : ℤ) := by exact_mod_cast hEq
      rw [hfval, hfval] at hc
      exact hc
    have hpq : (p₁ : ℤ) ≡ (p₂ : ℤ) [ZMOD (q : ℤ)] := by
      calc (p₁ : ℤ) = 1 * (p₁ : ℤ) := (one_mul _).symm
        _ ≡ (u * a) * (p₁ : ℤ) [ZMOD (q : ℤ)] := (hva.symm).mul_right _
        _ = a * (u * (p₁ : ℤ)) := by ring
        _ ≡ a * (u * (p₂ : ℤ)) [ZMOD (q : ℤ)] := hmodeq.mul_left _
        _ = (u * a) * (p₂ : ℤ) := by ring
        _ ≡ 1 * (p₂ : ℤ) [ZMOD (q : ℤ)] := hva.mul_right _
        _ = (p₂ : ℤ) := one_mul _
    have h1 : (p₁ : ℤ) % (q : ℤ) = (p₁ : ℤ) :=
      Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hp₁.1.2)
    have h2 : (p₂ : ℤ) % (q : ℤ) = (p₂ : ℤ) :=
      Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hp₂.1.2)
    have hfin : (p₁ : ℤ) = (p₂ : ℤ) := by
      have hpq' : (p₁ : ℤ) % (q : ℤ) = (p₂ : ℤ) % (q : ℤ) := hpq
      rw [h1, h2] at hpq'
      exact hpq'
    exact_mod_cast hfin
  · -- surjective
    intro r hr
    rw [Finset.mem_filter, Finset.mem_Ioo] at hr
    obtain ⟨⟨hr0, hrq⟩, hrP⟩ := hr
    set p : ℕ := ((a * (r : ℤ)) % (q : ℤ)).toNat with hpdef
    have hpmod_nonneg : (0 : ℤ) ≤ (a * (r : ℤ)) % (q : ℤ) := Int.emod_nonneg _ (by omega)
    have hpmod_lt : (a * (r : ℤ)) % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
    have hpval : ((p : ℕ) : ℤ) = (a * (r : ℤ)) % (q : ℤ) := by
      rw [hpdef]
      exact Int.toNat_of_nonneg hpmod_nonneg
    have hp_ar : (p : ℤ) ≡ a * (r : ℤ) [ZMOD (q : ℤ)] := by
      show (p : ℤ) % (q : ℤ) = (a * (r : ℤ)) % (q : ℤ)
      rw [hpval]
      exact Int.emod_emod_of_dvd _ dvd_rfl
    have hvp : u * (p : ℤ) ≡ (r : ℤ) [ZMOD (q : ℤ)] := by
      calc u * (p : ℤ)
          ≡ u * (a * (r : ℤ)) [ZMOD (q : ℤ)] := hp_ar.mul_left _
        _ = (u * a) * (r : ℤ) := by ring
        _ ≡ 1 * (r : ℤ) [ZMOD (q : ℤ)] := hva.mul_right _
        _ = (r : ℤ) := one_mul _
    have hfp : f p = r := by
      have hfpZ : ((f p : ℕ) : ℤ) = (r : ℤ) := by
        rw [hfval]
        have hrsmall : (r : ℤ) % (q : ℤ) = (r : ℤ) :=
          Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hrq)
        have hvp' : (u * (p : ℤ)) % (q : ℤ) = (r : ℤ) % (q : ℤ) := hvp
        rw [hrsmall] at hvp'
        exact hvp'
      exact_mod_cast hfpZ
    refine ⟨p, ?_, hfp⟩
    rw [Finset.mem_filter, Finset.mem_Ioo]
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · rcases Nat.eq_zero_or_pos p with h0 | hpos
      · exfalso
        have hr' : f p = r := hfp
        rw [h0] at hr'
        have hf0 : f 0 = 0 := by
          rw [hf]
          simp
        omega
      · exact hpos
    · have hlt : ((p : ℕ) : ℤ) < (q : ℤ) := by rw [hpval]; exact hpmod_lt
      exact_mod_cast hlt
    · show P (f p)
      rw [hfp]
      exact hrP

/-- **The residue-space narrow count**: `#{r ∈ (0,q) : 14Mr < q ∨ 14M(q−r) < q}
= 2·⌊(q−1)/(14M)⌋` for `1 ≤ M`. -/
theorem narrowFailN_count (M q : ℕ) (hM : 1 ≤ M) (hq : 0 < q) :
    ((Finset.Ioo 0 q).filter fun r =>
        14 * M * r < q ∨ 14 * M * (q - r) < q).card
      = 2 * ((q - 1) / (14 * M)) := by
  set b : ℕ := (q - 1) / (14 * M) with hb_def
  have hMq : 0 < 14 * M := by omega
  have hdm : 14 * M * b + (q - 1) % (14 * M) = q - 1 := by
    rw [hb_def]
    exact Nat.div_add_mod (q - 1) (14 * M)
  have hmodlt : (q - 1) % (14 * M) < 14 * M := Nat.mod_lt (q - 1) hMq
  have hb1 : 14 * M * b ≤ q - 1 := by linarith [hdm, Nat.zero_le ((q - 1) % (14 * M))]
  have hb2 : q - 1 < 14 * M * (b + 1) := by
    have hexp : 14 * M * (b + 1) = 14 * M * b + 14 * M := by ring
    linarith [hdm, hmodlt, hexp]
  have h14b : 14 * b ≤ q - 1 := by
    have h := Nat.mul_le_mul_right b (by omega : 14 ≤ 14 * M)
    omega
  -- characterize the two arcs
  have hlow : ∀ r, (r ∈ Finset.Ioo 0 q ∧ 14 * M * r < q) ↔ (1 ≤ r ∧ r ≤ b) := by
    intro r
    constructor
    · rintro ⟨hIoo, hlt⟩
      rw [Finset.mem_Ioo] at hIoo
      refine ⟨hIoo.1, ?_⟩
      by_contra hcon
      push Not at hcon
      have hmono : 14 * M * (b + 1) ≤ 14 * M * r :=
        Nat.mul_le_mul_left (14 * M) hcon
      omega
    · rintro ⟨h1, h2⟩
      have hmono : 14 * M * r ≤ 14 * M * b := Nat.mul_le_mul_left (14 * M) h2
      rw [Finset.mem_Ioo]
      refine ⟨⟨h1, ?_⟩, ?_⟩ <;> omega
  have hhigh : ∀ r, (r ∈ Finset.Ioo 0 q ∧ 14 * M * (q - r) < q)
      ↔ (q - b ≤ r ∧ r ≤ q - 1) := by
    intro r
    constructor
    · rintro ⟨hIoo, hlt⟩
      rw [Finset.mem_Ioo] at hIoo
      refine ⟨?_, by omega⟩
      by_contra hcon
      push Not at hcon
      -- r < q − b ⟹ q − r ≥ b + 1 ⟹ 14M(q−r) ≥ 14M(b+1) > q − 1 ≥ ...
      have hgeb : b + 1 ≤ q - r := by omega
      have hmono : 14 * M * (b + 1) ≤ 14 * M * (q - r) :=
        Nat.mul_le_mul_left (14 * M) hgeb
      omega
    · rintro ⟨h1, h2⟩
      have hsub : q - r ≤ b := by omega
      have hmono : 14 * M * (q - r) ≤ 14 * M * b := Nat.mul_le_mul_left (14 * M) hsub
      rw [Finset.mem_Ioo]
      refine ⟨⟨by omega, by omega⟩, by omega⟩
  -- the filter splits into the two arcs
  have hsplit : ((Finset.Ioo 0 q).filter fun r =>
      14 * M * r < q ∨ 14 * M * (q - r) < q)
      = Finset.Icc 1 b ∪ Finset.Icc (q - b) (q - 1) := by
    ext r
    simp only [Finset.mem_filter, Finset.mem_union, Finset.mem_Icc]
    constructor
    · rintro ⟨hIoo, hor⟩
      rcases hor with h | h
      · exact Or.inl ((hlow r).mp ⟨hIoo, h⟩)
      · exact Or.inr ((hhigh r).mp ⟨hIoo, h⟩)
    · rintro (h | h)
      · obtain ⟨hIoo, hlt⟩ := (hlow r).mpr h
        exact ⟨hIoo, Or.inl hlt⟩
      · obtain ⟨hIoo, hlt⟩ := (hhigh r).mpr h
        exact ⟨hIoo, Or.inr hlt⟩
  rw [hsplit]
  have hdisj : Disjoint (Finset.Icc 1 b) (Finset.Icc (q - b) (q - 1)) := by
    rw [Finset.disjoint_left]
    intro r hr1 hr2
    rw [Finset.mem_Icc] at hr1 hr2
    omega
  rw [Finset.card_union_of_disjoint hdisj, Nat.card_Icc, Nat.card_Icc]
  omega

/-- **THE LOCKED PAIR COUNT**: at coprime `(u, q)` the joint-failure count of
an exact-ratio pair with `1 ≤ M` is EXACTLY `2·⌊(q−1)/(14M)⌋`. -/
theorem locked_pair_count (u : ℤ) (M q : ℕ) (hM : 1 ≤ M) (hq : 0 < q)
    (hgcd : Int.gcd u (q : ℤ) = 1) :
    ((Finset.Ioo 0 q).filter fun p : ℕ =>
        14 * M * ((u * (p : ℤ) % (q : ℤ)).toNat) < q ∨
        14 * M * (q - ((u * (p : ℤ) % (q : ℤ)).toNat)) < q).card
      = 2 * ((q - 1) / (14 * M)) := by
  have htrans := card_mod_filter_eq u q hq hgcd
    (fun r => 14 * M * r < q ∨ 14 * M * (q - r) < q)
  rw [htrans]
  exact narrowFailN_count M q hM hq

/-! ## Axiom audit -/
#print axioms locked_pair_fail_iff
#print axioms locked_chain_fail_iff3
#print axioms card_mod_filter_eq
#print axioms narrowFailN_count
#print axioms locked_pair_count

end LRC14Concrete
end LonelyRunner
