/-
# LRCSchurPeel — the peeling recursion for the Schur count (E3-axis Freiman-stability ladder)

Peeling the **maximum** element `m` of `S` from the Schur count `E₃`:

  `E₃ S = E₃ (S.erase m) + repCount S m`,   `repCount S m = #{(a,b) ∈ S² : a+b = m}`,

because a Schur pair `(a,b)` with `a+b ∈ S` can never involve `m` (that would force `a+b > max S`), so
the *only* Schur pairs lost when `m` leaves are those **summing to** `m`.  With `repCount S m ≤ |S|-1`
(the map `(a,b) ↦ a` injects the representations into `S.erase m`) this yields the **deficit recursion**

  `deficit S + repCount S m = deficit (S.erase m) + (|S|-1)`,   `deficit S := C(|S|,2) − E₃ S`,

so `deficit` is monotone under peeling (`deficit_erase_le`), the **peel cost** `(|S|-1) − repCount S m`
is `≥ 0` and bounded by the total deficit (`peelCost_le_deficit`), and `deficit S = 0 ⟺ dilated`
(`deficit_eq_zero_iff_dilated`).  These are the rungs of the E3-side Freiman-stability ladder whose top
(`deficit 0`) is the dilated-interval rigidity `schurCount_eq_choose_iff_dilated`, and whose empirical
capstone `dist_to_dilated S ≤ deficit S` (verified exhaustively, see the S126 note) sits on this skeleton:
the total deficit is the sum of the per-peel costs, so a small deficit means few deficient peels.

The E3 axis (`repCount`/`E₃`, anchored at the origin: Schur incidences `a+b=c`) is the companion of the
burden axis (opus's translation-invariant `restrictedSum`); both feed THM-681's `W₀ > 0.08` branch.

kind-pasteur-2026-07-09-S126.
-/
import Mathlib
import TournamentH7.LRCSchurRigidity

namespace LRCSchurRigidity

open Finset

/-- **Ordered additive representations of `m` inside `S`:** `#{(a,b) ∈ S² : a+b = m}`. -/
def repCount (S : Finset ℕ) (m : ℕ) : ℕ := ((S ×ˢ S).filter (fun p => p.1 + p.2 = m)).card

/-- **The Schur deficit** `C(k,2) − E₃` — the number of unrealised 2-subsets (missing differences),
by `schurCount_add_sdiff_eq_choose`. -/
def deficit (S : Finset ℕ) : ℕ := S.card.choose 2 - E3 S

/-- A Schur pair never involves the maximum: if `a,b ∈ S`, `a+b ∈ S` and `m` bounds `S`, then
`a ≠ m` and `b ≠ m` (else `a+b > m ≥` everything in `S`, so `a+b ∉ S`). -/
lemma ne_max_of_sum_mem {S : Finset ℕ} (h0 : 0 ∉ S) {m : ℕ} (hmax : ∀ x ∈ S, x ≤ m)
    {a b : ℕ} (ha : a ∈ S) (hb : b ∈ S) (hab : a + b ∈ S) : a ≠ m ∧ b ≠ m := by
  have ha0 : 0 < a := Nat.pos_of_ne_zero (fun h => h0 (h ▸ ha))
  have hb0 : 0 < b := Nat.pos_of_ne_zero (fun h => h0 (h ▸ hb))
  have hle : a + b ≤ m := hmax _ hab
  exact ⟨by omega, by omega⟩

/-- `E₃ ≤ C(k,2)` — the Schur count never exceeds the number of 2-subsets. -/
theorem E3_le_choose (S : Finset ℕ) (h0 : 0 ∉ S) : E3 S ≤ S.card.choose 2 := by
  have := schurCount_add_sdiff_eq_choose S h0
  omega

/-- **Rung A — the E₃ peel recursion.**  Removing the maximum `m` drops the Schur count by exactly the
number of ordered representations of `m` as a sum of two elements of `S`. -/
theorem schurCount_erase_max {S : Finset ℕ} (h0 : 0 ∉ S) {m : ℕ} (hmS : m ∈ S)
    (hmax : ∀ x ∈ S, x ≤ m) :
    E3 S = E3 (S.erase m) + repCount S m := by
  -- The Schur pairs of `S` split as {sum = m}  ⊔  {a Schur pair of S.erase m}.
  have hdisj : Disjoint ((S ×ˢ S).filter (fun p => p.1 + p.2 = m))
      ((S.erase m ×ˢ S.erase m).filter (fun p => p.1 + p.2 ∈ S.erase m)) := by
    rw [Finset.disjoint_left]
    intro q hqA hqB
    simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_erase] at hqA hqB
    exact hqB.2.1 hqA.2
  have hunion : (S ×ˢ S).filter (fun p => p.1 + p.2 = m) ∪
      (S.erase m ×ˢ S.erase m).filter (fun p => p.1 + p.2 ∈ S.erase m)
      = (S ×ˢ S).filter (fun p => p.1 + p.2 ∈ S) := by
    ext q
    simp only [Finset.mem_union, Finset.mem_filter, Finset.mem_product, Finset.mem_erase]
    constructor
    · rintro (⟨⟨ha, hb⟩, hm⟩ | ⟨⟨⟨_, ha⟩, ⟨_, hb⟩⟩, ⟨_, hab⟩⟩)
      · exact ⟨⟨ha, hb⟩, by rw [hm]; exact hmS⟩
      · exact ⟨⟨ha, hb⟩, hab⟩
    · rintro ⟨⟨ha, hb⟩, hab⟩
      by_cases hm : q.1 + q.2 = m
      · exact Or.inl ⟨⟨ha, hb⟩, hm⟩
      · obtain ⟨hane, hbne⟩ := ne_max_of_sum_mem h0 hmax ha hb hab
        exact Or.inr ⟨⟨⟨hane, ha⟩, ⟨hbne, hb⟩⟩, ⟨hm, hab⟩⟩
  have hcard : E3 S = ((S ×ˢ S).filter (fun p => p.1 + p.2 = m)).card
      + ((S.erase m ×ˢ S.erase m).filter (fun p => p.1 + p.2 ∈ S.erase m)).card := by
    rw [show E3 S = ((S ×ˢ S).filter (fun p => p.1 + p.2 ∈ S)).card from rfl, ← hunion,
      Finset.card_union_of_disjoint hdisj]
  rw [hcard]
  simp only [E3, repCount]
  omega

/-- **Rung B — the representation bound.**  The maximum has at most `|S|-1` ordered representations,
because `(a,b) ↦ a` injects `{(a,b) : a+b=m}` into `S.erase m` (each `a < m`, and `a` determines `b`). -/
theorem repCount_le {S : Finset ℕ} (h0 : 0 ∉ S) {m : ℕ} (hmS : m ∈ S) :
    repCount S m ≤ S.card - 1 := by
  rw [← Finset.card_erase_of_mem hmS, repCount]
  refine Finset.card_le_card_of_injOn (fun p => p.1) ?_ ?_
  · rintro ⟨a, b⟩ hp
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp
    obtain ⟨⟨ha, hb⟩, hab⟩ := hp
    have hb0 : 0 < b := Nat.pos_of_ne_zero (fun h => h0 (h ▸ hb))
    simp only [Finset.mem_coe, Finset.mem_erase]
    exact ⟨by omega, ha⟩
  · rintro ⟨a, b⟩ hp ⟨c, d⟩ hq hpq
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp hq
    obtain ⟨_, hab⟩ := hp; obtain ⟨_, hcd⟩ := hq
    have hac : a = c := hpq
    simp only [Prod.mk.injEq]
    exact ⟨hac, by omega⟩

/-- Pascal step `C(n+1,2) = C(n,2) + n`, robust against numeral matching. -/
private lemma choose_two_succ (n : ℕ) : (n + 1).choose 2 = n.choose 2 + n := by
  show (n + 1).choose (1 + 1) = n.choose (1 + 1) + n
  rw [Nat.choose_succ_succ, Nat.choose_one_right, Nat.succ_eq_add_one]
  omega

/-- **The deficit peel recursion.**  `deficit S + repCount S m = deficit (S.erase m) + (|S|-1)`:
removing the maximum lowers the number of 2-subsets by `|S|-1` and the realised ones (Schur pairs) by
`repCount S m`, so the deficit changes by the **peel cost** `(|S|-1) − repCount S m`. -/
theorem deficit_erase_max {S : Finset ℕ} (h0 : 0 ∉ S) {m : ℕ} (hmS : m ∈ S)
    (hmax : ∀ x ∈ S, x ≤ m) :
    deficit S + repCount S m = deficit (S.erase m) + (S.card - 1) := by
  have hE3 := schurCount_erase_max h0 hmS hmax
  have hle := E3_le_choose S h0
  have hle' := E3_le_choose (S.erase m) (fun h => h0 (Finset.mem_of_mem_erase h))
  have hcard : (S.erase m).card = S.card - 1 := Finset.card_erase_of_mem hmS
  have hk : 1 ≤ S.card := Finset.card_pos.mpr ⟨m, hmS⟩
  have hpascal : S.card.choose 2 = (S.card - 1).choose 2 + (S.card - 1) := by
    obtain ⟨n, hn⟩ : ∃ n, S.card = n + 1 := ⟨S.card - 1, by omega⟩
    rw [hn, Nat.add_sub_cancel, choose_two_succ]
  unfold deficit
  rw [hcard] at hle' ⊢
  omega

/-- **Deficit monotonicity under peeling:** removing the maximum never increases the deficit. -/
theorem deficit_erase_le {S : Finset ℕ} (h0 : 0 ∉ S) {m : ℕ} (hmS : m ∈ S)
    (hmax : ∀ x ∈ S, x ≤ m) : deficit (S.erase m) ≤ deficit S := by
  have h := deficit_erase_max h0 hmS hmax
  have hr := repCount_le h0 hmS
  omega

/-- **The top peel cost is bounded by the total deficit:** `(|S|-1) − repCount S m ≤ deficit S`.
The whole ladder's rung structure — each peel's deficiency is at most the accumulated deficit. -/
theorem peelCost_le_deficit {S : Finset ℕ} (h0 : 0 ∉ S) {m : ℕ} (hmS : m ∈ S)
    (hmax : ∀ x ∈ S, x ≤ m) :
    (S.card - 1) - repCount S m ≤ deficit S := by
  have h := deficit_erase_max h0 hmS hmax
  omega

/-- **The base rung — `deficit = 0 ⟺ dilated interval`.**  Zero deficit (every 2-subset realised)
is exactly the dilated-interval rigidity. -/
theorem deficit_eq_zero_iff_dilated (S : Finset ℕ) (h0 : 0 ∉ S) (hne : S.Nonempty) :
    deficit S = 0 ↔ DilatedInterval S := by
  rw [← schurCount_eq_choose_iff_dilated S h0 hne]
  unfold deficit
  have := E3_le_choose S h0
  omega

/-! ### Convenience wrappers anchored at the actual maximum `S.max' hne`. -/

theorem schurCount_erase_max' {S : Finset ℕ} (h0 : 0 ∉ S) (hne : S.Nonempty) :
    E3 S = E3 (S.erase (S.max' hne)) + repCount S (S.max' hne) :=
  schurCount_erase_max h0 (S.max'_mem hne) (fun x hx => S.le_max' x hx)

theorem repCount_le' {S : Finset ℕ} (h0 : 0 ∉ S) (hne : S.Nonempty) :
    repCount S (S.max' hne) ≤ S.card - 1 :=
  repCount_le h0 (S.max'_mem hne)

theorem deficit_erase_le' {S : Finset ℕ} (h0 : 0 ∉ S) (hne : S.Nonempty) :
    deficit (S.erase (S.max' hne)) ≤ deficit S :=
  deficit_erase_le h0 (S.max'_mem hne) (fun x hx => S.le_max' x hx)

end LRCSchurRigidity
