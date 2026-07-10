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
(`deficit_eq_zero_iff_dilated`).  The whole ladder assembles as `deficit S = totalPeelCost S`
(`deficit_eq_totalPeelCost`): the deficit is exactly the sum of the per-peel costs, so a small deficit
means few deficient peels down the max-peeling chain.  Its top (`deficit 0`) is the dilated-interval
rigidity `schurCount_eq_choose_iff_dilated`, and a full peel is a local reflection symmetry
(`repCount_max_eq_iff`).

The *quantitative* capstone `dist_to_dilated S ≤ deficit S` (S can be made a dilated interval by changing
`≤ deficit` elements) is a Freiman-stability statement that is **NOT proved here** — and is in fact FALSE
for `|S| ≤ 4` (e.g. `S = {1,4,5}`: deficit `1`, distance `2`), holding only for `|S| ≥ 5` (verified
exhaustively to `N = 5k`, 0 failures).  That `|S| ≥ 5` threshold coincides exactly with the burden axis's
`LRCFreimanAP.ap_of_min_burden` (false for `|S| ≤ 4`, MISTAKE-133): both stability statements fail on
small sets for the same accidental-additive-structure reason.  For LRC(14) the relevant regime is `k=13`.

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

/-- **Full-peel characterization.**  The maximum `m` is *fully represented* (`repCount = |S|-1`, so the
peel cost is `0`) **iff** `S` below `m` is closed under the reflection `a ↦ m - a`: every `a ∈ S` with
`a ≠ m` has its complement `m - a ∈ S`.  This is the local reflection symmetry whose failure at some peel
is exactly what a positive deficit records. -/
theorem repCount_max_eq_iff {S : Finset ℕ} (h0 : 0 ∉ S) {m : ℕ} (hmS : m ∈ S)
    (hmax : ∀ x ∈ S, x ≤ m) :
    repCount S m = S.card - 1 ↔ ∀ a ∈ S, a ≠ m → m - a ∈ S := by
  set F := (S ×ˢ S).filter (fun p => p.1 + p.2 = m) with hF
  have himg : F.image (fun p => p.1) ⊆ S.erase m := by
    intro x hx
    rw [Finset.mem_image] at hx
    obtain ⟨⟨a, b⟩, hp, rfl⟩ := hx
    rw [hF, Finset.mem_filter, Finset.mem_product] at hp
    obtain ⟨⟨ha, hb⟩, hab⟩ := hp
    have hb0 : 0 < b := Nat.pos_of_ne_zero (fun h => h0 (h ▸ hb))
    exact Finset.mem_erase.mpr ⟨by omega, ha⟩
  have hinj : Set.InjOn (fun p => p.1) (F : Set (ℕ × ℕ)) := by
    rintro ⟨a, b⟩ hp ⟨c, d⟩ hq hpq
    simp only [hF, Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp hq
    obtain ⟨_, hab⟩ := hp; obtain ⟨_, hcd⟩ := hq
    have hac : a = c := hpq
    simp only [Prod.mk.injEq]; exact ⟨hac, by omega⟩
  have hcardimg : (F.image (fun p => p.1)).card = repCount S m := by
    rw [repCount, ← hF, Finset.card_image_of_injOn hinj]
  rw [← Finset.card_erase_of_mem hmS]
  constructor
  · intro hcard a haS hane
    have heq : F.image (fun p => p.1) = S.erase m :=
      Finset.eq_of_subset_of_card_le himg (by rw [hcardimg]; omega)
    have hain : a ∈ F.image (fun p => p.1) := by rw [heq]; exact Finset.mem_erase.mpr ⟨hane, haS⟩
    rw [Finset.mem_image] at hain
    obtain ⟨⟨x, y⟩, hp, hxa⟩ := hain
    rw [hF, Finset.mem_filter, Finset.mem_product] at hp
    obtain ⟨⟨_, hy⟩, hxy⟩ := hp
    have hxa' : x = a := hxa
    have hya : m - a = y := by omega
    rw [hya]; exact hy
  · intro hclosed
    have hsup : S.erase m ⊆ F.image (fun p => p.1) := by
      intro a ha
      rw [Finset.mem_erase] at ha
      obtain ⟨hane, haS⟩ := ha
      have hma : m - a ∈ S := hclosed a haS hane
      have haM : a ≤ m := hmax a haS
      rw [Finset.mem_image]
      refine ⟨(a, m - a), ?_, rfl⟩
      rw [hF, Finset.mem_filter, Finset.mem_product]
      exact ⟨⟨haS, hma⟩, by omega⟩
    rw [← hcardimg, Finset.Subset.antisymm himg hsup]

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

/-- **The total peel cost** — the sum of the per-peel deficiencies `(|T|-1) − repCount T (max T)` down
the maximum-peeling chain `S ⊃ S.erase (max) ⊃ …`. -/
def totalPeelCost (S : Finset ℕ) : ℕ :=
  if h : S.Nonempty then
    (S.card - 1 - repCount S (S.max' h)) + totalPeelCost (S.erase (S.max' h))
  else 0
termination_by S.card
decreasing_by exact Finset.card_erase_lt_of_mem (S.max'_mem h)

/-- **The E3-stability ladder, assembled — `deficit S = totalPeelCost S`.**  The deficit is exactly the
accumulated peel cost down the maximum-peeling chain.  Each rung's cost `(|T|-1) − repCount T (max)` is
`≥ 0` (`repCount_le`) and `= 0` iff that peel is *full* (the max is fully represented); so `deficit = 0
⟺` every peel is full `⟺` dilated, and — the quantitative content of the ladder — a small deficit means
few deficient peels down the chain.  This is the E3-axis analogue of accumulating the burden along
opus's `restrictedSum` chain. -/
theorem deficit_eq_totalPeelCost : ∀ (S : Finset ℕ), 0 ∉ S → deficit S = totalPeelCost S := by
  intro S
  induction S using Finset.strongInductionOn with
  | _ S ih =>
    intro h0
    rw [totalPeelCost]
    split
    · rename_i h
      have hmem := S.max'_mem h
      have hmax : ∀ x ∈ S, x ≤ S.max' h := fun x hx => S.le_max' x hx
      have hrec := ih (S.erase (S.max' h)) (Finset.erase_ssubset hmem)
        (fun hh => h0 (Finset.mem_of_mem_erase hh))
      have hpeel := deficit_erase_max h0 hmem hmax
      have hrep := repCount_le h0 hmem
      rw [← hrec]; omega
    · rename_i h
      rw [Finset.not_nonempty_iff_eq_empty] at h
      subst h
      simp [deficit, E3]

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

/-! ### The E₃ diagonal split — separating the endgame's `W₀`-carrier (doublings) from the free part.

By THM-682(d) (monad), the sole carriers of the exact-load `W₀` in the LRC(14) final rung are the
**doublings** `v_b = 2 v_a`; the genuine Schur triples are nearly weightless (line weight `0.0027`).  In
the additive-relation ↔ cycle-length ladder (THM-446) a Schur triple `a+b=c` is a **triangle** (the
smallest odd cycle, the OCF atom) and a doubling `d ↦ 2d` is the **dyadic** rung.  The doubling content
is exactly the *diagonal* `(a,a)` of `E₃`, so the endgame's operative additive content lives on this
diagonal — the 2-adic axis where the dyadic dispatches (LEM-019/021) and the 2-adic tower operate. -/

/-- **The doubling count** — the diagonal of `E₃`: the number of `a ∈ S` with `2a ∈ S`.  A doubling is a
step up a geometric-ratio-2 chain; `doublingCount` counts the edges of the doubling forest. -/
def doublingCount (S : Finset ℕ) : ℕ := (S.filter (fun a => 2 * a ∈ S)).card

/-- **The E₃ diagonal split.**  `E₃ S = doublingCount S + (off-diagonal Schur count)`: the Schur count
splits into the diagonal doublings `(a,a)` (`2a ∈ S`, the dyadic/`W₀`-carrier content) and the
off-diagonal genuine Schur triples `(a,b)`, `a ≠ b`, `a + b ∈ S` (the triangle/3-cycle content, nearly
weightless in the final rung). -/
theorem schurCount_eq_doubling_add_offDiag (S : Finset ℕ) :
    E3 S = doublingCount S
      + ((S ×ˢ S).filter (fun p => p.1 + p.2 ∈ S ∧ p.1 ≠ p.2)).card := by
  have hdiag : (((S ×ˢ S).filter (fun p => p.1 + p.2 ∈ S)).filter (fun p => p.1 = p.2))
             = (S.filter (fun a => 2 * a ∈ S)).image (fun a => (a, a)) := by
    ext ⟨x, y⟩
    simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_image, Prod.mk.injEq]
    constructor
    · rintro ⟨⟨⟨hx, _⟩, hxy⟩, heq⟩
      subst heq
      exact ⟨x, ⟨hx, by rwa [show x + x = 2 * x from by ring] at hxy⟩, rfl, rfl⟩
    · rintro ⟨a, ⟨haS, h2a⟩, rfl, rfl⟩
      exact ⟨⟨⟨haS, haS⟩, by rwa [show a + a = 2 * a from by ring]⟩, rfl⟩
  have hsplit := Finset.filter_card_add_filter_neg_card_eq_card
    (s := (S ×ˢ S).filter (fun p => p.1 + p.2 ∈ S)) (p := fun p => p.1 = p.2)
  rw [hdiag, Finset.card_image_of_injOn (fun a _ b _ h => congrArg Prod.fst h),
    Finset.filter_filter] at hsplit
  exact hsplit.symm

end LRCSchurRigidity
