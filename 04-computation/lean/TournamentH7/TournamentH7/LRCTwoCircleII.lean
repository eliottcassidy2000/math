/-
  TournamentH7.LRCTwoCircleII — THE TWO-CIRCLE DEEP CERTIFICATE, PART II
  (death-star-2026-07-17-S53, HYP-7270): deep ⟹ circles, completing the iff.

  THE ARCHITECTURE: case analysis on the smallest failing speed `a`.
  * a = 1: the hub lock nests all failures; some failing speed is ≥ 6, so the
    integer circle has width `q/84` (`deep_smallest_one`).
  * a = 2: `w₂ = 1` is forced (divisor descent + speed-1-safe); parity locks
    kill odd speeds ≤ 11; the 13-branch ray estimate `13·|2p−q| > 6q/7`
    kills coexistence with any even ≥ 4, so 13 is out; the failing set is the
    six evens and speed 12 delivers the half circle (`deep_smallest_two`).
  * 3 ≤ a ≤ 8: THE COLLAPSE — every failing partner `l` satisfies the branch
    congruence `(a/g) ∣ ((l/g)·w_a − ε)` with `ε = 0` forced at reduced
    weight ≤ 14 (the lock) and `ε ∈ {−1,0,1}` otherwise (the branch bound);
    the number of compatible partners is ≤ 4 for EVERY coprime witness — one
    kernel `decide` over `a ∈ [3,8], u ∈ (0,a)` (`compat_card_le`).  Hence
    `|F| ≤ 5 < 6` (`deep_smallest_middle`).
  * a ≥ 9: at most 5 speeds remain (`deep_smallest_large`).

  Assembly: `deep_implies_circles`, and with Part I, `deep_iff_circles` —
  the deep set of the canonical family is EXACTLY the two resonance circles,
  in-kernel.  (Recon: exact over 1185 moduli, every multiplier.)

  Kernel-pure: no `sorry`, no `native_decide` (`decide` on the finite
  congruence table only).
-/
import Mathlib
import TournamentH7.LRCTwoCircle

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Witness extraction from a band failure. -/
theorem witness_of_not_inBand (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13)
    (hq : 0 < q) (h : ¬ inBand v q p i) :
    ∃ w : ℤ, 14 * |v i * (p : ℤ) - w * q| < q := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  unfold inBand at h
  push Not at h
  set r : ℤ := (v i * (p : ℤ)) % (q : ℤ) with hrdef
  set m : ℤ := (v i * (p : ℤ)) / (q : ℤ) with hmdef
  have hdm : (q : ℤ) * m + r = v i * (p : ℤ) := Int.ediv_add_emod _ _
  have hr0 : (0 : ℤ) ≤ r := Int.emod_nonneg _ (by omega)
  have hrq : r < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
  rcases lt_or_ge (14 * r) (q : ℤ) with hlow | hhigh
  · -- low residue: witness m
    refine ⟨m, ?_⟩
    have hval : v i * (p : ℤ) - m * q = r := by linarith
    rw [hval, abs_of_nonneg hr0]
    exact hlow
  · -- high residue: witness m + 1
    have hup : 13 * (q : ℤ) < 14 * r := h hhigh
    refine ⟨m + 1, ?_⟩
    have hval : v i * (p : ℤ) - (m + 1) * q = r - q := by linarith
    rw [hval, abs_of_neg (by linarith : r - (q : ℤ) < 0)]
    linarith

/-- Witness range: a failing speed `a ≥ 1` at `0 < p < q` has `0 ≤ w ≤ a`. -/
theorem witness_range (a w : ℤ) (q p : ℕ) (ha : 1 ≤ a)
    (hp : 0 < p) (hpq : p < q)
    (h : 14 * |a * (p : ℤ) - w * q| < q) :
    0 ≤ w ∧ w ≤ a := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by
    have : (p : ℤ) < q := by exact_mod_cast hpq
    have : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp
    linarith [hpq]
  have hpZ : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp
  have hpqZ : (p : ℤ) < (q : ℤ) := by exact_mod_cast hpq
  have habs := abs_lt.mp (by nlinarith [h] : |a * (p : ℤ) - w * q| < q)
  constructor
  · -- w ≥ 0: else w ≤ −1: a·p − w·q ≥ a·p + q > q
    by_contra hcon
    push Not at hcon
    have hw1 : w ≤ -1 := by omega
    nlinarith [habs.2, hpZ, ha]
  · -- w ≤ a: else w ≥ a+1: w·q − a·p ≥ (a+1)q − a·p > q·... ≥ q
    by_contra hcon
    push Not at hcon
    have hw1 : a + 1 ≤ w := by omega
    nlinarith [habs.1, hpqZ, ha, hpZ]

/-- Lock form (reduced weight ≤ 14): `a' ∣ l'·w_a`. -/
theorem partner_lock (g : ℤ) (a' l' : ℕ) (wa wl : ℤ) (q p : ℕ)
    (ha' : 1 ≤ a') (hl' : 1 ≤ l') (hsum : a' + l' ≤ 14)
    (hA : 14 * |(g * a') * (p : ℤ) - wa * q| < q)
    (hL : 14 * |(g * l') * (p : ℤ) - wl * q| < q) :
    (a' : ℤ) ∣ (l' : ℤ) * wa := by
  have hlock : (l' : ℤ) * wa = (a' : ℤ) * wl :=
    rational_lock_weight14 (g * a') (g * l') wa wl a' l' q p ha' hl' hsum
      (by ring) hA hL
  exact ⟨wl, hlock⟩

/-- Branch form (any reduced weight ≤ 27): `a' ∣ (l'·w_a − ε)` for some
`|ε| ≤ 1`. -/
theorem partner_branch (g : ℤ) (a' l' : ℕ) (wa wl : ℤ) (q p : ℕ)
    (ha' : 1 ≤ a') (hl' : 1 ≤ l') (hsum : a' + l' ≤ 27) (hq : 0 < q)
    (hA : 14 * |(g * a') * (p : ℤ) - wa * q| < q)
    (hL : 14 * |(g * l') * (p : ℤ) - wl * q| < q) :
    ∃ ε : ℤ, |ε| ≤ 1 ∧ (a' : ℤ) ∣ ((l' : ℤ) * wa - ε) := by
  have hbr := rational_branch_bound (g * a') (g * l') wa wl a' l' q p
    ha' hl' hsum hq (by ring) hA hL
  refine ⟨(l' : ℤ) * wa - (a' : ℤ) * wl, hbr, ⟨wl, by ring⟩⟩

/-- The compatibility predicate: partner `l` of the smallest failing speed
`a` with coprime witness residue `u` is arithmetically admissible. -/
def compat (a l u : ℕ) : Prop :=
  let g := Nat.gcd a l
  let a' := a / g
  let l' := l / g
  if a' + l' ≤ 14 then a' ∣ l' * u
  else (a' ∣ l' * u) ∨ (a' ∣ l' * u + 1) ∨ (a' ∣ l' * u + (a' - 1))

instance (a l u : ℕ) : Decidable (compat a l u) := by
  unfold compat
  exact inferInstance

/-- **THE COLLAPSE TABLE** (kernel decide): for every smallest speed
`a ∈ [3,8]` and coprime witness residue `u ∈ (0,a)`, at most 4 partners in
`[a+1, 13]` are compatible. -/
theorem compat_card_le :
    ∀ a ∈ Finset.Icc 3 8, ∀ u ∈ Finset.Ico 1 a, Nat.gcd u a = 1 →
      ((Finset.Icc (a+1) 13).filter fun l => compat a l u).card ≤ 4 := by
  decide

/-- Local abbreviation: speed `k` fails at `p/q` (witness form). -/
def Fails (q p k : ℕ) : Prop :=
  ∃ w : ℤ, 14 * |(k : ℤ) * (p : ℤ) - w * q| < q

/-- Minimality forces the witness into `(0, a)` and coprime to `a`. -/
theorem witness_normalized (q p a : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (ha : 2 ≤ a) (wa : ℤ)
    (hA : 14 * |(a : ℤ) * (p : ℤ) - wa * q| < q)
    (hsafe : ∀ b : ℕ, 1 ≤ b → b < a → ¬ Fails q p b) :
    ∃ u : ℕ, (u : ℤ) = wa ∧ 1 ≤ u ∧ u < a ∧ Nat.gcd u a = 1 := by
  have haZ : (1 : ℤ) ≤ (a : ℤ) := by exact_mod_cast Nat.one_le_of_lt ha
  obtain ⟨hw0, hwa⟩ := witness_range (a : ℤ) wa q p haZ hp hpq hA
  -- exclude wa = 0 (descend to speed 1, witness 0)
  have hne0 : wa ≠ 0 := by
    intro h0
    apply hsafe 1 le_rfl (by omega)
    refine ⟨0, ?_⟩
    have := divisor_descent (a : ℤ) wa (a : ℤ) q p (by linarith)
      dvd_rfl (h0 ▸ dvd_zero _) hA
    have hdiv : (a : ℤ) / a = 1 := Int.ediv_self (by linarith)
    rwa [hdiv, h0, Int.zero_ediv] at this
  -- exclude wa = a (descend to speed 1, witness 1)
  have hnea : wa ≠ (a : ℤ) := by
    intro h0
    apply hsafe 1 le_rfl (by omega)
    refine ⟨1, ?_⟩
    have := divisor_descent (a : ℤ) wa (a : ℤ) q p (by linarith)
      dvd_rfl (h0 ▸ dvd_rfl) hA
    have hdiv : (a : ℤ) / a = 1 := Int.ediv_self (by linarith)
    rwa [hdiv, h0, hdiv] at this
  -- exclude gcd(wa, a) = d > 1 (descend to speed a/d < a)
  set u : ℕ := wa.toNat with hu_def
  have huZ : (u : ℤ) = wa := Int.toNat_of_nonneg hw0
  have hu1 : 1 ≤ u := by omega
  have hua : u < a := by omega
  refine ⟨u, huZ, hu1, hua, ?_⟩
  by_contra hcon
  set d : ℕ := Nat.gcd u a with hd_def
  have hd1 : 1 < d := by
    have hdpos : 0 < d := Nat.gcd_pos_of_pos_left a (by omega)
    omega
  have hda : d ∣ a := Nat.gcd_dvd_right u a
  have hdu : d ∣ u := Nat.gcd_dvd_left u a
  have hd_le : d ≤ a := Nat.le_of_dvd (by omega) hda
  apply hsafe (a / d) (Nat.one_le_div_iff (by omega) |>.mpr hd_le)
    (Nat.div_lt_self (by omega) hd1)
  have hdivcast : ∀ (x y : ℕ), y ∣ x → 0 < y →
      ((x : ℕ) : ℤ) / (y : ℤ) = ((x / y : ℕ) : ℤ) := by
    rintro x y ⟨c, rfl⟩ hy
    rw [Nat.mul_div_cancel_left c hy]
    push_cast
    exact Int.mul_ediv_cancel_left _ (by exact_mod_cast hy.ne')
  refine ⟨((u / d : ℕ) : ℤ), ?_⟩
  have hdesc := divisor_descent (a : ℤ) wa (d : ℤ) q p
    (by exact_mod_cast Nat.pos_of_ne_zero (by omega))
    (Int.natCast_dvd_natCast.mpr hda)
    (by rw [← huZ]; exact Int.natCast_dvd_natCast.mpr hdu) hA
  rw [← huZ] at hdesc
  rwa [hdivcast a d hda (by omega), hdivcast u d hdu (by omega)] at hdesc

open Classical in
/-- The middle collapse: smallest failing speed `a ∈ [3,8]` caps the failing
set at 5 — the compat injection plus the kernel table. -/
theorem middle_case_bound (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (a : ℕ) (ha3 : 3 ≤ a) (ha8 : a ≤ 8)
    (hafail : Fails q p a)
    (hsafe : ∀ b : ℕ, 1 ≤ b → b < a → ¬ Fails q p b) :
    ((Finset.Icc 1 13).filter (Fails q p)).card ≤ 5 := by
  obtain ⟨wa, hA⟩ := hafail
  obtain ⟨u, huZ, hu1, hua, hgcd⟩ :=
    witness_normalized q p a hq hp hpq (by omega) wa hA hsafe
  -- every failing partner l > a is compat
  have hcompat : ∀ l : ℕ, a < l → l ≤ 13 → Fails q p l → compat a l u := by
    intro l hal hl13 ⟨wl, hL⟩
    set g : ℕ := Nat.gcd a l with hg_def
    have hgpos : 0 < g := Nat.gcd_pos_of_pos_left l (by omega)
    have hga : g ∣ a := Nat.gcd_dvd_left a l
    have hgl : g ∣ l := Nat.gcd_dvd_right a l
    set a' : ℕ := a / g with ha'_def
    set l' : ℕ := l / g with hl'_def
    have ha'1 : 1 ≤ a' := Nat.one_le_div_iff hgpos |>.mpr (Nat.le_of_dvd (by omega) hga)
    have hl'1 : 1 ≤ l' := Nat.one_le_div_iff hgpos |>.mpr (Nat.le_of_dvd (by omega) hgl)
    have hcop : Nat.Coprime a' l' := Nat.coprime_div_gcd_div_gcd hgpos
    have hsum27 : a' + l' ≤ 27 := by
      have h1 : a' ≤ a := Nat.div_le_self a g
      have h2 : l' ≤ l := Nat.div_le_self l g
      omega
    have haeq : (a : ℤ) = (g : ℤ) * a' := by
      rw [ha'_def]
      exact_mod_cast (Nat.mul_div_cancel' hga).symm
    have hleq : (l : ℤ) = (g : ℤ) * l' := by
      rw [hl'_def]
      exact_mod_cast (Nat.mul_div_cancel' hgl).symm
    have hA' : 14 * |((g : ℤ) * a') * (p : ℤ) - wa * q| < q := by rwa [← haeq]
    have hL' : 14 * |((g : ℤ) * l') * (p : ℤ) - wl * q| < q := by rwa [← hleq]
    unfold compat
    simp only [← hg_def, ← ha'_def, ← hl'_def]
    by_cases hw : a' + l' ≤ 14
    · rw [if_pos hw]
      have hdvd := partner_lock (g : ℤ) a' l' wa wl q p ha'1 hl'1 hw hA' hL'
      rw [← huZ] at hdvd
      exact_mod_cast hdvd
    · rw [if_neg hw]
      obtain ⟨ε, hε, hdvd⟩ := partner_branch (g : ℤ) a' l' wa wl q p
        ha'1 hl'1 hsum27 hq hA' hL'
      rw [← huZ] at hdvd
      have hε3 : ε = -1 ∨ ε = 0 ∨ ε = 1 := by
        rcases abs_le.mp hε with ⟨h1, h2⟩
        omega
      rcases hε3 with rfl | rfl | rfl
      · right; left
        have hbr : (a' : ℤ) ∣ ((l' * u : ℕ) : ℤ) + 1 := by
          have heq : ((l' * u : ℕ) : ℤ) + 1 = (l' : ℤ) * (u : ℤ) - (-1) := by
            push_cast; ring
          rw [heq]
          exact hdvd
        exact_mod_cast hbr
      · left
        have hbr : (a' : ℤ) ∣ ((l' * u : ℕ) : ℤ) := by
          have heq : ((l' * u : ℕ) : ℤ) = (l' : ℤ) * (u : ℤ) - 0 := by
            push_cast; ring
          rw [heq]
          exact hdvd
        exact_mod_cast hbr
      · right; right
        have hlu1 : 1 ≤ l' * u := Nat.one_le_iff_ne_zero.mpr (by positivity)
        have hstep : (a' : ℤ) ∣ ((l' * u : ℕ) : ℤ) - 1 := by
          have heq : ((l' * u : ℕ) : ℤ) - 1 = (l' : ℤ) * (u : ℤ) - 1 := by
            push_cast; ring
          rw [heq]
          exact hdvd
        have hbr : (a' : ℤ) ∣ ((l' * u + (a' - 1) : ℕ) : ℤ) := by
          have hcast : ((l' * u + (a' - 1) : ℕ) : ℤ)
              = (((l' * u : ℕ) : ℤ) - 1) + a' := by
            have ha'ge : 1 ≤ a' := ha'1
            push_cast [Nat.cast_sub ha'ge]
            omega
          rw [hcast]
          exact dvd_add hstep dvd_rfl
        exact_mod_cast hbr
  -- inject the failing set into {a} ∪ compat-partners
  have hsub : (Finset.Icc 1 13).filter (Fails q p)
      ⊆ insert a ((Finset.Icc (a+1) 13).filter fun l => compat a l u) := by
    intro l hl
    rw [Finset.mem_filter, Finset.mem_Icc] at hl
    obtain ⟨⟨hl1, hl13⟩, hlfail⟩ := hl
    rcases lt_trichotomy l a with hlt | rfl | hgt
    · exact absurd hlfail (hsafe l hl1 hlt)
    · exact Finset.mem_insert_self _ _
    · apply Finset.mem_insert_of_mem
      rw [Finset.mem_filter, Finset.mem_Icc]
      exact ⟨⟨by omega, hl13⟩, hcompat l hgt hl13 hlfail⟩
  calc ((Finset.Icc 1 13).filter (Fails q p)).card
      ≤ (insert a ((Finset.Icc (a+1) 13).filter fun l => compat a l u)).card :=
        Finset.card_le_card hsub
    _ ≤ ((Finset.Icc (a+1) 13).filter fun l => compat a l u).card + 1 :=
        Finset.card_insert_le _ _
    _ ≤ 4 + 1 := by
        have := compat_card_le a (Finset.mem_Icc.mpr ⟨ha3, ha8⟩) u
          (Finset.mem_Ico.mpr ⟨hu1, hua⟩) hgcd
        omega
    _ = 5 := rfl

open Classical in
/-- The large case: all speeds ≤ 8 safe caps the failing set at 5. -/
theorem large_case_bound (q p : ℕ)
    (hsafe : ∀ b : ℕ, 1 ≤ b → b ≤ 8 → ¬ Fails q p b) :
    ((Finset.Icc 1 13).filter (Fails q p)).card ≤ 5 := by
  have hsub : (Finset.Icc 1 13).filter (Fails q p) ⊆ Finset.Icc 9 13 := by
    intro l hl
    rw [Finset.mem_filter, Finset.mem_Icc] at hl
    obtain ⟨⟨h1, h13⟩, hf⟩ := hl
    rw [Finset.mem_Icc]
    refine ⟨?_, h13⟩
    by_contra hcon
    exact hsafe l h1 (by omega) hf
  calc ((Finset.Icc 1 13).filter (Fails q p)).card
      ≤ (Finset.Icc 9 13).card := Finset.card_le_card hsub
    _ = 5 := rfl

open Classical in
/-- The hub case: speed 1 fails and six speeds fail ⟹ the integer circle at
width `q/84`. -/
theorem hub_case_circleI (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (h1fail : Fails q p 1)
    (hsix : 6 ≤ ((Finset.Icc 1 13).filter (Fails q p)).card) :
    84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q := by
  obtain ⟨s, hS⟩ := h1fail
  have hs01 : s = 0 ∨ s = 1 := by
    obtain ⟨h0, h1⟩ := witness_range 1 s q p le_rfl hp hpq (by simpa using hS)
    omega
  -- some failing speed is ≥ 6
  have hbig : ∃ l : ℕ, 6 ≤ l ∧ l ≤ 13 ∧ Fails q p l := by
    by_contra hcon
    push Not at hcon
    have hsub : (Finset.Icc 1 13).filter (Fails q p) ⊆ Finset.Icc 1 5 := by
      intro l hl
      rw [Finset.mem_filter, Finset.mem_Icc] at hl
      obtain ⟨⟨ha, hb⟩, hf⟩ := hl
      rw [Finset.mem_Icc]
      refine ⟨ha, ?_⟩
      by_contra hc
      exact absurd hf (by
        have := hcon l (by omega) hb
        exact this)
    have := Finset.card_le_card hsub
    simp [Nat.card_Icc] at this
    omega
  obtain ⟨l, hl6, hl13, wl, hL⟩ := hbig
  -- the hub lock pins w_l = l·s
  have hlock : (l : ℤ) * s = (1 : ℤ) * wl := by
    have := rational_lock_weight14 1 (l : ℤ) s wl 1 l q p le_rfl
      (by omega) (by omega) (by push_cast; ring)
      (by simpa using hS) (by simpa using hL)
    simpa using this
  -- l·|p − s·q| < q/14 with l ≥ 6 gives the 84-width
  have hXl : (l : ℤ) * (p : ℤ) - wl * q = (l : ℤ) * ((p : ℤ) - s * q) := by
    have hwl : wl = (l : ℤ) * s := by linarith [hlock]
    rw [hwl]; ring
  rw [hXl, abs_mul, abs_of_pos (by exact_mod_cast (by omega : 0 < l))] at hL
  have hl6Z : (6 : ℤ) ≤ (l : ℤ) := by exact_mod_cast hl6
  have h84 : 84 * |(p : ℤ) - s * q| < q := by
    nlinarith [hL, abs_nonneg ((p : ℤ) - s * q)]
  rcases hs01 with rfl | rfl
  · left
    have : |(p : ℤ) - 0 * q| = p := by
      simp [abs_of_pos (by exact_mod_cast hp : (0:ℤ) < (p:ℤ))]
    rw [this] at h84
    exact h84
  · right
    have : |(p : ℤ) - 1 * q| = (q : ℤ) - p := by
      rw [abs_of_neg (by
        have : (p : ℤ) < q := by exact_mod_cast hpq
        linarith : (p : ℤ) - 1 * q < 0)]
      ring
    rw [this] at h84
    exact h84

open Classical in
/-- The even case: speed 2 fails, speed 1 safe, six speeds fail ⟹ the half
circle at width `q/84`. -/
theorem even_case_circleII (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (h2fail : Fails q p 2) (h1safe : ¬ Fails q p 1)
    (hsix : 6 ≤ ((Finset.Icc 1 13).filter (Fails q p)).card) :
    84 * |2 * (p : ℤ) - q| < q := by
  obtain ⟨w2, hW⟩ := h2fail
  have hsafe : ∀ b : ℕ, 1 ≤ b → b < 2 → ¬ Fails q p b := by
    intro b hb1 hb2
    have hb : b = 1 := by omega
    rw [hb]
    exact h1safe
  obtain ⟨u, huZ, hu1, hu2, _⟩ :=
    witness_normalized q p 2 hq hp hpq le_rfl w2 hW hsafe
  have hu : u = 1 := by omega
  have hw2 : w2 = 1 := by rw [hu] at huZ; exact_mod_cast huZ.symm
  subst hw2
  -- odd speeds 3..11 cannot fail (parity lock)
  have hodd : ∀ l : ℕ, 3 ≤ l → l ≤ 11 → l % 2 = 1 → ¬ Fails q p l := by
    rintro l h3 h11 hoddl ⟨wl, hL⟩
    have hlock := rational_lock_weight14 2 (l : ℤ) 1 wl 2 l q p
      (by omega) (by omega) (by omega)
      (by push_cast; ring) (by simpa using hW) (by simpa using hL)
    have hpar : (l : ℤ) = 2 * wl := by
      push_cast at hlock
      linarith [hlock]
    omega
  -- even speeds 2m (2 ≤ m ≤ 6) fail with witness m, giving m·|2p−q| bounds
  have heven_bound : ∀ m : ℕ, 1 ≤ m → m ≤ 6 → Fails q p (2 * m) →
      14 * (m : ℤ) * |2 * (p : ℤ) - q| < q := by
    rintro m hm1 hm6 ⟨wm, hL⟩
    have hlock := rational_lock_weight14 2 (2 * m : ℤ) 1 wm 1 m q p
      le_rfl (by omega) (by omega)
      (by push_cast; ring) (by simpa using hW)
      (by
        have hcast : ((2 * m : ℕ) : ℤ) = 2 * (m : ℤ) := by push_cast; ring
        rw [← hcast]
        simpa using hL)
    have hwm : wm = (m : ℤ) := by
      have : (m : ℤ) * 1 = 1 * wm := hlock
      linarith
    subst hwm
    have hX : ((2 * m : ℕ) : ℤ) * (p : ℤ) - (m : ℤ) * q
        = (m : ℤ) * (2 * (p : ℤ) - q) := by push_cast; ring
    rw [hX, abs_mul, abs_of_pos (by exact_mod_cast hm1 : (0:ℤ) < (m:ℤ))] at hL
    calc 14 * (m : ℤ) * |2 * (p : ℤ) - q|
        = 14 * ((m : ℤ) * |2 * (p : ℤ) - q|) := by ring
      _ < q := hL
  by_cases h13 : Fails q p 13
  · -- 13 alongside 2: the ray estimate kills every even ≥ 4, |F| ≤ 2, absurd
    exfalso
    obtain ⟨w13, hL13⟩ := h13
    have hbr := rational_branch_bound 2 13 1 w13 2 13 q p
      (by omega) (by omega) (by omega) hq
      (by push_cast) (by simpa using hW) (by simpa using hL13)
    have hw13 : w13 = 6 ∨ w13 = 7 := by
      have habs := abs_le.mp hbr
      push_cast at habs
      omega
    have hL13' : 14 * |13 * (p : ℤ) - w13 * q| < q := by simpa using hL13
    have hray : 6 * (q : ℤ) < 91 * |2 * (p : ℤ) - q| := by
      rcases hw13 with rfl | rfl
      · rcases abs_cases (2 * (p : ℤ) - q) with ⟨h1, _⟩ | ⟨h1, _⟩ <;>
          rcases abs_cases (13 * (p : ℤ) - 6 * q) with ⟨h2, _⟩ | ⟨h2, _⟩ <;>
          rw [h1] <;> rw [h2] at hL13' <;> linarith [hL13']
      · rcases abs_cases (2 * (p : ℤ) - q) with ⟨h1, _⟩ | ⟨h1, _⟩ <;>
          rcases abs_cases (13 * (p : ℤ) - 7 * q) with ⟨h2, _⟩ | ⟨h2, _⟩ <;>
          rw [h1] <;> rw [h2] at hL13' <;> linarith [hL13']
    have hnoeven : ∀ m : ℕ, 2 ≤ m → m ≤ 6 → ¬ Fails q p (2 * m) := by
      intro m hm2 hm6 hf
      have hb := heven_bound m (by omega) hm6 hf
      have hm2Z : (2 : ℤ) ≤ (m : ℤ) := by exact_mod_cast hm2
      nlinarith [hb, hray, abs_nonneg (2 * (p : ℤ) - q)]
    have hsub : (Finset.Icc 1 13).filter (Fails q p) ⊆ ({2, 13} : Finset ℕ) := by
      intro l hl
      rw [Finset.mem_filter, Finset.mem_Icc] at hl
      obtain ⟨⟨h1, h13'⟩, hf⟩ := hl
      simp only [Finset.mem_insert, Finset.mem_singleton]
      interval_cases l
      · exact absurd hf h1safe
      · left; rfl
      · exact absurd hf (hodd 3 (by omega) (by omega) rfl)
      · exact absurd hf (hnoeven 2 le_rfl (by omega))
      · exact absurd hf (hodd 5 (by omega) (by omega) rfl)
      · exact absurd hf (hnoeven 3 (by omega) (by omega))
      · exact absurd hf (hodd 7 (by omega) (by omega) rfl)
      · exact absurd hf (hnoeven 4 (by omega) (by omega))
      · exact absurd hf (hodd 9 (by omega) (by omega) rfl)
      · exact absurd hf (hnoeven 5 (by omega) (by omega))
      · exact absurd hf (hodd 11 (by omega) (by omega) rfl)
      · exact absurd hf (hnoeven 6 (by omega) le_rfl)
      · right; rfl
    have := Finset.card_le_card hsub
    have hc2 : ({2, 13} : Finset ℕ).card ≤ 2 := by decide
    omega
  · -- 13 safe: the failing set sits in the six evens; 12 must fail
    have h12 : Fails q p 12 := by
      by_contra h12safe
      have hsub : (Finset.Icc 1 13).filter (Fails q p)
          ⊆ ({2, 4, 6, 8, 10} : Finset ℕ) := by
        intro l hl
        rw [Finset.mem_filter, Finset.mem_Icc] at hl
        obtain ⟨⟨h1, h13'⟩, hf⟩ := hl
        simp only [Finset.mem_insert, Finset.mem_singleton]
        interval_cases l
        · exact absurd hf h1safe
        · left; rfl
        · exact absurd hf (hodd 3 (by omega) (by omega) rfl)
        · right; left; rfl
        · exact absurd hf (hodd 5 (by omega) (by omega) rfl)
        · right; right; left; rfl
        · exact absurd hf (hodd 7 (by omega) (by omega) rfl)
        · right; right; right; left; rfl
        · exact absurd hf (hodd 9 (by omega) (by omega) rfl)
        · right; right; right; right; rfl
        · exact absurd hf (hodd 11 (by omega) (by omega) rfl)
        · exact absurd hf h12safe
        · exact absurd hf h13
      have := Finset.card_le_card hsub
      have hc5 : ({2, 4, 6, 8, 10} : Finset ℕ).card ≤ 5 := by decide
      omega
    obtain ⟨w12, hL12⟩ := h12
    have hlock := rational_lock_weight14 2 12 1 w12 1 6 q p
      le_rfl (by omega) (by omega)
      (by push_cast) (by simpa using hW) (by simpa using hL12)
    have hw12 : w12 = 6 := by
      push_cast at hlock
      linarith [hlock]
    subst hw12
    have hL12' : 14 * |12 * (p : ℤ) - 6 * q| < q := by simpa using hL12
    have hX : (12 : ℤ) * (p : ℤ) - 6 * q = 6 * (2 * (p : ℤ) - q) := by ring
    rw [hX, abs_mul, abs_of_pos (by norm_num : (0:ℤ) < 6)] at hL12'
    linarith [hL12']



open Classical in
/-- Speed–index bridge: speed `i+1` fails iff index `i` is out of band. -/
theorem fails_iff_not_inBand (q p : ℕ) (hq : 0 < q) (i : Fin 13) :
    Fails q p ((i : ℕ) + 1) ↔ ¬ inBand canonical q p i := by
  have hcan : canonical i = (((i : ℕ) + 1 : ℕ) : ℤ) := by
    unfold canonical; push_cast; ring
  constructor
  · rintro ⟨w, hw⟩
    apply not_inBand_of_witness canonical q p i hq w
    rw [hcan]
    exact hw
  · intro h
    obtain ⟨w, hw⟩ := witness_of_not_inBand canonical q p i hq h
    refine ⟨w, ?_⟩
    rw [← hcan]
    exact hw

open Classical in
/-- `bandCount` equals the failing-speed count. -/
theorem bandCount_eq_fails_card (q p : ℕ) (hq : 0 < q) :
    bandCount canonical q p = ((Finset.Icc 1 13).filter (Fails q p)).card := by
  unfold bandCount
  apply Finset.card_bij (fun (i : Fin 13) _ => (i : ℕ) + 1)
  · intro i hi
    rw [Finset.mem_filter] at hi
    rw [Finset.mem_filter, Finset.mem_Icc]
    exact ⟨⟨by omega, by omega⟩, (fails_iff_not_inBand q p hq i).mpr hi.2⟩
  · intro i _ j _ h
    exact Fin.ext (by omega)
  · intro l hl
    rw [Finset.mem_filter, Finset.mem_Icc] at hl
    obtain ⟨⟨h1, h13⟩, hf⟩ := hl
    refine ⟨⟨l - 1, by omega⟩, ?_, by
      show (l - 1) + 1 = l
      omega⟩
    rw [Finset.mem_filter]
    refine ⟨Finset.mem_univ _, ?_⟩
    apply (fails_iff_not_inBand q p hq ⟨l - 1, by omega⟩).mp
    show Fails q p ((l - 1) + 1)
    have hll : (l - 1) + 1 = l := by omega
    rw [hll]
    exact hf

open Classical in
/-- **PART II: deep ⟹ circles.** -/
theorem deep_implies_circles (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hdeep : 6 ≤ bandCount canonical q p) :
    (84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q) ∨ 84 * |2 * (p : ℤ) - q| < q := by
  rw [bandCount_eq_fails_card q p hq] at hdeep
  by_cases h1 : Fails q p 1
  · exact Or.inl (hub_case_circleI q p hq hp hpq h1 hdeep)
  by_cases h2 : Fails q p 2
  · exact Or.inr (even_case_circleII q p hq hp hpq h2 h1 hdeep)
  exfalso
  set S : Finset ℕ := (Finset.Icc 3 8).filter (Fails q p) with hS
  by_cases hSne : S.Nonempty
  · set a : ℕ := S.min' hSne with ha
    have haS : a ∈ S := Finset.min'_mem _ _
    rw [hS, Finset.mem_filter, Finset.mem_Icc] at haS
    obtain ⟨⟨ha3, ha8⟩, hafail⟩ := haS
    have hbound := middle_case_bound q p hq hp hpq a ha3 ha8 hafail (by
      intro b hb1 hba
      rcases Nat.lt_or_ge b 3 with hb3 | hb3
      · interval_cases b
        · exact h1
        · exact h2
      · intro hbf
        have hbS : b ∈ S := by
          rw [hS, Finset.mem_filter, Finset.mem_Icc]
          exact ⟨⟨hb3, by omega⟩, hbf⟩
        have := Finset.min'_le S b hbS
        omega)
    omega
  · have hallsafe : ∀ b : ℕ, 1 ≤ b → b ≤ 8 → ¬ Fails q p b := by
      intro b hb1 hb8 hbf
      rcases Nat.lt_or_ge b 3 with hb3 | hb3
      · interval_cases b
        · exact h1 hbf
        · exact h2 hbf
      · exact hSne ⟨b, by
          rw [hS, Finset.mem_filter, Finset.mem_Icc]
          exact ⟨⟨hb3, hb8⟩, hbf⟩⟩
    have := large_case_bound q p hallsafe
    omega

/-- **THE TWO-CIRCLE THEOREM** (with Part I): the deep set of the canonical
family is EXACTLY the two resonance circles. -/
theorem deep_iff_circles (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q) :
    6 ≤ bandCount canonical q p ↔
      ((84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q) ∨
        84 * |2 * (p : ℤ) - q| < q) := by
  constructor
  · exact deep_implies_circles q p hq hp hpq
  · rintro (hI | hII)
    · exact circleI_deep q p hq hp hpq hI
    · exact circleII_deep q p hq hII

/-! ## Axiom audit (Part II core pieces) -/
#print axioms witness_of_not_inBand
#print axioms witness_range
#print axioms partner_lock
#print axioms partner_branch
#print axioms compat_card_le
#print axioms witness_normalized
#print axioms middle_case_bound
#print axioms large_case_bound
#print axioms hub_case_circleI
#print axioms even_case_circleII
#print axioms bandCount_eq_fails_card
#print axioms deep_implies_circles
#print axioms deep_iff_circles

end LRC14Concrete
end LonelyRunner
