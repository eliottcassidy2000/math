import TournamentH7.LRCSparseBranch

/-!
# The exact one-dimensional lattice behind the positive sparse branch

`LRCSparseBranch.branch_one_iff` reduces the positive Bezout branch to two
strict bands in one integer `Z`.  This file performs the remaining local
intersection and count.  For `0 < i <= j`, `j - i <= 13`, and
`j * X - i * Y = 1`, the first band's lower wall and the second band's upper
wall are the active walls.  Consequently the admissible integers form one
explicit `Finset.Ioo`, whose cardinality is the advertised difference of
Euclidean floors (with `Int.toNat`, which also handles the zero-width
`i + j = 14` boundary correctly).

The residue-bijection step is also carried out: the positive branch has the
exact truncated floor count, the negative branch has the same count by
reflection, and the three disjoint Bezout branches give the complete finite-
`q` joint-failure formula for every canonically ordered reduced pair.

Tournament-analysis audit: the faithful vertices here are Bezout branches
and integer lattice points.  Their observable is wall membership, and the
active-wall order is transitive, so a tournament orientation adds no data.
Quotienting to runners would forget the branch residue and the interval
intercept.  The challenged assumption is that both bands contribute four
independent walls: Bezout order makes two of them redundant.

No `sorry`; no `native_decide`.
-/

namespace LonelyRunner
namespace LRC14Concrete

open Finset

noncomputable section

/-- The integer interval cut out by the active lower and upper walls of the
positive Bezout branch. -/
def sparsePositiveLattice (i j X Y q : ℤ) : Finset ℤ :=
  Finset.Ioo
    (((14 * X - 1) * q) / (14 * i))
    ((((14 * Y + 1) * q - 1) / (14 * j)) + 1)

/-- Under the canonical ordering hypotheses, the two Bezout bands are exactly
one open integer interval. -/
theorem mem_sparsePositiveLattice_iff
    (i j X Y q Z : ℤ)
    (hi : 0 < i) (hij : i ≤ j) (hgap : j - i ≤ 13)
    (hq : 0 < q) (hbez : j * X - i * Y = 1) :
    Z ∈ sparsePositiveLattice i j X Y q ↔
      14 * |i * Z - X * q| < q ∧
      14 * |j * Z - Y * q| < q := by
  have hj : 0 < j := lt_of_lt_of_le hi hij
  have hCi : 0 < 14 * i := by positivity
  have hCj : 0 < 14 * j := by positivity
  have hlowerCoefficient :
      0 < j * (14 * X - 1) - i * (14 * Y - 1) := by
    nlinarith
  have hupperCoefficient :
      0 < j * (14 * X + 1) - i * (14 * Y + 1) := by
    nlinarith
  constructor
  · intro hZ
    rw [sparsePositiveLattice, Finset.mem_Ioo] at hZ
    have hlower : (14 * X - 1) * q < Z * (14 * i) :=
      (Int.ediv_lt_iff_lt_mul hCi).mp hZ.1
    have hupperDiv :
        Z ≤ ((14 * Y + 1) * q - 1) / (14 * j) := by
      omega
    have hupperWeak :
        Z * (14 * j) ≤ (14 * Y + 1) * q - 1 :=
      (Int.le_ediv_iff_mul_le hCj).mp hupperDiv
    have hupper : Z * (14 * j) < (14 * Y + 1) * q := by
      omega

    have hlowerCross :
        i * ((14 * Y - 1) * q) < j * ((14 * X - 1) * q) := by
      nlinarith [mul_pos hlowerCoefficient hq]
    have hlowerScaled :
        j * ((14 * X - 1) * q) < j * (Z * (14 * i)) :=
      mul_lt_mul_of_pos_left hlower hj
    have hlowerSecond : (14 * Y - 1) * q < Z * (14 * j) := by
      nlinarith

    have hupperCross :
        i * ((14 * Y + 1) * q) < j * ((14 * X + 1) * q) := by
      nlinarith [mul_pos hupperCoefficient hq]
    have hupperScaled :
        i * (Z * (14 * j)) < i * ((14 * Y + 1) * q) :=
      mul_lt_mul_of_pos_left hupper hi
    have hupperFirst : Z * (14 * i) < (14 * X + 1) * q := by
      nlinarith

    constructor
    · have habs : |14 * (i * Z - X * q)| < q := by
        rw [abs_lt]
        constructor <;> nlinarith
      calc
        14 * |i * Z - X * q| = |14 * (i * Z - X * q)| := by
          rw [abs_mul]
          norm_num
        _ < q := habs
    · have habs : |14 * (j * Z - Y * q)| < q := by
        rw [abs_lt]
        constructor <;> nlinarith
      calc
        14 * |j * Z - Y * q| = |14 * (j * Z - Y * q)| := by
          rw [abs_mul]
          norm_num
        _ < q := habs
  · rintro ⟨hiBand, hjBand⟩
    have hiAbs : |14 * (i * Z - X * q)| < q := by
      calc
        |14 * (i * Z - X * q)| = 14 * |i * Z - X * q| := by
          rw [abs_mul]
          norm_num
        _ < q := hiBand
    have hjAbs : |14 * (j * Z - Y * q)| < q := by
      calc
        |14 * (j * Z - Y * q)| = 14 * |j * Z - Y * q| := by
          rw [abs_mul]
          norm_num
        _ < q := hjBand
    rw [abs_lt] at hiAbs hjAbs
    have hlower : (14 * X - 1) * q < Z * (14 * i) := by
      nlinarith [hiAbs.1]
    have hupper : Z * (14 * j) ≤ (14 * Y + 1) * q - 1 := by
      nlinarith [hjAbs.2]
    rw [sparsePositiveLattice, Finset.mem_Ioo]
    constructor
    · exact (Int.ediv_lt_iff_lt_mul hCi).mpr hlower
    · have hupperDiv :
          Z ≤ ((14 * Y + 1) * q - 1) / (14 * j) :=
        (Int.le_ediv_iff_mul_le hCj).mpr hupper
      omega

/-- Exact local floor formula for the positive sparse branch lattice.  The
`toNat` is essential at the zero-width boundary; without it, the raw floor
difference can be `-1` when the common wall is itself integral. -/
theorem sparsePositiveLattice_card (i j X Y q : ℤ) :
    (sparsePositiveLattice i j X Y q).card =
      ((((14 * Y + 1) * q - 1) / (14 * j)) -
        (((14 * X - 1) * q) / (14 * i))).toNat := by
  rw [sparsePositiveLattice, Int.card_Ioo]
  congr 1
  ring

/-- Composition with `LRCSparseBranch.branch_one_iff`: positive-branch
multipliers are precisely the residue classes represented by the explicit
lattice interval.  No counting or choice of a canonical representative is
hidden in this statement. -/
theorem branch_one_iff_lattice_modEq
    (g : ℤ) (i' j' : ℕ) (X Y : ℤ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hij : i' ≤ j')
    (hgap : j' - i' ≤ 13) (hq : 0 < q)
    (hcop : Nat.Coprime i' j')
    (hbez : (j' : ℤ) * X - (i' : ℤ) * Y = 1) :
    (Exists fun wa : ℤ => Exists fun wb : ℤ =>
        14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
        14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
        (j' : ℤ) * wa - (i' : ℤ) * wb = 1) ↔
      ∃ Z ∈ sparsePositiveLattice (i' : ℤ) (j' : ℤ) X Y q,
        Z ≡ g * (p : ℤ) [ZMOD (q : ℤ)] := by
  rw [branch_one_iff g i' j' X Y q p hi hj hcop hbez]
  have hiZ : (0 : ℤ) < (i' : ℤ) := by exact_mod_cast hi
  have hijZ : (i' : ℤ) ≤ (j' : ℤ) := by exact_mod_cast hij
  have hgapZ : (j' : ℤ) - (i' : ℤ) ≤ 13 := by
    have hijSub : j' - i' + i' = j' := Nat.sub_add_cancel hij
    exact_mod_cast hgap
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  constructor
  · rintro ⟨t, hiBand, hjBand⟩
    let Z : ℤ := g * (p : ℤ) - t * q
    refine ⟨Z, ?_, ?_⟩
    · exact (mem_sparsePositiveLattice_iff
          (i' : ℤ) (j' : ℤ) X Y q Z hiZ hijZ hgapZ hqZ hbez).2
        ⟨hiBand, hjBand⟩
    · rw [Int.modEq_iff_dvd]
      refine ⟨t, ?_⟩
      dsimp [Z]
      ring
  · rintro ⟨Z, hZ, hmod⟩
    have hbands := (mem_sparsePositiveLattice_iff
      (i' : ℤ) (j' : ℤ) X Y q Z hiZ hijZ hgapZ hqZ hbez).1 hZ
    rw [Int.modEq_iff_dvd] at hmod
    obtain ⟨t, ht⟩ := hmod
    have hZeq : Z = g * (p : ℤ) - t * q := by
      nlinarith
    refine ⟨t, ?_, ?_⟩
    · simpa [hZeq] using hbands.1
    · simpa [hZeq] using hbands.2

open Classical in
/-- A finite set of nonzero, pairwise incongruent integer representatives has
the same cardinality as its nonzero residue image. -/
theorem residueRepresentatives_card
    (S : Finset ℤ) (q : ℕ) (hq : 0 < q)
    (hnonzero : ∀ Z ∈ S, ¬ (q : ℤ) ∣ Z)
    (hinjective : ∀ Z₁ ∈ S, ∀ Z₂ ∈ S,
      Z₁ ≡ Z₂ [ZMOD (q : ℤ)] → Z₁ = Z₂) :
    ((Finset.Ioo 0 q).filter fun r : ℕ =>
      ∃ Z ∈ S, Z ≡ (r : ℤ) [ZMOD (q : ℤ)]).card = S.card := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  symm
  apply Finset.card_bij
      (fun Z (_ : Z ∈ S) => (Z % (q : ℤ)).toNat)
  · intro Z hZ
    rw [Finset.mem_filter, Finset.mem_Ioo]
    have hmodNonneg : (0 : ℤ) ≤ Z % (q : ℤ) :=
      Int.emod_nonneg _ (by omega)
    have hmodLt : Z % (q : ℤ) < (q : ℤ) :=
      Int.emod_lt_of_pos _ hqZ
    have hcast : (((Z % (q : ℤ)).toNat : ℕ) : ℤ) = Z % (q : ℤ) :=
      Int.toNat_of_nonneg hmodNonneg
    have htoNatLt : (Z % (q : ℤ)).toNat < q := by omega
    refine ⟨⟨?_, htoNatLt⟩, Z, hZ, ?_⟩
    · by_contra hnot
      have hzero : Z % (q : ℤ) = 0 := by omega
      exact hnonzero Z hZ (Int.dvd_of_emod_eq_zero hzero)
    · show Z % (q : ℤ) =
          ((((Z % (q : ℤ)).toNat : ℕ) : ℤ) % (q : ℤ))
      rw [hcast, Int.emod_emod_of_dvd _ dvd_rfl]
  · intro Z₁ hZ₁ Z₂ hZ₂ heq
    apply hinjective Z₁ hZ₁ Z₂ hZ₂
    show Z₁ % (q : ℤ) = Z₂ % (q : ℤ)
    have hmodNonneg₁ : (0 : ℤ) ≤ Z₁ % (q : ℤ) :=
      Int.emod_nonneg _ (by omega)
    have hmodNonneg₂ : (0 : ℤ) ≤ Z₂ % (q : ℤ) :=
      Int.emod_nonneg _ (by omega)
    have hcast₁ : (((Z₁ % (q : ℤ)).toNat : ℕ) : ℤ) = Z₁ % (q : ℤ) :=
      Int.toNat_of_nonneg hmodNonneg₁
    have hcast₂ : (((Z₂ % (q : ℤ)).toNat : ℕ) : ℤ) = Z₂ % (q : ℤ) :=
      Int.toNat_of_nonneg hmodNonneg₂
    omega
  · intro r hr
    rw [Finset.mem_filter, Finset.mem_Ioo] at hr
    obtain ⟨⟨hr0, hrq⟩, Z, hZ, hmod⟩ := hr
    refine ⟨Z, hZ, ?_⟩
    have hrSmall : (r : ℤ) % (q : ℤ) = (r : ℤ) :=
      Int.emod_eq_of_lt (by exact_mod_cast hr0.le) (by exact_mod_cast hrq)
    have hmodNonneg : (0 : ℤ) ≤ Z % (q : ℤ) :=
      Int.emod_nonneg _ (by omega)
    have hcast : (((Z % (q : ℤ)).toNat : ℕ) : ℤ) = Z % (q : ℤ) :=
      Int.toNat_of_nonneg hmodNonneg
    change Z % (q : ℤ) = (r : ℤ) % (q : ℤ) at hmod
    have hmod' : Z % (q : ℤ) = (r : ℤ) := by
      simpa [hrSmall] using hmod
    omega

/-- Every two points of a canonical sparse lattice lie less than one modulus
apart.  This is the exact fact that prevents two interval points from mapping
to the same residue. -/
theorem sparsePositiveLattice_abs_sub_lt
    (i j X Y q : ℤ)
    (hi : 0 < i) (hij : i ≤ j) (hgap : j - i ≤ 13)
    (hsum : i + j ≤ 27) (hq : 0 < q)
    (hbez : j * X - i * Y = 1)
    {Z₁ Z₂ : ℤ}
    (hZ₁ : Z₁ ∈ sparsePositiveLattice i j X Y q)
    (hZ₂ : Z₂ ∈ sparsePositiveLattice i j X Y q) :
    |Z₁ - Z₂| < q := by
  have hj : 0 < j := lt_of_lt_of_le hi hij
  have hCi : 0 < 14 * i := by positivity
  have hCj : 0 < 14 * j := by positivity
  have activeBounds : ∀ Z ∈ sparsePositiveLattice i j X Y q,
      (14 * X - 1) * q < Z * (14 * i) ∧
      Z * (14 * j) < (14 * Y + 1) * q := by
    intro Z hZ
    rw [sparsePositiveLattice, Finset.mem_Ioo] at hZ
    have hlower := (Int.ediv_lt_iff_lt_mul hCi).mp hZ.1
    have hupperDiv : Z ≤ ((14 * Y + 1) * q - 1) / (14 * j) := by
      omega
    have hupperWeak := (Int.le_ediv_iff_mul_le hCj).mp hupperDiv
    exact ⟨hlower, by omega⟩
  have hsumCoefficient :
      i * (14 * Y + 1) - j * (14 * X - 1) = i + j - 14 := by
    nlinarith
  have hproduct : 1 ≤ i * j := by
    nlinarith [mul_pos hi hj]
  have orderedDistance : ∀ A B : ℤ,
      A ∈ sparsePositiveLattice i j X Y q →
      B ∈ sparsePositiveLattice i j X Y q →
      A ≤ B → B - A < q := by
    intro A B hA hB hAB
    obtain ⟨hAlower, -⟩ := activeBounds A hA
    obtain ⟨-, hBupper⟩ := activeBounds B hB
    have hAlowerScaled :
        j * ((14 * X - 1) * q) < j * (A * (14 * i)) :=
      mul_lt_mul_of_pos_left hAlower hj
    have hBupperScaled :
        i * (B * (14 * j)) < i * ((14 * Y + 1) * q) :=
      mul_lt_mul_of_pos_left hBupper hi
    have hscaled : 14 * i * j * (B - A) < (i + j - 14) * q := by
      nlinarith
    have hcoefficient : i + j - 14 ≤ 13 := by omega
    have hright : (i + j - 14) * q ≤ 13 * q :=
      mul_le_mul_of_nonneg_right hcoefficient (le_of_lt hq)
    have hnonnegative : 0 ≤ B - A := sub_nonneg.mpr hAB
    have hleft : 14 * (B - A) ≤ 14 * i * j * (B - A) := by
      nlinarith [mul_nonneg (sub_nonneg.mpr hproduct) hnonnegative]
    nlinarith
  rcases le_total Z₁ Z₂ with h12 | h21
  · rw [abs_of_nonpos (sub_nonpos.mpr h12)]
    simpa only [neg_sub] using orderedDistance Z₁ Z₂ hZ₁ hZ₂ h12
  · rw [abs_of_nonneg (sub_nonneg.mpr h21)]
    exact orderedDistance Z₂ Z₁ hZ₂ hZ₁ h21

open Classical in
/-- The nonzero residue classes represented by the sparse lattice have exactly
the floor-difference cardinality. -/
theorem sparsePositiveResidues_card
    (i' j' : ℕ) (X Y : ℤ) (q : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hij : i' ≤ j')
    (hgap : j' - i' ≤ 13) (hsum : i' + j' ≤ 27)
    (hq : 0 < q) (hbez : (j' : ℤ) * X - (i' : ℤ) * Y = 1) :
    ((Finset.Ioo 0 q).filter fun r : ℕ =>
      ∃ Z ∈ sparsePositiveLattice (i' : ℤ) (j' : ℤ) X Y q,
        Z ≡ (r : ℤ) [ZMOD (q : ℤ)]).card =
      ((((14 * Y + 1) * q - 1) / (14 * (j' : ℤ))) -
        (((14 * X - 1) * q) / (14 * (i' : ℤ)))).toNat := by
  have hiZ : (0 : ℤ) < (i' : ℤ) := by exact_mod_cast hi
  have hijZ : (i' : ℤ) ≤ (j' : ℤ) := by exact_mod_cast hij
  have hgapZ : (j' : ℤ) - (i' : ℤ) ≤ 13 := by
    exact_mod_cast hgap
  have hsumZ : (i' : ℤ) + (j' : ℤ) ≤ 27 := by
    exact_mod_cast hsum
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  rw [residueRepresentatives_card
    (sparsePositiveLattice (i' : ℤ) (j' : ℤ) X Y q) q hq]
  · exact sparsePositiveLattice_card (i' : ℤ) (j' : ℤ) X Y q
  · intro Z hZ hdvd
    obtain ⟨m, hm⟩ := hdvd
    have hbands := (mem_sparsePositiveLattice_iff
      (i' : ℤ) (j' : ℤ) X Y q Z hiZ hijZ hgapZ hqZ hbez).1 hZ
    rw [hm] at hbands
    apply branch_no_qmultiple i' j' X Y m q hq hi hj hbez
    · simpa [mul_comm, mul_left_comm, mul_assoc] using hbands.1
    · simpa [mul_comm, mul_left_comm, mul_assoc] using hbands.2
  · intro Z₁ hZ₁ Z₂ hZ₂ hmod
    have habs := sparsePositiveLattice_abs_sub_lt
      (i' : ℤ) (j' : ℤ) X Y q hiZ hijZ hgapZ hsumZ hqZ hbez
      hZ₁ hZ₂
    have habs' : |Z₂ - Z₁| < (q : ℤ) := by
      simpa [abs_sub_comm] using habs
    have hdvd : (q : ℤ) ∣ Z₂ - Z₁ :=
      Int.modEq_iff_dvd.mp hmod
    have hzero := Int.eq_zero_of_abs_lt_dvd hdvd habs'
    omega

open Classical in
/-- Exact positive-branch multiplier count.  This is the finite-`q` floor
formula left open after THM-964; the mirror theorem in `LRCSparseBranch` gives
the same count for the negative branch. -/
theorem branch_one_count
    (g : ℤ) (i' j' : ℕ) (X Y : ℤ) (q : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hij : i' ≤ j')
    (hgap : j' - i' ≤ 13) (hsum : i' + j' ≤ 27)
    (hq : 0 < q) (hcop : Nat.Coprime i' j')
    (hgcd : Int.gcd g (q : ℤ) = 1)
    (hbez : (j' : ℤ) * X - (i' : ℤ) * Y = 1) :
    ((Finset.Ioo 0 q).filter fun p : ℕ =>
      ∃ wa wb : ℤ,
        14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
        14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
        (j' : ℤ) * wa - (i' : ℤ) * wb = 1).card =
      ((((14 * Y + 1) * q - 1) / (14 * (j' : ℤ))) -
        (((14 * X - 1) * q) / (14 * (i' : ℤ)))).toNat := by
  let S : Finset ℤ :=
    sparsePositiveLattice (i' : ℤ) (j' : ℤ) X Y q
  let P : ℕ → Prop := fun r =>
    ∃ Z ∈ S, Z ≡ (r : ℤ) [ZMOD (q : ℤ)]
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hfilter :
      ((Finset.Ioo 0 q).filter fun p : ℕ =>
        ∃ wa wb : ℤ,
          14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
          14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
          (j' : ℤ) * wa - (i' : ℤ) * wb = 1) =
        ((Finset.Ioo 0 q).filter fun p : ℕ =>
          P ((g * (p : ℤ) % (q : ℤ)).toNat)) := by
    apply Finset.filter_congr
    intro p _
    rw [branch_one_iff_lattice_modEq g i' j' X Y q p
      hi hj hij hgap hq hcop hbez]
    change (∃ Z ∈ S, Z ≡ g * (p : ℤ) [ZMOD (q : ℤ)]) ↔
      ∃ Z ∈ S,
        Z ≡ (((g * (p : ℤ) % (q : ℤ)).toNat : ℕ) : ℤ)
          [ZMOD (q : ℤ)]
    have hmodNonneg :
        (0 : ℤ) ≤ g * (p : ℤ) % (q : ℤ) :=
      Int.emod_nonneg _ (by omega)
    have hcast :
        ((((g * (p : ℤ) % (q : ℤ)).toNat : ℕ) : ℤ)) =
          g * (p : ℤ) % (q : ℤ) :=
      Int.toNat_of_nonneg hmodNonneg
    have hresidue :
        g * (p : ℤ) ≡
          (((g * (p : ℤ) % (q : ℤ)).toNat : ℕ) : ℤ)
            [ZMOD (q : ℤ)] := by
      show (g * (p : ℤ)) % (q : ℤ) =
        ((((g * (p : ℤ) % (q : ℤ)).toNat : ℕ) : ℤ) % (q : ℤ))
      rw [hcast, Int.emod_emod_of_dvd _ dvd_rfl]
    constructor
    · rintro ⟨Z, hZ, hZmod⟩
      exact ⟨Z, hZ, hZmod.trans hresidue⟩
    · rintro ⟨Z, hZ, hZmod⟩
      exact ⟨Z, hZ, hZmod.trans hresidue.symm⟩
  rw [hfilter]
  rw [card_mod_filter_eq g q hq hgcd P]
  exact sparsePositiveResidues_card i' j' X Y q
    hi hj hij hgap hsum hq hbez

/-- Joint-failure multipliers for one reduced rational pair. -/
noncomputable def rationalPairFailMultipliers
    (g : ℤ) (i' j' q : ℕ) : Finset ℕ := by
  classical
  exact (Finset.Ioo 0 q).filter fun p : ℕ =>
    ∃ wa wb : ℤ,
      14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
      14 * |(g * j') * (p : ℤ) - wb * q| < q

/-- The subfamily whose unique witnesses have Bezout residue `k`. -/
noncomputable def rationalPairBranchMultipliers
    (g : ℤ) (i' j' q : ℕ) (k : ℤ) : Finset ℕ := by
  classical
  exact (Finset.Ioo 0 q).filter fun p : ℕ =>
    ∃ wa wb : ℤ,
      14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
      14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
      (j' : ℤ) * wa - (i' : ℤ) * wb = k

/-- Different Bezout residues give disjoint multiplier families because both
runner witnesses are unique. -/
theorem rationalPairBranchMultipliers_disjoint
    (g : ℤ) (i' j' q : ℕ) (hq : 0 < q) {k l : ℤ} (hkl : k ≠ l) :
    Disjoint (rationalPairBranchMultipliers g i' j' q k)
      (rationalPairBranchMultipliers g i' j' q l) := by
  classical
  rw [Finset.disjoint_left]
  intro p hpk hpl
  rw [rationalPairBranchMultipliers, Finset.mem_filter] at hpk hpl
  obtain ⟨-, wa₁, wb₁, ha₁, hb₁, hk⟩ := hpk
  obtain ⟨-, wa₂, wb₂, ha₂, hb₂, hl⟩ := hpl
  have hwa : wa₁ = wa₂ :=
    witness_unique (g * i') wa₁ wa₂ q p hq ha₁ ha₂
  have hwb : wb₁ = wb₂ :=
    witness_unique (g * j') wb₁ wb₂ q p hq hb₁ hb₂
  exact hkl <| by
    calc
      k = (j' : ℤ) * wa₁ - (i' : ℤ) * wb₁ := hk.symm
      _ = (j' : ℤ) * wa₂ - (i' : ℤ) * wb₂ := by rw [hwa, hwb]
      _ = l := hl

/-- Under the three-branch bound, the joint-failure family is the disjoint
union of the zero, positive, and negative Bezout branches. -/
theorem rationalPairFailMultipliers_eq_threeBranches
    (g : ℤ) (i' j' q : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hsum : i' + j' ≤ 27)
    (hq : 0 < q) :
    rationalPairFailMultipliers g i' j' q =
      rationalPairBranchMultipliers g i' j' q 0 ∪
        (rationalPairBranchMultipliers g i' j' q 1 ∪
          rationalPairBranchMultipliers g i' j' q (-1)) := by
  classical
  ext p
  simp only [rationalPairFailMultipliers, rationalPairBranchMultipliers,
    Finset.mem_filter, Finset.mem_union]
  constructor
  · rintro ⟨hp, wa, wb, ha, hb⟩
    have hk := rational_branch_bound
      (g * i') (g * j') wa wb i' j' q p hi hj hsum hq (by ring) ha hb
    have hk' : -(1 : ℤ) ≤ (j' : ℤ) * wa - (i' : ℤ) * wb ∧
        (j' : ℤ) * wa - (i' : ℤ) * wb ≤ 1 :=
      (abs_le.mp hk)
    have hcases :
        (j' : ℤ) * wa - (i' : ℤ) * wb = 0 ∨
        (j' : ℤ) * wa - (i' : ℤ) * wb = 1 ∨
        (j' : ℤ) * wa - (i' : ℤ) * wb = -1 := by
      omega
    rcases hcases with hzero | hpos | hneg
    · exact Or.inl ⟨hp, wa, wb, ha, hb, hzero⟩
    · exact Or.inr (Or.inl ⟨hp, wa, wb, ha, hb, hpos⟩)
    · exact Or.inr (Or.inr ⟨hp, wa, wb, ha, hb, hneg⟩)
  · rintro (hzero | hpos | hneg)
    · obtain ⟨hp, wa, wb, ha, hb, -⟩ := hzero
      exact ⟨hp, wa, wb, ha, hb⟩
    · obtain ⟨hp, wa, wb, ha, hb, -⟩ := hpos
      exact ⟨hp, wa, wb, ha, hb⟩
    · obtain ⟨hp, wa, wb, ha, hb, -⟩ := hneg
      exact ⟨hp, wa, wb, ha, hb⟩

open Classical in
/-- Complete finite-`q` joint-failure count for every canonically ordered
reduced pair in the three-branch range.  This upgrades the recon-exact full
pair ledger after THM-964 to a kernel-checked theorem. -/
theorem rationalPairFailMultipliers_card
    (g : ℤ) (i' j' : ℕ) (X Y : ℤ) (q : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hij : i' ≤ j')
    (hgap : j' - i' ≤ 13) (hsum : i' + j' ≤ 27)
    (hq : 0 < q) (hcop : Nat.Coprime i' j')
    (hgcd : Int.gcd g (q : ℤ) = 1)
    (hbez : (j' : ℤ) * X - (i' : ℤ) * Y = 1) :
    (rationalPairFailMultipliers g i' j' q).card =
      2 * ((q - 1) / (14 * max i' j')) +
      2 * ((((14 * Y + 1) * q - 1) / (14 * (j' : ℤ))) -
        (((14 * X - 1) * q) / (14 * (i' : ℤ)))).toNat := by
  let B₀ := rationalPairBranchMultipliers g i' j' q 0
  let B₁ := rationalPairBranchMultipliers g i' j' q 1
  let Bneg := rationalPairBranchMultipliers g i' j' q (-1)
  have hdisjoint₀₁ : Disjoint B₀ B₁ :=
    rationalPairBranchMultipliers_disjoint g i' j' q hq (by norm_num)
  have hdisjoint₀neg : Disjoint B₀ Bneg :=
    rationalPairBranchMultipliers_disjoint g i' j' q hq (by norm_num)
  have hdisjoint₁neg : Disjoint B₁ Bneg :=
    rationalPairBranchMultipliers_disjoint g i' j' q hq (by norm_num)
  have hzero : B₀.card = 2 * ((q - 1) / (14 * max i' j')) := by
    simpa [B₀, rationalPairBranchMultipliers] using
      branch_zero_count g i' j' q hq hi hj hcop hgcd
  have hpositive : B₁.card =
      ((((14 * Y + 1) * q - 1) / (14 * (j' : ℤ))) -
        (((14 * X - 1) * q) / (14 * (i' : ℤ)))).toNat := by
    simpa [B₁, rationalPairBranchMultipliers] using
      branch_one_count g i' j' X Y q hi hj hij hgap hsum hq hcop hgcd hbez
  have hmirror : B₁.card = Bneg.card := by
    simpa [B₁, Bneg, rationalPairBranchMultipliers] using
      branch_mirror_card g i' j' q
  rw [rationalPairFailMultipliers_eq_threeBranches g i' j' q hi hj hsum hq,
    Finset.card_union_of_disjoint
      (Finset.disjoint_union_right.mpr ⟨hdisjoint₀₁, hdisjoint₀neg⟩),
    Finset.card_union_of_disjoint hdisjoint₁neg,
    hzero, hpositive, ← hmirror, hpositive]
  omega

/-! ## Axiom audit -/

#print axioms mem_sparsePositiveLattice_iff
#print axioms sparsePositiveLattice_card
#print axioms branch_one_iff_lattice_modEq
#print axioms residueRepresentatives_card
#print axioms sparsePositiveLattice_abs_sub_lt
#print axioms sparsePositiveResidues_card
#print axioms branch_one_count
#print axioms rationalPairBranchMultipliers_disjoint
#print axioms rationalPairFailMultipliers_eq_threeBranches
#print axioms rationalPairFailMultipliers_card

end

end LRC14Concrete
end LonelyRunner
