/-
  TournamentH7.StaircaseTileModel — Concrete tile-coordinate model
  (project-novel, oracle-2026-05-29-S4)

  ─── What this module provides ────────────────────────────────────────
  The project canon distinguishes two different involutions on tilings:
   • *Tile-complement* (`tilde T` in `Tilings.lean`):  flip every tile bit.
   • *Grid reflection* (THM-280):  apply the coordinate map
     (x, y) ↦ (n+1−y, n+1−x), permuting which tile gets which bit.

  This module introduces explicit `Tile n` coordinates and the
  `Tiling n := Tile n → Bool` representation, enabling formal statement
  of THM-280.

  ─── Definitions ──────────────────────────────────────────────────────
   • `StTile n` — a tile, identified with a non-consecutive pair
     `(hi, lo)` of vertices with `lo + 2 ≤ hi`.
   • `StTiling n := StTile n → Bool` — a tiling.
   • `StTiling.toTournament` — induced tournament.
   • `StTiling.reflect` — the grid-reflection involution on tilings.
   • `StTiling.complement` — the tile-complement involution (flip all
     bits), connecting to the `Tournament.tilde` notion.

  ─── THM-280 statement (project-novel, opus-2026-04-03-S27) ──────────
  `T(reflect b) ≅ σ (T(b)^op)` where σ = vertex-reversal `v ↦ n − 1 − v`.

  ─── Status ──────────────────────────────────────────────────────────
  This module gives DEFINITIONS only.  The THM-280 proof is left as a
  future-work axiom; see canonical proof in `THM-280-grid-reflection-...md`.
-/

import TournamentH7.Basic
import TournamentH7.Tilings
import TournamentH7.GridReflection
import Mathlib.Tactic

namespace Tournament

variable {n : ℕ}

/-! ### Staircase tile coordinates -/

/-- A *staircase tile* on n vertices: a non-consecutive vertex pair
    `(hi, lo) : Fin n × Fin n` with `lo.val + 2 ≤ hi.val`. -/
structure StTile (n : ℕ) where
  hi  : Fin n
  lo  : Fin n
  gap : lo.val + 2 ≤ hi.val

namespace StTile
variable {n : ℕ}

@[simp] lemma lt (t : StTile n) : t.lo.val < t.hi.val := by
  have := t.gap; omega

@[simp] lemma ne (t : StTile n) : t.lo ≠ t.hi := by
  intro h; have := t.gap
  have : t.lo.val = t.hi.val := congrArg Fin.val h; omega

end StTile

/-! ### Staircase tilings -/

/-- A staircase tiling: an assignment of a Bool to each tile.
    Convention (project canon):
       false  ⟹  arc hi → lo (downward, "default")
       true   ⟹  arc lo → hi (upward) -/
def StTiling (n : ℕ) : Type := StTile n → Bool

/-! ### The grid reflection -/

/-- The grid reflection on tile coordinates: `(hi, lo) ↦ (n−1−lo, n−1−hi)`.
    (This is the 0-indexed version of the canonical map
     `(x, y) ↦ (n+1−y, n+1−x)`.) -/
def StTile.reflect (t : StTile n) : StTile n where
  hi := ⟨n - 1 - t.lo.val, by
    have h := t.lo.is_lt
    by_cases hn : n = 0
    · exfalso; exact absurd h (by simp [hn])
    · omega⟩
  lo := ⟨n - 1 - t.hi.val, by
    have h := t.hi.is_lt
    by_cases hn : n = 0
    · exfalso; exact absurd h (by simp [hn])
    · omega⟩
  gap := by
    have h_lt := t.lt
    have h_hi : t.hi.val < n := t.hi.is_lt
    have h_lo : t.lo.val < n := t.lo.is_lt
    have h_gap := t.gap
    show (n - 1 - t.hi.val) + 2 ≤ n - 1 - t.lo.val
    omega

/-- StTile equality from hi and lo equality (the gap proof is propositionally
    equal automatically). -/
@[ext] lemma StTile.ext {s t : StTile n} (h_hi : s.hi = t.hi)
    (h_lo : s.lo = t.lo) : s = t := by
  cases s; cases t; congr

/-! ### Finite/equality infrastructure for concrete tiles -/

/-- Forget a staircase tile to its legal coordinate pair. -/
def StTile.toGapPair (t : StTile n) :
    {p : Fin n × Fin n // p.2.val + 2 ≤ p.1.val} :=
  ⟨(t.hi, t.lo), t.gap⟩

/-- Build a staircase tile from a legal coordinate pair. -/
def StTile.ofGapPair
    (p : {p : Fin n × Fin n // p.2.val + 2 ≤ p.1.val}) : StTile n where
  hi := p.1.1
  lo := p.1.2
  gap := p.2

@[simp] lemma StTile.toGapPair_ofGapPair
    (p : {p : Fin n × Fin n // p.2.val + 2 ≤ p.1.val}) :
    StTile.toGapPair (StTile.ofGapPair p) = p := by
  cases p
  rfl

@[simp] lemma StTile.ofGapPair_toGapPair (t : StTile n) :
    StTile.ofGapPair (StTile.toGapPair t) = t := by
  cases t
  rfl

/-- Staircase tiles are equivalent to the finite subtype of legal coordinate pairs. -/
def StTile.equivGapPair (n : ℕ) :
    StTile n ≃ {p : Fin n × Fin n // p.2.val + 2 ≤ p.1.val} where
  toFun := StTile.toGapPair
  invFun := StTile.ofGapPair
  left_inv := StTile.ofGapPair_toGapPair
  right_inv := StTile.toGapPair_ofGapPair

instance StTile.instDecidableEq (n : ℕ) : DecidableEq (StTile n) := by
  classical
  exact fun s t =>
    decidable_of_iff (StTile.toGapPair s = StTile.toGapPair t) (by
      constructor
      · intro h
        have := congrArg StTile.ofGapPair h
        simpa using this
      · intro h
        exact congrArg StTile.toGapPair h)

noncomputable instance StTile.instFintype (n : ℕ) : Fintype (StTile n) := by
  classical
  exact Fintype.ofEquiv _ (StTile.equivGapPair n).symm

namespace StTiling

noncomputable instance instFintype (n : ℕ) : Fintype (StTiling n) := by
  dsimp [StTiling]
  infer_instance

noncomputable instance instDecidableEq (n : ℕ) : DecidableEq (StTiling n) := by
  classical
  dsimp [StTiling]
  infer_instance

end StTiling

/-- The grid reflection is an involution. -/
lemma StTile.reflect_reflect (t : StTile n) : t.reflect.reflect = t := by
  have h_hi : t.hi.val < n := t.hi.is_lt
  have h_lo : t.lo.val < n := t.lo.is_lt
  apply StTile.ext
  · apply Fin.ext
    show n - 1 - (n - 1 - t.hi.val) = t.hi.val
    omega
  · apply Fin.ext
    show n - 1 - (n - 1 - t.lo.val) = t.lo.val
    omega

/-- The reflected tiling: `b.reflect t = b (t.reflect)`. -/
def StTiling.reflect (b : StTiling n) : StTiling n :=
  fun t => b t.reflect

/-- The grid reflection is an involution on tilings. -/
lemma StTiling.reflect_reflect (b : StTiling n) : b.reflect.reflect = b := by
  funext t
  show b t.reflect.reflect = b t
  rw [StTile.reflect_reflect]

/-! ### The tile-complement (distinct from reflection!) -/

/-- The tile-complement of a tiling: flip every bit. -/
def StTiling.complement (b : StTiling n) : StTiling n := fun t => !(b t)

lemma StTiling.complement_complement (b : StTiling n) :
    b.complement.complement = b := by
  funext t; simp [complement]

/-- Reflection and complement commute. -/
lemma StTiling.reflect_complement (b : StTiling n) :
    b.reflect.complement = b.complement.reflect := by
  funext t; rfl

/-! ### Tournament from a tiling

  The arc relation: for i, j ∈ Fin n, i ≠ j:
   • if consecutive (|i.val − j.val| = 1): base path arc (higher to lower)
   • if non-consecutive: determined by b at the tile (hi, lo) = (max, min).
-/

/-- The arc relation induced by a staircase tiling. -/
def StTiling.arc (b : StTiling n) (i j : Fin n) : Bool :=
  if i = j then false
  else if _h₁ : i.val = j.val + 1 then true
  else if _h₂ : j.val = i.val + 1 then false
  else if h₃ : j.val + 2 ≤ i.val then
    !(b ⟨i, j, h₃⟩)
  else if h₄ : i.val + 2 ≤ j.val then
    b ⟨j, i, h₄⟩
  else false  -- unreachable when i ≠ j

@[simp] lemma StTiling.arc_self (b : StTiling n) (i : Fin n) :
    b.arc i i = false := by
  simp [StTiling.arc]

lemma StTiling.arc_succ_down (b : StTiling n) (i : Fin n)
    (h : i.val + 1 < n) :
    b.arc ⟨i.val + 1, h⟩ i = true := by
  have hne : (⟨i.val + 1, h⟩ : Fin n) ≠ i := by
    intro heq
    have := congrArg Fin.val heq
    simp at this
  simp [StTiling.arc, hne]

lemma StTiling.arc_succ_up (b : StTiling n) (i : Fin n)
    (h : i.val + 1 < n) :
    b.arc i ⟨i.val + 1, h⟩ = false := by
  have hne : i ≠ (⟨i.val + 1, h⟩ : Fin n) := by
    intro heq
    have := congrArg Fin.val heq
    simp at this
  have hnot : ¬ i.val = (⟨i.val + 1, h⟩ : Fin n).val + 1 := by
    simp
    omega
  simp [StTiling.arc, hne, hnot]

lemma StTiling.arc_nonconsec_down (b : StTiling n) {hi lo : Fin n}
    (hgap : lo.val + 2 ≤ hi.val) :
    b.arc hi lo = !(b ⟨hi, lo, hgap⟩) := by
  have hne : hi ≠ lo := by
    intro heq
    have := congrArg Fin.val heq
    omega
  have hnot₁ : ¬ hi.val = lo.val + 1 := by omega
  have hnot₂ : ¬ lo.val = hi.val + 1 := by omega
  simp [StTiling.arc, hne, hnot₁, hnot₂, hgap]

lemma StTiling.arc_nonconsec_up (b : StTiling n) {lo hi : Fin n}
    (hgap : lo.val + 2 ≤ hi.val) :
    b.arc lo hi = b ⟨hi, lo, hgap⟩ := by
  have hne : lo ≠ hi := by
    intro heq
    have := congrArg Fin.val heq
    omega
  have hnot₁ : ¬ lo.val = hi.val + 1 := by omega
  have hnot₂ : ¬ hi.val = lo.val + 1 := by omega
  have hnot₃ : ¬ hi.val + 2 ≤ lo.val := by omega
  simp [StTiling.arc, hne, hnot₁, hnot₂, hnot₃, hgap]

lemma StTiling.arc_tile_down (b : StTiling n) (t : StTile n) :
    b.arc t.hi t.lo = !(b t) := by
  cases t
  exact StTiling.arc_nonconsec_down b ‹_›

lemma StTiling.arc_tile_up (b : StTiling n) (t : StTile n) :
    b.arc t.lo t.hi = b t := by
  cases t
  exact StTiling.arc_nonconsec_up b ‹_›

private lemma StTiling.arc_exactly_one (b : StTiling n) {i j : Fin n}
    (hij : i ≠ j) :
    (b.arc i j = true ∧ b.arc j i = false) ∨
    (b.arc i j = false ∧ b.arc j i = true) := by
  have hval_ne : i.val ≠ j.val := fun h => hij (Fin.ext h)
  rcases lt_or_gt_of_ne hval_ne with hlt | hgt
  · by_cases hsucc : i.val + 1 = j.val
    · have hsucc_lt : i.val + 1 < n := by
        rw [hsucc]
        exact j.is_lt
      have hj_eq : j = ⟨i.val + 1, hsucc_lt⟩ := Fin.ext hsucc.symm
      subst hj_eq
      right
      exact ⟨StTiling.arc_succ_up b i hsucc_lt,
        StTiling.arc_succ_down b i hsucc_lt⟩
    · have hgap : i.val + 2 ≤ j.val := by omega
      let t : StTile n := ⟨j, i, hgap⟩
      have hij_arc : b.arc i j = b t := StTiling.arc_nonconsec_up b hgap
      have hji_arc : b.arc j i = !(b t) := StTiling.arc_nonconsec_down b hgap
      cases ht : b t
      · right
        rw [hij_arc, hji_arc, ht]
        simp
      · left
        rw [hij_arc, hji_arc, ht]
        simp
  · by_cases hsucc : j.val + 1 = i.val
    · have hsucc_lt : j.val + 1 < n := by
        rw [hsucc]
        exact i.is_lt
      have hi_eq : i = ⟨j.val + 1, hsucc_lt⟩ := Fin.ext hsucc.symm
      subst hi_eq
      left
      exact ⟨StTiling.arc_succ_down b j hsucc_lt,
        StTiling.arc_succ_up b j hsucc_lt⟩
    · have hgap : j.val + 2 ≤ i.val := by omega
      let t : StTile n := ⟨i, j, hgap⟩
      have hij_arc : b.arc i j = !(b t) := StTiling.arc_nonconsec_down b hgap
      have hji_arc : b.arc j i = b t := StTiling.arc_nonconsec_up b hgap
      cases ht : b t
      · left
        rw [hij_arc, hji_arc, ht]
        simp
      · right
        rw [hij_arc, hji_arc, ht]
        simp

/-- The tournament induced by a concrete staircase tiling. -/
def StTiling.toTournament (b : StTiling n) : Tournament n where
  arc := b.arc
  irrefl := StTiling.arc_self b
  total := by
    intro i j hij
    rcases StTiling.arc_exactly_one b hij with h | h
    · exact Or.inl h.1
    · exact Or.inr h.2
  asym := by
    intro i j hboth
    by_cases hij : i = j
    · subst hij
      simp at hboth
    · rcases StTiling.arc_exactly_one b hij with h | h
      · rw [h.2] at hboth
        exact Bool.false_eq_true.mp hboth.2
      · rw [h.1] at hboth
        exact Bool.false_eq_true.mp hboth.1

@[simp] lemma StTiling.toTournament_arc (b : StTiling n) (i j : Fin n) :
    b.toTournament.arc i j = b.arc i j := rfl

theorem StTiling.toTournament_hasBasePath (b : StTiling n) :
    HasBasePath b.toTournament := by
  intro i h
  exact StTiling.arc_succ_down b i h

/-! ### THM-280 statement (axiomatised; proof requires further work) -/

/-- **Axiom (THM-280, project-novel, opus-2026-04-03-S27).**

    The tournament from the *grid-reflected tiling* `b.reflect` is
    isomorphic to the vertex-reversed op-tournament of the original
    `T(b)`.  At the arc level:

        (StTiling.reflect b).arc i j =
            (b.arc) (vertexReversal n j) (vertexReversal n i).

    Proof (project canon): direct case-by-case analysis on consecutive
    vs non-consecutive arcs and the action of grid reflection on tile
    coordinates.  Lean formalisation deferred.

    Note: this is DIFFERENT from `tilde T`, which corresponds to
    `StTiling.complement b`. -/
axiom thm280_arc (b : StTiling n) (i j : Fin n) :
    (b.reflect).arc i j = b.arc (vertexReversal n j) (vertexReversal n i)

end Tournament
