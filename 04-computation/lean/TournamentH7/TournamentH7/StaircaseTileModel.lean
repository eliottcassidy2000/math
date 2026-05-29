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
  else if h₁ : i.val = j.val + 1 then true
  else if h₂ : j.val = i.val + 1 then false
  else if h₃ : j.val + 2 ≤ i.val then
    !(b ⟨i, j, h₃⟩)
  else if h₄ : i.val + 2 ≤ j.val then
    b ⟨j, i, h₄⟩
  else false  -- unreachable when i ≠ j

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
