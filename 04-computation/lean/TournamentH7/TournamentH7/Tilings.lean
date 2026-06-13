/-
  TournamentH7.Tilings — Base-path tournaments and the tile-complement
                          score formula (oracle-2026-05-11, novel).

  See file header in original commit; restructured for clarity.

  Defines:
   • `IsConsec i j` — symmetric "consecutive in base path" predicate.
   • `HasBasePath T` — predicate that T contains every consecutive arc.
   • `Tournament.tilde T` — tile-complement: same on consecutive arcs,
     negated on non-consecutive arcs.
   • `tilde_arc_basepath`, `tilde_arc_nonconsec` — simp lemmas.

  Project canon reference: oracle-2026-05-11-S1 memo.
-/

import TournamentH7.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Tactic.Linarith

open Finset

namespace Tournament

variable {n : ℕ}

/-! ### Consecutive pairs -/

/-- `i, j : Fin n` are *consecutive* iff their `Fin.val`s differ by 1. -/
def IsConsec (i j : Fin n) : Prop :=
  i.val + 1 = j.val ∨ j.val + 1 = i.val

instance (i j : Fin n) : Decidable (IsConsec i j) := by
  unfold IsConsec; infer_instance

@[symm] lemma IsConsec.symm {i j : Fin n} (h : IsConsec i j) : IsConsec j i := by
  rcases h with h | h
  · exact Or.inr h
  · exact Or.inl h

@[simp] lemma isConsec_comm {i j : Fin n} : IsConsec i j ↔ IsConsec j i :=
  ⟨IsConsec.symm, IsConsec.symm⟩

/-- `T` *has the base path* iff every consecutive arc `i+1 → i` is present. -/
def HasBasePath (T : Tournament n) : Prop :=
  ∀ (i : Fin n) (h : i.val + 1 < n),
    T.arc ⟨i.val + 1, h⟩ i = true

/-! ### Tile-complement T̃ -/

/-- The arc relation of T̃: same on consecutive pairs, negated on the rest. -/
def tildeArc (T : Tournament n) (i j : Fin n) : Bool :=
  if i = j then false
  else if IsConsec i j then T.arc i j
  else !(T.arc i j)

@[simp] lemma tildeArc_self (T : Tournament n) (i : Fin n) :
    tildeArc T i i = false := by simp [tildeArc]

lemma tildeArc_consec (T : Tournament n) {i j : Fin n}
    (hij : i ≠ j) (hcons : IsConsec i j) :
    tildeArc T i j = T.arc i j := by
  unfold tildeArc; simp [hij, hcons]

lemma tildeArc_nonconsec (T : Tournament n) {i j : Fin n}
    (hij : i ≠ j) (hnc : ¬ IsConsec i j) :
    tildeArc T i j = !(T.arc i j) := by
  unfold tildeArc; simp [hij, hnc]

/-- The crucial fact: tildeArc inherits T's "one of (i→j, j→i) is true" via
    symmetry of the IsConsec predicate. -/
private lemma tildeArc_exactly_one (T : Tournament n) {i j : Fin n}
    (hij : i ≠ j) :
    (tildeArc T i j = true ∧ tildeArc T j i = false) ∨
    (tildeArc T i j = false ∧ tildeArc T j i = true) := by
  have hji : j ≠ i := fun H => hij H.symm
  -- Use Boolean dichotomy on T.arc i j
  by_cases hcons : IsConsec i j
  · -- consecutive: tilde = T
    have hconsji : IsConsec j i := hcons.symm
    rw [tildeArc_consec T hij hcons, tildeArc_consec T hji hconsji]
    rcases T.total i j hij with hl | hr
    · left
      refine ⟨hl, ?_⟩
      cases h_eq : T.arc j i with
      | true =>  exact absurd (T.asym i j ⟨hl, h_eq⟩) (fun x => x)
      | false => rfl
    · right
      refine ⟨?_, hr⟩
      cases h_eq : T.arc i j with
      | true =>  exact absurd (T.asym i j ⟨h_eq, hr⟩) (fun x => x)
      | false => rfl
  · -- non-consecutive: tilde = !T
    have hncji : ¬ IsConsec j i := fun H => hcons H.symm
    rw [tildeArc_nonconsec T hij hcons, tildeArc_nonconsec T hji hncji]
    rcases T.total i j hij with hl | hr
    · right
      refine ⟨by rw [hl]; rfl, ?_⟩
      cases h_eq : T.arc j i with
      | true =>  exact absurd (T.asym i j ⟨hl, h_eq⟩) (fun x => x)
      | false => simp
    · left
      refine ⟨?_, by rw [hr]; rfl⟩
      cases h_eq : T.arc i j with
      | true =>  exact absurd (T.asym i j ⟨h_eq, hr⟩) (fun x => x)
      | false => simp

/-- Tile-complement T̃ as a `Tournament`. -/
def tilde (T : Tournament n) : Tournament n where
  arc := tildeArc T
  irrefl := tildeArc_self T
  total := by
    intro i j hij
    rcases tildeArc_exactly_one T hij with ⟨h1, _⟩ | ⟨_, h2⟩
    · exact Or.inl h1
    · exact Or.inr h2
  asym := by
    intro i j ⟨h1, h2⟩
    by_cases hij : i = j
    · simp [hij] at h1
    · rcases tildeArc_exactly_one T hij with ⟨_, h⟩ | ⟨h, _⟩
      · rw [h] at h2; exact absurd h2 (by decide)
      · rw [h] at h1; exact absurd h1 (by decide)

/-! ### Basic simp lemmas for `tilde` -/

@[simp] lemma tilde_arc (T : Tournament n) (i j : Fin n) :
    (tilde T).arc i j = tildeArc T i j := rfl

lemma tilde_arc_consec (T : Tournament n) {i j : Fin n}
    (hij : i ≠ j) (hcons : IsConsec i j) :
    (tilde T).arc i j = T.arc i j := by
  rw [tilde_arc]; exact tildeArc_consec T hij hcons

lemma tilde_arc_nonconsec (T : Tournament n) {i j : Fin n}
    (hij : i ≠ j) (hnc : ¬ IsConsec i j) :
    (tilde T).arc i j = !(T.arc i j) := by
  rw [tilde_arc]; exact tildeArc_nonconsec T hij hnc

/-- Tile-complement is an involution. -/
lemma tilde_tilde (T : Tournament n) (i j : Fin n) :
    (tilde (tilde T)).arc i j = T.arc i j := by
  by_cases hij : i = j
  · simp [hij, tilde, tildeArc, T.irrefl]
  · by_cases hcons : IsConsec i j
    · rw [tilde_arc_consec _ hij hcons, tilde_arc_consec _ hij hcons]
    · rw [tilde_arc_nonconsec _ hij hcons, tilde_arc_nonconsec _ hij hcons]
      cases T.arc i j <;> simp

/-! ### Score-formula bookkeeping -/

/-- The out-neighbours of `v` in T as a `Finset`. -/
def outNbrs (T : Tournament n) (v : Fin n) : Finset (Fin n) :=
  Finset.univ.filter (fun w => T.arc v w = true)

lemma mem_outNbrs (T : Tournament n) {v w : Fin n} :
    w ∈ outNbrs T v ↔ T.arc v w = true := by
  unfold outNbrs; simp

@[simp] lemma outDegree_eq_card_outNbrs (T : Tournament n) (v : Fin n) :
    T.outDegree v = (outNbrs T v).card := rfl

/-! ### The score formula for T̃  (novel; oracle-2026-05-11-S1)

  Proof sketch (counting argument):

  Each vertex v ∈ Fin n has (n − 1) neighbours, partitioned into:
    • consecutive neighbours C(v) ⊆ {v − 1, v + 1};
    • non-consecutive neighbours N(v).

  Cardinalities (for n ≥ 2):
    • v.val = 0       ⟹ |C(v)| = 1 (just v + 1),  |N(v)| = n − 2.
    • v.val = n − 1   ⟹ |C(v)| = 1 (just v − 1),  |N(v)| = n − 2.
    • 0 < v.val < n−1 ⟹ |C(v)| = 2,               |N(v)| = n − 3.

  In a tournament T with the base path (HasBasePath T):
    • From v we have base-path out-arcs:
        − v.val = 0:     0 (no v − 1 to point to),
        − v.val = n − 1: 1 (v → v − 1),
        − interior:      1 (v → v − 1).
    • The other consec arc (when present, v + 1 → v) is IN, not OUT.
    • Non-consec out-arcs from v: `s(v) − (1 if v.val ≠ 0 else 0)`.

  T̃ keeps consecutive arcs unchanged and flips non-consec arcs, so:
    • T̃'s base-path out-arcs from v = T's (above).
    • T̃'s non-consec out-arcs from v = |N(v)| − (T's non-consec out-arcs).

  Putting it together:
    s̃(0)     = 0 + (n − 2)       − s(0)       =  (n − 2) − s(0)
    s̃(n−1)  = 1 + (n − 2)       − (s(n−1) − 1) = n − s(n−1)
    s̃(v)     = 1 + (n − 3)       − (s(v) − 1)    = (n − 1) − s(v)   (interior)

  The full counting proof is left as an axiom for now (the proof is a
  routine but tedious bookkeeping; see canonical version in
  `00-navigation/SESSION-LOG.md`, oracle-2026-05-11-S1). -/

/-! ### Toward de-axiomatising the score formula

    Key identity: for any v, w ≠ v:
      • non-consecutive: indicator(T.arc v w) + indicator(tilde T.arc v w) = 1.
      • consecutive: indicator(T.arc v w) + indicator(tilde T.arc v w) = 2·indicator(T.arc v w).
    Summing over w ≠ v gives:
      outDeg T v + outDeg tilde T v = #nonconsec(v) + 2·#consec_out_T(v).

    For v.val = 0: #nonconsec = n - 2, #consec_out = 0 ⟹ sum = n - 2.
    For v.val = n-1: #nonconsec = n - 2, #consec_out = 1 ⟹ sum = n.
    For interior: #nonconsec = n - 3, #consec_out = 1 ⟹ sum = n - 1.

    The partition formulas are straightforward but lengthy.  We KEEP the
    axioms for now and add a `Verify`-side reproof sketch as a comment. -/

/-! ### Helper for the score formula -/

/-- For HasBasePath T and n ≥ 2, T.arc 0 1 = false. -/
private lemma T_arc_at_zero_to_one (T : Tournament n) (hbp : HasBasePath T)
    (hn : 2 ≤ n) :
    T.arc ⟨0, by omega⟩ ⟨1, by omega⟩ = false := by
  have h0 : 0 < n := by omega
  have h1 : 1 < n := by omega
  have h_succ : (⟨0, h0⟩ : Fin n).val + 1 < n := by show 0 + 1 < n; omega
  have h_bp : T.arc ⟨(⟨0, h0⟩ : Fin n).val + 1, h_succ⟩ ⟨0, h0⟩ = true :=
    hbp ⟨0, h0⟩ h_succ
  -- h_bp has type T.arc ⟨1, _⟩ ⟨0, _⟩ = true (after simp).
  have h_bp' : T.arc ⟨1, h1⟩ ⟨0, h0⟩ = true := by
    convert h_bp using 2
  cases h : T.arc ⟨0, h0⟩ ⟨1, h1⟩ with
  | false => rfl
  | true =>
    exact absurd (T.asym ⟨0, h0⟩ ⟨1, h1⟩ ⟨h, h_bp'⟩) (fun x => x)

/-- **Axiom (oracle-2026-05-11-S1, project-novel).** Score formula for T̃,
    base-path sink case.

    Proof outline (deferred): at vertex 0, the only consecutive neighbor
    is vertex 1; the base path gives 1 → 0 so T.arc 0 1 = false.  Hence
    T's out-neighbors at vertex 0 lie among {w : w.val ≥ 2}.  For non-
    consecutive w, exactly one of T.arc 0 w and tilde T.arc 0 w is true.
    Hence the two `outNbrs` sets are *complementary* in the size-(n-2) set
    {w : w.val ≥ 2}. -/
axiom tilde_score_sink (T : Tournament n) (hbp : HasBasePath T)
    (hn : 2 ≤ n) (v : Fin n) (hv : v.val = 0) :
    (tilde T).outDegree v + T.outDegree v = n - 2

/-- **Axiom (oracle-2026-05-11-S1, project-novel).** Score formula for T̃,
    base-path source case. -/
axiom tilde_score_source (T : Tournament n) (h : HasBasePath T)
    (hn : 2 ≤ n) (v : Fin n) (hv : v.val = n - 1) :
    (tilde T).outDegree v + T.outDegree v = n

/-- **Axiom (oracle-2026-05-11-S1, project-novel).** Score formula for T̃,
    interior vertex case. -/
axiom tilde_score_interior (T : Tournament n) (h : HasBasePath T)
    (hn : 2 ≤ n) (v : Fin n) (hv_lo : 0 < v.val) (hv_hi : v.val + 1 < n) :
    (tilde T).outDegree v + T.outDegree v = n - 1

/-! ### Consequence: regular tournaments are never SF -/

/-- A tournament is *regular* iff every score equals (n − 1)/2 (only
    possible for odd n). -/
def IsRegular (T : Tournament n) : Prop :=
  ∀ v : Fin n, 2 * T.outDegree v = n - 1

/-- **Theorem (project novel, oracle-2026-05-11-S1).**

    A regular tournament T cannot equal its tile-complement T̃ vertex-wise
    on out-degrees — i.e. it is not *self-flip* (SF) in the project's sense
    via the score-formula obstruction. Specifically, T̃'s score sequence at
    any interior vertex satisfies `s̃(v) + s(v) = n − 1`, so if T is regular
    then `s̃(v) = s(v)`, which is consistent — BUT at endpoints we get
    `s̃(0) + s(0) = n − 2 ≠ n − 1`, so `s̃(0) ≠ s(0)` and the score
    sequence of T̃ differs from T at vertex 0. Hence T ≇ T̃, so T is not SF. -/
theorem regular_not_SF (T : Tournament n) (hbp : HasBasePath T)
    (hn : 3 ≤ n) (hreg : IsRegular T)
    (v0 : Fin n) (hv0 : v0.val = 0) :
    (tilde T).outDegree v0 ≠ T.outDegree v0 := by
  have heq := tilde_score_sink T hbp (by omega) v0 hv0
  have hreg0 := hreg v0
  intro hsame
  -- heq:   (tilde T).outDegree v0 + T.outDegree v0 = n - 2
  -- hreg0: 2 * T.outDegree v0 = n - 1
  -- hsame: (tilde T).outDegree v0 = T.outDegree v0
  -- ⟹ 2 * T.outDegree v0 = n - 2 = n - 1: contradiction (since 2 | n-1 here).
  -- Note: in ℕ subtraction is truncated, but n ≥ 3 means we're safe.
  omega

end Tournament
