/-
  TournamentH7.LRCQWindow — THE Q-WINDOW LEMMA and THE WITNESS-LADDER RIGIDITY
  (death-star-2026-07-17-S41, HYP-7190; the two opens THM-947 named).

  * `failWitness_pos_of_window` — THE Q-WINDOW LEMMA: for `q ≤ 14·v i`, every band
    failure of runner `i` at `p ≥ 1` has witness `n ≥ 1` — the near-zero
    catastrophe (all runners bad at `p = 1` via `n = 0`) cannot occur inside the
    window; THM-947's pair constraints become near-proportionality of POSITIVE
    integers.
  * `failWitness_le` — the ceiling `n ≤ v i` for `p < q`: witnesses live in
    `[1, v]` inside the window.
  * `witness_ladder` — THE RIGIDITY ENGINE: simultaneously bad runners with
    `3·v_i ≤ v_j ≤ 13·v_i` (the residual's own ratio window: ladder steps ≥ 3,
    compression ≤ 13) and `n_i ≥ 1` force `n_j ≥ 3·n_i` — the witness vector
    INHERITS the ladder.  The `13` is sharp for this argument: the constraint's
    slack is `(v_i + v_j)/(14·v_i) < 1` exactly when `v_j < 13·v_i` (at 13,
    equality still lands `n_j > 3n_i − 1`).
  * `witness_chain_geometric` — iterated along a chain: `k+1` simultaneously bad
    runners with consecutive ratios in `[3, 13]` force `n_top ≥ 3^k·n_bottom`;
    seven force `n_top ≥ 729`.
  * `witness_pin` — THE PIN: witness `n ≥ N` forces `14·N·q − q < 14·v·p` — deep
    ladder overlaps live only FAR from zero, quantitatively (seven ladder bads:
    `p/q > 10205/(14·v_top)`).  On trapped strata this is the formal face of
    "deep coverage requires structure": the witness vector must be geometric, and
    the family's own magnitudes must accommodate it inside `p < q`.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCArcWire

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **THE Q-WINDOW LEMMA**: for `q ≤ 14·v i`, every band failure at `p ≥ 1` has
witness `≥ 1`. -/
theorem failWitness_pos_of_window (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ) (i : Fin 13)
    (hq : 0 < q) (hp : 0 < p) (hvpos : 0 < v i) (hwin : (q : ℤ) ≤ 14 * v i)
    (hbad : ¬ inBand v q p i) :
    1 ≤ failWitness v q p i := by
  have hw := bad_at_witness v q p i hq hbad
  set n : ℤ := failWitness v q p i with hn
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hpZ : (1 : ℤ) ≤ (p : ℤ) := by exact_mod_cast hp
  have hd : v i * (p : ℤ) - n * (q : ℤ) ≤ |v i * (p : ℤ) - n * (q : ℤ)| := le_abs_self _
  have hlow : 14 * (v i * (p : ℤ)) - 14 * (n * (q : ℤ)) < (q : ℤ) := by
    nlinarith [hw, hd]
  have hvp : (q : ℤ) * (p : ℤ) ≤ 14 * (v i * (p : ℤ)) := by
    have := mul_le_mul_of_nonneg_right hwin (by linarith : (0:ℤ) ≤ (p : ℤ))
    nlinarith [this]
  have hnq : (q : ℤ) * (p : ℤ) - (q : ℤ) < 14 * (n * (q : ℤ)) := by linarith
  by_contra hle
  push Not at hle
  have hn0 : n ≤ 0 := by omega
  have h1 : 14 * (n * (q : ℤ)) ≤ 0 := by
    have : n * (q : ℤ) ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hn0 hqZ.le
    linarith
  have h2 : (0 : ℤ) ≤ (q : ℤ) * (p : ℤ) - (q : ℤ) := by nlinarith
  linarith

/-- The ceiling: for `p < q`, the witness is at most `v i`. -/
theorem failWitness_le (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ) (i : Fin 13)
    (hq : 0 < q) (hpq : p < q) (hvpos : 0 < v i)
    (hbad : ¬ inBand v q p i) :
    failWitness v q p i ≤ v i := by
  have hw := bad_at_witness v q p i hq hbad
  set n : ℤ := failWitness v q p i with hn
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hpqZ : (p : ℤ) < (q : ℤ) := by exact_mod_cast hpq
  have hd : n * (q : ℤ) - v i * (p : ℤ) ≤ |v i * (p : ℤ) - n * (q : ℤ)| := by
    rw [abs_sub_comm]
    exact le_abs_self _
  have hhigh : 14 * (n * (q : ℤ)) - 14 * (v i * (p : ℤ)) < (q : ℤ) := by
    nlinarith [hw, hd]
  by_contra hgt
  push Not at hgt
  have h1 : v i + 1 ≤ n := by omega
  have h2 : (v i + 1) * (q : ℤ) ≤ n * (q : ℤ) := mul_le_mul_of_nonneg_right h1 hqZ.le
  have h3 : v i * (p : ℤ) ≤ v i * ((q : ℤ) - 1) :=
    mul_le_mul_of_nonneg_left (by omega) hvpos.le
  nlinarith [hhigh, h2, h3, hqZ, hvpos]

/-- **THE RIGIDITY ENGINE**: simultaneously bad runners with ratio in `[3, 13]`
and a positive lower witness force `n_j ≥ 3·n_i` — the witness vector inherits
the ladder. -/
theorem witness_ladder (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ) (i j : Fin 13)
    (hq : 0 < q) (hvi : 0 < v i) (hlad : 3 * v i ≤ v j) (hcomp : v j ≤ 13 * v i)
    (hbadi : ¬ inBand v q p i) (hbadj : ¬ inBand v q p j)
    (hni : 1 ≤ failWitness v q p i) :
    3 * failWitness v q p i ≤ failWitness v q p j := by
  have hvj : 0 < v j := by linarith
  have hcon := seven_overlap_pair_constraint v q p i j hq hvi.ne' hvj.ne'
    hbadi hbadj
  set ni : ℤ := failWitness v q p i with hni'
  set nj : ℤ := failWitness v q p j with hnj'
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  -- divide the constraint by q: 14·|v_i·n_j − v_j·n_i| < v_i + v_j
  have habs : 14 * |v i * nj - v j * ni| < v i + v j := by
    have h1 : 14 * (q : ℤ) * |v i * nj - v j * ni| < (|v i| + |v j|) * (q : ℤ) := hcon
    rw [abs_of_pos hvi, abs_of_pos hvj] at h1
    by_contra hge
    push Not at hge
    have h2 : (v i + v j) * (q : ℤ) ≤ (14 * |v i * nj - v j * ni|) * (q : ℤ) :=
      mul_le_mul_of_nonneg_right hge hqZ.le
    nlinarith [h1, h2]
  have hd : v j * ni - v i * nj ≤ |v i * nj - v j * ni| := by
    rw [abs_sub_comm]
    exact le_abs_self _
  -- 14·v_i·n_j > 14·v_j·n_i − (v_i + v_j) ≥ 42·v_i·n_i − 14·v_i (using v_j ≤ 13 v_i)
  have hkey : 14 * (v i * nj) > 42 * (v i * ni) - 14 * v i := by
    have h1 : 14 * (v j * ni) - 14 * (v i * nj) < v i + v j := by nlinarith [habs, hd]
    have h2 : 42 * (v i * ni) ≤ 14 * (v j * ni) := by
      have := mul_le_mul_of_nonneg_right hlad (by linarith : (0:ℤ) ≤ ni)
      nlinarith [this]
    have h3 : v i + v j ≤ 14 * v i := by linarith
    linarith
  -- divide by v_i: n_j > 3·n_i − 1, hence ≥ 3·n_i
  by_contra hlt
  push Not at hlt
  have h1 : nj ≤ 3 * ni - 1 := by omega
  have h2 : v i * nj ≤ v i * (3 * ni - 1) := mul_le_mul_of_nonneg_left h1 hvi.le
  nlinarith [hkey, h2, hvi]

/-- **The geometric chain**: a list of simultaneously bad runners with consecutive
ratios in `[3, 13]` (indexed form; build it from any chain structure) forces the
witness of position `k` to exceed `3^k` times the first witness. -/
theorem witness_chain_geometric (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ)
    (idx : List (Fin 13)) (hlen : 1 ≤ idx.length)
    (hpos : ∀ i ∈ idx, 0 < v i)
    (hbad : ∀ i ∈ idx, ¬ inBand v q p i)
    (hchain : ∀ (m : ℕ) (hm1 : m + 1 < idx.length),
      3 * v (idx.get ⟨m, by omega⟩) ≤ v (idx.get ⟨m + 1, hm1⟩) ∧
      v (idx.get ⟨m + 1, hm1⟩) ≤ 13 * v (idx.get ⟨m, by omega⟩))
    (hq : 0 < q)
    (hfirst : ∀ h : 0 < idx.length, 1 ≤ failWitness v q p (idx.get ⟨0, h⟩)) :
    ∀ (k : ℕ) (hk : k < idx.length),
      3 ^ k * failWitness v q p (idx.get ⟨0, by omega⟩)
        ≤ failWitness v q p (idx.get ⟨k, hk⟩) := by
  intro k
  induction k with
  | zero =>
    intro hk
    simp
  | succ m ih =>
    intro hk
    have hm : m < idx.length := by omega
    have hstep := ih hm
    -- the chain link between positions m and m+1
    have hlink : 3 * v (idx.get ⟨m, hm⟩) ≤ v (idx.get ⟨m + 1, hk⟩) ∧
        v (idx.get ⟨m + 1, hk⟩) ≤ 13 * v (idx.get ⟨m, hm⟩) := hchain m hk
    have hni : 1 ≤ failWitness v q p (idx.get ⟨m, hm⟩) := by
      have h0 := hfirst (by omega)
      have h3 : (1 : ℤ) ≤ 3 ^ m * failWitness v q p (idx.get ⟨0, by omega⟩) := by
        have hpow : (1 : ℤ) ≤ 3 ^ m := one_le_pow₀ (by norm_num)
        nlinarith [h0, hpow]
      linarith [hstep]
    have hlad := witness_ladder v q p (idx.get ⟨m, hm⟩) (idx.get ⟨m + 1, hk⟩) hq
      (hpos _ (List.get_mem _ _)) hlink.1 hlink.2
      (hbad _ (List.get_mem _ _)) (hbad _ (List.get_mem _ _)) hni
    have hpow3 : (3 : ℤ) ^ (m + 1) * failWitness v q p (idx.get ⟨0, by omega⟩)
        = 3 * (3 ^ m * failWitness v q p (idx.get ⟨0, by omega⟩)) := by ring
    calc (3 : ℤ) ^ (m + 1) * failWitness v q p (idx.get ⟨0, by omega⟩)
        = 3 * (3 ^ m * failWitness v q p (idx.get ⟨0, by omega⟩)) := hpow3
      _ ≤ 3 * failWitness v q p (idx.get ⟨m, hm⟩) := by linarith [hstep]
      _ ≤ failWitness v q p (idx.get ⟨m + 1, hk⟩) := hlad

/-- **THE PIN**: witness `n ≥ N` forces `(14·N − 1)·q < 14·v·p` — deep overlaps
live far from zero, in exact integers. -/
theorem witness_pin (v : Fin 13 → ℤ) (q : ℕ) (p : ℕ) (i : Fin 13) (N : ℤ)
    (hq : 0 < q) (hbad : ¬ inBand v q p i)
    (hN : N ≤ failWitness v q p i) :
    (14 * N - 1) * (q : ℤ) < 14 * (v i * (p : ℤ)) := by
  have hw := bad_at_witness v q p i hq hbad
  set n : ℤ := failWitness v q p i with hn
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hd : n * (q : ℤ) - v i * (p : ℤ) ≤ |v i * (p : ℤ) - n * (q : ℤ)| := by
    rw [abs_sub_comm]
    exact le_abs_self _
  have h1 : 14 * (n * (q : ℤ)) - 14 * (v i * (p : ℤ)) < (q : ℤ) := by
    nlinarith [hw, hd]
  have h2 : (14 * N) * (q : ℤ) ≤ 14 * (n * (q : ℤ)) := by
    have := mul_le_mul_of_nonneg_right hN hqZ.le
    nlinarith [this]
  nlinarith [h1, h2]

/-! ## Axiom audit -/
#print axioms failWitness_pos_of_window
#print axioms failWitness_le
#print axioms witness_ladder
#print axioms witness_chain_geometric
#print axioms witness_pin

end LRC14Concrete
end LonelyRunner
