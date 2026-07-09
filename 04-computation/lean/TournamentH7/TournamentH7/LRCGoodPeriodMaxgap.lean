/-
  TournamentH7.LRCGoodPeriodMaxgap — the maxgap good-period certificate (klein-2026-07-09-S200)

  The dissociated branch of THM-527-A closes on the MAXGAP margin (kps-S94/LEM-013), NOT the
  arc-count (klein-S200: 7-structured co-offsets — differences ≡ 0 mod 7 — spike `#arcs` while
  leaving the best maxgap large; so the arc-count is the wrong invariant, the maxgap is the right
  one).  This file gives the DECIDABLE maxgap good-period predicate and a `native_decide`
  certificate on the exact cluster where the arc-count route fails.

  A good period `j` has `maxgap{e·j mod Vmax : e ∈ E} > Vmax/7`, i.e. `7 · maxCircGap > Vmax`
  (all integer arithmetic — no rationals, no equidistribution).  This is the LEM-013 mechanism,
  fully computable; it is the template for the finite-check nodes.
-/
import Mathlib

namespace LRC14

/-- The maximal circular gap of the residues `{e · j mod Vmax : e ∈ E}` on `ZMod Vmax`,
returned as a natural number in `[0, Vmax]` (units of `1/Vmax`).  Sort the residues, append the
first `+ Vmax` to close the circle, and take the largest consecutive difference. -/
def maxCircGap (E : List ℕ) (Vmax j : ℕ) : ℕ :=
  let ps := (E.map (fun e => e * j % Vmax)).mergeSort (· ≤ ·)
  match ps with
  | []        => Vmax
  | p0 :: _   =>
    let cyc := ps ++ [p0 + Vmax]
    (List.zipWith (fun a b => b - a) cyc cyc.tail).foldl max 0

/-- `j` is a **good period** for cluster `E` and ruler `Vmax` when the phases leave a circular
gap `> Vmax/7`, i.e. `7 · maxCircGap > Vmax`.  (THM-527: this forces `M(S) ≥ 1/14`.) -/
def IsGoodPeriod (E : List ℕ) (Vmax j : ℕ) : Prop := Vmax < 7 * maxCircGap E Vmax j

instance (E : List ℕ) (Vmax j : ℕ) : Decidable (IsGoodPeriod E Vmax j) :=
  inferInstanceAs (Decidable (Vmax < 7 * maxCircGap E Vmax j))

/-- A cluster `E` **has a good period** on ruler `Vmax` if some `j ∈ {1,…,Vmax−1}` is good.
Decidable, hence `native_decide`-checkable — the LEM-013 existence, made a finite computation. -/
def HasGoodPeriod (E : List ℕ) (Vmax : ℕ) : Prop :=
  ∃ j ∈ Finset.Ioo 0 Vmax, IsGoodPeriod E Vmax j

instance (E : List ℕ) (Vmax : ℕ) : Decidable (HasGoodPeriod E Vmax) :=
  inferInstanceAs (Decidable (∃ j ∈ Finset.Ioo 0 Vmax, _))

/-- **The maxgap route works where the arc-count route fails (klein-S200).**
The primitive, hard, dissociated (longest-AP = 4) cluster
`E = {0,7,14,21,26,29,37,44,51,58,67,75,82}` has `#arcs/spread = 0.878 > D3(E) = 0.629`, so
mac-mini's Route-(c) arc-count certificate `c < D3` FAILS (MISTAKE-128) — the four co-offsets
`≡ 0 mod 7` spike the arc-count.  Yet a good period EXISTS (LEM-013): here `native_decide`
certifies it directly on the ruler `Vmax = 91 = 7·13` where the cluster is hardest. -/
theorem worst7Struct_hasGoodPeriod :
    HasGoodPeriod [0,7,14,21,26,29,37,44,51,58,67,75,82] 91 := by
  native_decide

/-- The severe large-spread 7-structured cluster (spread 406, `#arcs/spread = 0.768`, again
`c > D3`) also has a good period, on its hard ruler — the arc-count spike does not obstruct
existence. -/
theorem worst7StructLarge_hasGoodPeriod :
    HasGoodPeriod [0,7,63,126,151,189,252,305,315,362,378,385,406] 458 := by
  native_decide

end LRC14
