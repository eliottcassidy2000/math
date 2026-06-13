/-
  TournamentH7.SCFraction — SC tiling counts and fraction (THM-330 corollary)

  ─── What this module records ──────────────────────────────────────────
  Project canon (THM-330 Corollary 3): the number of SC tilings on n
  vertices among all 2^{C(n-1, 2)} tilings.

  The values are:
      n |  m=C(n-1,2)  |  Total tilings  |  SC tilings  |  Non-SC
      3 |       1      |        2        |       1      |    1
      4 |       3      |        8        |       5      |    3
      5 |       6      |       64        |      50      |   14
      6 |      10      |     1024        |     903      |  121
      7 |      15      |    32768        |   30773      | 1995
      8 |      21      |  2097152        | 2032504      | 64648

  ─── HYP-1739 sequence ────────────────────────────────────────────────
  The non-SC counts 1, 3, 14, 121, 1995, 64648 form a sequence whose
  OEIS lookup gives no canonical match — this is a *new* sequence (or
  at least one not recorded in OEIS as of 2026-05-27).

  ─── Lean encoding ────────────────────────────────────────────────────
  We axiomatise the SC tiling counts and provide trivial derived
  identities (Total = SC + NonSC).
-/

import TournamentH7.Basic
import TournamentH7.StaircaseModel
import Mathlib.Data.Nat.Choose.Basic

namespace Tournament

/-! ### SC tiling counts -/

/-- `SCtilings n` = number of SC tilings on n vertices (project canon). -/
opaque SCtilings : ℕ → ℕ

@[simp] axiom SCtilings_3 : SCtilings 3 = 1
@[simp] axiom SCtilings_4 : SCtilings 4 = 5
@[simp] axiom SCtilings_5 : SCtilings 5 = 50
@[simp] axiom SCtilings_6 : SCtilings 6 = 903
@[simp] axiom SCtilings_7 : SCtilings 7 = 30773
@[simp] axiom SCtilings_8 : SCtilings 8 = 2032504

/-- `NonSCtilings n` = total tilings minus SC tilings. -/
opaque NonSCtilings : ℕ → ℕ

@[simp] axiom NonSCtilings_3 : NonSCtilings 3 = 1
@[simp] axiom NonSCtilings_4 : NonSCtilings 4 = 3
@[simp] axiom NonSCtilings_5 : NonSCtilings 5 = 14
@[simp] axiom NonSCtilings_6 : NonSCtilings 6 = 121
@[simp] axiom NonSCtilings_7 : NonSCtilings 7 = 1995
@[simp] axiom NonSCtilings_8 : NonSCtilings 8 = 64648

/-! ### Total tilings = 2^{C(n-1, 2)} -/

/-- Total tilings on n vertices. -/
def totalTilings (n : ℕ) : ℕ := 2 ^ (Nat.choose (n - 1) 2)

example : totalTilings 3 = 2 := by unfold totalTilings; decide
example : totalTilings 4 = 8 := by unfold totalTilings; decide
example : totalTilings 5 = 64 := by unfold totalTilings; decide
example : totalTilings 6 = 1024 := by unfold totalTilings; decide
example : totalTilings 7 = 32768 := by unfold totalTilings; decide

/-! ### Partition identity: Total = SC + NonSC -/

/-- **Axiom (partition).**  Every tiling is either SC or non-SC. -/
axiom totalTilings_eq_partition (n : ℕ) (hn : 3 ≤ n) :
    totalTilings n = SCtilings n + NonSCtilings n

/-! ### Some derived facts -/

example : SCtilings 3 + NonSCtilings 3 = totalTilings 3 := by
  rw [SCtilings_3, NonSCtilings_3]
  show 2 = totalTilings 3
  unfold totalTilings; decide

example : SCtilings 4 + NonSCtilings 4 = totalTilings 4 := by
  rw [SCtilings_4, NonSCtilings_4]
  show 8 = totalTilings 4
  unfold totalTilings; decide

example : SCtilings 5 + NonSCtilings 5 = totalTilings 5 := by
  rw [SCtilings_5, NonSCtilings_5]
  show 64 = totalTilings 5
  unfold totalTilings; decide

/-! ### Asymptotic fact -/

/-- **Observation.**  The ratio `SCtilings n / totalTilings n` tends to 1
    as n → ∞.  This is the project's "SC fraction grows toward 2^m"
    observation. -/
example : 100 * SCtilings 8 > 96 * totalTilings 8 := by
  rw [SCtilings_8]
  show 100 * 2032504 > 96 * totalTilings 8
  unfold totalTilings
  -- 100 * 2,032,504 = 203,250,400; 96 * 2,097,152 = 201,326,592. So yes.
  decide

end Tournament
