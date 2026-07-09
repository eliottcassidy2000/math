/-
  TournamentH7.LRCLedgerConsumer — the pair-sum liveness consumer (kind-pasteur-2026-07-09-S117).

  mac-mini's THM-668 addendum gives *liveness certificates* for a pair-sum ruler `q = v_i+v_j`: the
  C1 gcd-exact ledger (`|B_l| = g(2⌊m/g⌋+1)`, `g = gcd(v_l mod q, q)`, `m = ⌈q/14⌉−1`, ± class merge,
  now THM-672/674), C2 divisor descent, C3 six-pair prime — each **bounds the number of BLOCKED
  multipliers** `p` (those for which some runner `v_l·p mod q` leaves the band `[q/14,13q/14]`).  A
  ruler is *live* when that blocked count is `< q−1`, so a good multiplier survives.

  This file supplies the shared **consumer**: the union-bound step that turns any such blocked-count
  bound into loneliness.  Over the nonzero multipliers `{1,…,q−1}` (indexed by `Finset.range (q−1)`,
  kept in `ℕ` so the count is `native_decide`-computable), if fewer than `q−1` fail to fire the ruler,
  some `p` fires — every residue lands in the band — and `mreach_ge_of_pairsum_band` (kps-S114,
  THM-668) gives `Mreach ≥ 1/14`.  The certificates (mac-mini's ledger / klein's signed box) discharge
  the count; this is the half that reaches the `Mreach` socket — so a liveness census becomes a
  loneliness theorem.  Self-contained apart from `LRCPairSumDispatch`.
-/
import Mathlib
import TournamentH7.LRCPairSumDispatch

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The multiplier `p` **fires** the ruler `q`: every runner's residue is in the safe band
`[q/14, 13q/14]` — i.e. `p/q` is a lonely instant. -/
def fires (v : Fin 13 → ℤ) (q p : ℤ) : Prop :=
  ∀ i, q ≤ 14 * ((v i * p) % q) ∧ 14 * ((v i * p) % q) ≤ 13 * q

instance (v : Fin 13 → ℤ) (q p : ℤ) : Decidable (fires v q p) :=
  Fintype.decidableForallFintype

/-- **The pair-sum liveness consumer (union bound → loneliness).**  The nonzero multipliers `{1,…,q−1}`
are `{(p+1) : p ∈ range (q−1)}`.  If the number of them that FAIL to fire the ruler `q` is strictly
below the total `q−1`, then some multiplier fires and `Mreach ≥ 1/14`.  This is the half of every
liveness certificate (C0/C1/C2/C3) that turns a blocked-count bound into loneliness through
`mreach_ge_of_pairsum_band`; the gcd-exact ledger supplies the count.  The count is over `Finset.range`
(computable), so a concrete liveness census is `native_decide`-checkable (see `demo_c1_lonely`). -/
theorem mreach_ge_of_blocked_lt (v : Fin 13 → ℤ) (q : ℤ) (hq : 0 < q) (N : ℕ) (hN : (N : ℤ) = q - 1)
    (hcount : ((Finset.range N).filter (fun p : ℕ => ¬ fires v q ((p : ℤ) + 1))).card < N) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  set B := (Finset.range N).filter (fun p : ℕ => ¬ fires v q ((p : ℤ) + 1)) with hB
  have hBsub : B ⊆ Finset.range N := Finset.filter_subset _ _
  have hne : B ≠ Finset.range N := by
    intro h; rw [h, Finset.card_range] at hcount; exact absurd hcount (lt_irrefl _)
  have hss : B ⊂ Finset.range N := Finset.ssubset_iff_subset_ne.mpr ⟨hBsub, hne⟩
  obtain ⟨p, hpR, hpB⟩ := Finset.exists_of_ssubset hss
  have hfires : fires v q ((p : ℤ) + 1) := by
    by_contra h
    exact hpB (Finset.mem_filter.mpr ⟨hpR, h⟩)
  exact mreach_ge_of_pairsum_band v ((p : ℤ) + 1) q hq hfires

/-- **Consumer, phrased as an existence hypothesis.**  If *some* multiplier fires the ruler (the
certificate's conclusion stated directly), then `Mreach ≥ 1/14` — the trivial socket. -/
theorem mreach_ge_of_live_ruler (v : Fin 13 → ℤ) (q : ℤ) (hq : 0 < q)
    (hlive : ∃ p : ℤ, fires v q p) : (1 : ℝ) / 14 ≤ Mreach v := by
  obtain ⟨p, hp⟩ := hlive
  exact mreach_ge_of_pairsum_band v p q hq hp

/-- **Demonstration: certificate census → loneliness theorem, machine-checked.**  For the covering set
`{1,2,3,4,5,6,8,10,11,12,13,14,18}` at the pair-sum ruler `q = 16` (`= 2+14`), `native_decide` verifies
that only `13 < 15 = q−1` multipliers are blocked — the C1 liveness condition — and the consumer turns
that count directly into `Mreach ≥ 1/14`, WITHOUT exhibiting the firing multiplier.  This is the shape a
general ledger bound plugs into. -/
theorem demo_c1_lonely :
    (1 : ℝ) / 14 ≤ Mreach ![1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 14, 18] := by
  refine mreach_ge_of_blocked_lt _ 16 (by norm_num) 15 (by norm_num) ?_
  native_decide

end LRC14Concrete
end LonelyRunner
