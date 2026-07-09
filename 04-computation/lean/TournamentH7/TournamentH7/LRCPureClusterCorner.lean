/-
  TournamentH7.LRCPureClusterCorner — the P = ∅ (pure-cluster) corner of the realization node
  is ALREADY CLOSED, by subsumption (death-star-2026-07-09-S1, HYP-5710 part 1).

  CONTEXT.  The fleet's sole remaining Part-A node is the Kronecker/equidistribution realization
  (klein-S207, kps-S112): turn a positive-measure good event of the co-offset cluster into a real
  lonely time.  Four concurrent claims attack it (boxeph HYP-5708 LEM-014 P-separated wide regime;
  monad HYP-5717 clamped grid-port + P∪L criterion; mac-mini HYP-5720 witness slack; this session's
  HYP-5710 pure-cluster continuum sweep).  This file records that for PURE clusters — every runner
  in ONE cluster at the ruler scale, P = ∅ — **no realization machinery is needed at all**:

      co-offsets e_i = Vmax − v_i ∈ [0, D]  with  13·D ≤ 12·Vmax
        ⟹  all speeds lie in the band [Vmax − D, Vmax] of ratio ≤ 13
        ⟹  lonely at the EXPLICIT witness  t = 1/(2·Vmax − D)     (kps-S28 `spread13_lonely`)
        ⟹  Mreach v ≥ 1/14                                        (kps-S10 `Mreach_ge_of_lonely`).

  I claimed HYP-5710 to build this corner via an endpoint-interval + IVT sweep (and klein-S205's
  `Mreach_ge_of_driftGap` instantiates at j = 1, a = D/V, g = 1 − D/V for Vmax ≳ 2.76·D); the
  claim's own verify-first step found BOTH subsumed by the sharp band window `spread13_lonely`
  (kps-2026-07-03-S28), with a BETTER constant: Vmax ≥ 13·D/12 ≈ 1.083·D, and no drift, no gap,
  no measure.  mac-mini-S64 pointed here ("r ≤ 13 already covered by kps-S28"); this file makes the
  pointer a theorem in the co-offset vocabulary and the `Mreach` socket the realization work uses.

  DELINEATION (what remains open).  The residual of the realization node is exactly the
  speeds-ratio > 13 case: max|v| > 13·min|v|, i.e. a genuine SCALE GAP in the speed spectrum —
  the multi-scale P ∪ L composition (boxeph LEM-014, monad HYP-5717, the coarse reduction).
  A pure cluster can never be that case; effort spent realizing good periods for P = ∅ clusters
  (at any diameter D, any Vmax ≥ 13·D/12) is spent on a closed corner.

  Kernel-pure target: no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCSpread13
import TournamentH7.LRC14Dispatch

namespace LonelyRunner
namespace LRC14

/-- **Pure clusters are lonely at an explicit witness.**  If every co-offset
`e i = Vmax − v i` lies in `[0, D]` and `13·D ≤ 12·Vmax`, the speeds all lie in the
band `[Vmax − D, Vmax]` of ratio `≤ 13`, so `spread13_lonely` (kps-S28) applies with
the explicit witness `t = 1/((Vmax − D) + Vmax)`.  No good period, no drift, no
density floor: the P = ∅ corner of the realization node needs none of them. -/
theorem pure_cluster_lonely (v e : Fin 13 → ℤ) (Vmax D : ℤ)
    (hbind : ∀ i, v i = Vmax - e i)
    (he0 : ∀ i, 0 ≤ e i) (heD : ∀ i, e i ≤ D)
    (hV : 0 < Vmax) (hband : 13 * D ≤ 12 * Vmax) :
    Lonely 14 v (1 / (((Vmax - D : ℤ) : ℝ) + ((Vmax : ℤ) : ℝ))) := by
  have hDltV : D < Vmax := by linarith
  have hlo : ∀ i, Vmax - D ≤ |v i| := by
    intro i
    have h1 : Vmax - D ≤ v i := by
      have := heD i
      rw [hbind i]; omega
    rw [abs_of_pos (by omega)]
    exact h1
  have hhi : ∀ i, |v i| ≤ Vmax := by
    intro i
    have h1 : v i ≤ Vmax := by
      have := he0 i
      rw [hbind i]; omega
    have h2 : 0 < v i := by
      have := heD i
      rw [hbind i]; omega
    rw [abs_of_pos h2]
    exact h1
  exact spread13_lonely v (Vmax - D) Vmax (by omega) hlo hhi (by omega)

/-- **The pure-cluster corner of the realization node, in the reach socket.**
Co-offsets in `[0, D]` with `13·D ≤ 12·Vmax` force `Mreach v ≥ 1/14` — the exact
conclusion shape of `hpartA`/`hembed`, with no good-period input.  Composition
work on the Kronecker realization (P ∪ L) can consume this as the L-leg base
case and restrict attention to speeds-ratio > 13 shapes. -/
theorem mreach_ge_of_pure_cluster (v e : Fin 13 → ℤ) (Vmax D : ℤ)
    (hbind : ∀ i, v i = Vmax - e i)
    (he0 : ∀ i, 0 ≤ e i) (heD : ∀ i, e i ≤ D)
    (hV : 0 < Vmax) (hband : 13 * D ≤ 12 * Vmax) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v :=
  Mreach_ge_of_lonely ⟨_, pure_cluster_lonely v e Vmax D hbind he0 heD hV hband⟩

/-- **Band reading, no co-offset bookkeeping**: all speeds in `[a, 13a]` (`a > 0`) —
"no scale gap anywhere" — forces the reach.  Contrapositive of the delineation:
a family on which the realization node is genuinely open must have a speed pair
of ratio > 13. -/
theorem mreach_ge_of_ratio_band (v : Fin 13 → ℤ) (a : ℤ) (ha : 0 < a)
    (hlo : ∀ i, a ≤ |v i|) (hhi : ∀ i, |v i| ≤ 13 * a) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v :=
  Mreach_ge_of_lonely ⟨_, spread13_lonely v a (13 * a) ha hlo hhi (by omega)⟩

/-! ## Axiom audit -/
#print axioms pure_cluster_lonely
#print axioms mreach_ge_of_pure_cluster
#print axioms mreach_ge_of_ratio_band

end LRC14
end LonelyRunner
