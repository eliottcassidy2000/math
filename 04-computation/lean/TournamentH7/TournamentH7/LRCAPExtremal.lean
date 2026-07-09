/-
  TournamentH7.LRCAPExtremal — the AP equality extremal and the density-floor → reach wiring
  (kind-pasteur-2026-07-09-S110).

  Two pieces closing the loop on the window finite check (kps-S109):

  1. `mreach_AP_ge` — the arithmetic-progression runner set `v_i = i` (`i = 1..13`) is lonely with
     `Mreach ≥ 1/14`, witnessed at `τ = 1/14`: `nearInt(i/14) = min(i,14−i)/14 ≥ 1/14`.  This is the
     LRC(14) EQUALITY extremal (the AP is the conjectured minimiser, `M = 1/14`; here the `≥` half —
     the loneliness — which is the direction LRC asserts).

  2. `mreach_ge_of_rhoStar_pos` — the density floor `ρ*_{1/7}(E) ≥ m_P = 14249/252252 > 0` (THM-530
     Bonferroni + THM-661, PROVED) WIRES to `Mreach ≥ 1/14` through the reformulation bridge (a
     point of positive good measure is a lonely time).  This is the CONTINUUM route — no ruler grid,
     no slow-fast drift, no bounded window — so it covers the window (and every cluster) once the
     reformulation bridge is supplied.  The bridge is the sole remaining hypothesis, shared with
     Part A but stated at the continuum (kps-S107/S108: run it on the smooth surrogate).

  Builds on `LRCHembedScaleSep` (`le_nearInt_of_forall_int`, kps-S106) and `LRCReachWitness`
  (`Mreach_ge_of_lonely_instant`, kps-S99b).
-/
import Mathlib
import TournamentH7.LRCHembedScaleSep

namespace LonelyRunner
namespace LRC14Concrete

/-- **`nearInt(k/n) ≥ 1/n` when `n ∤ k`.**  A non-integer multiple of `1/n` is at least `1/n` from
every integer (the balanced remainder has absolute value `≥ 1`). -/
theorem one_div_le_nearInt_of_not_dvd (n : ℕ) (hn : 0 < n) (k : ℤ) (h : ¬ (n : ℤ) ∣ k) :
    (1 : ℝ) / (n : ℝ) ≤ nearInt ((k : ℝ) / (n : ℝ)) := by
  apply le_nearInt_of_forall_int
  intro m
  have hnR : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn
  have hne : k - (n : ℤ) * m ≠ 0 := by
    intro hz; exact h ⟨m, by linarith [hz]⟩
  have h1R : (1 : ℝ) ≤ |((k - (n : ℤ) * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast Int.one_le_abs hne
  have hkey : (k : ℝ) / (n : ℝ) - (m : ℝ) = ((k - (n : ℤ) * m : ℤ) : ℝ) / (n : ℝ) := by
    push_cast; field_simp
  rw [hkey, abs_div, abs_of_pos hnR]
  gcongr

/-- **The AP is lonely: `Mreach ≥ 1/14`.**  For `v_i = i` (`i = 1..13`), the time `τ = 1/14` clears
every runner by `≥ 1/14`, so `Mreach v ≥ 1/14`.  The equality `M = 1/14` (the AP is the extremal) —
the `≤` half is the pigeonhole tightness, separate; this is the loneliness `≥` half. -/
theorem mreach_AP_ge :
    (1 : ℝ) / 14 ≤ Mreach (fun i : Fin 13 => (i.val : ℤ) + 1) := by
  apply Mreach_ge_of_lonely_instant
  refine ⟨1 / 14, fun i => ?_⟩
  have hnd : ¬ ((14 : ℕ) : ℤ) ∣ ((i.val : ℤ) + 1) := by
    have h2 : (i.val : ℤ) < 13 := by exact_mod_cast i.isLt
    have h0 : (0 : ℤ) ≤ (i.val : ℤ) := Int.natCast_nonneg _
    omega
  have hkey := one_div_le_nearInt_of_not_dvd 14 (by norm_num) ((i.val : ℤ) + 1) hnd
  have e1 : (((i.val : ℤ) + 1 : ℤ) : ℝ) * (1 / 14)
      = (((i.val : ℤ) + 1 : ℤ) : ℝ) / ((14 : ℕ) : ℝ) := by push_cast; ring
  have e2 : (1 : ℝ) / 14 = (1 : ℝ) / ((14 : ℕ) : ℝ) := by push_cast; ring
  rw [e1, e2]
  exact hkey

/-- **The density floor WIRES to the reach.**  Given the density floor `m_P ≤ ρ*` (THM-661) with
`0 < m_P` (THM-530, `m_P = 14249/252252`), and the reformulation bridge `0 < ρ* → ∃ lonely time`,
the runner set satisfies `Mreach ≥ 1/14`.  This covers the window (and all clusters) at the
continuum — no grid, no drift.  The bridge is the sole remaining Part-A hypothesis. -/
theorem mreach_ge_of_rhoStar_pos (v : Fin 13 → ℤ) (rhoStar mP : ℝ)
    (hfloor : mP ≤ rhoStar) (hmP : 0 < mP)
    (hrefl : 0 < rhoStar → ∃ τ : ℝ, ∀ i, (1 : ℝ) / 14 ≤ nearInt ((v i : ℝ) * τ)) :
    (1 : ℝ) / 14 ≤ Mreach v :=
  Mreach_ge_of_lonely_instant v (hrefl (lt_of_lt_of_le hmP hfloor))

end LRC14Concrete
end LonelyRunner
