/-
  TournamentH7.LRCEvenDescent — THE TOWER STEP OF THE CONFINEMENT DESCENT, FORMALIZED
  (klein-2026-07-03-S121, HYP-4069).

  THM-612 (mac-mini) reduces the covering-min lower-bound RIGIDITY to CONFINEMENT: a primitive
  tight family `S` (`M(S)=1/n`, hiding denominator `q*=14m`) has its `m`-divisible part
  `E = m·U` carry the local geometry (Lemma B: `f_S(t)=g_U(m·t)` near `t*`), so `U` is `n`-lonely
  at `m·t*` and (Corollary, via LRC≤13) strictly loose.  The DESCENT engine is here made a clean,
  citable Lean lemma: **the `m`-divisible sub-family of an `n`-lonely family descends to `U`
  `n`-lonely one scale down at `m·t`.**  Built on `LonelyRunner.lonely_scale`.

  SCOPE (HONEST — this does NOT close the rigidity): this is the tower STEP (elementary, proved).
  The confinement's open residual is *bounding `v_max(U)`* — and mac-mini-S34 showed that residual
  IS the covering-min itself (`q*>14` tight ⟹ the family covers every `q≤14` ⟹ it is a primitive
  covering family with `M=1/14`, forbidden by the covering-min HYP-4060).  So the even-part descent
  is NOT an independent lever: bounding the even part is LRC(14)-equivalent.  opus-S62 has it only at
  EVIDENCE (the `f=2` tightening gap grows with `u_max`).  This file formalizes the descent's
  reduction step; the LRC-hard core (the covering-min) is untouched and open.
-/
import TournamentH7.LonelyRunner

namespace LonelyRunner.EvenDescent

/-- A sub-family of an `n`-lonely family is `n`-lonely (loneliness is downward-closed on the
  index set): reindex by any `e : κ → ι`. -/
theorem lonely_subfamily {ι κ : Type*} (n : ℕ) (v : ι → ℤ) (t : ℝ) (e : κ → ι)
    (h : Lonely n v t) : Lonely n (v ∘ e) t :=
  fun j m => h (e j) m

/-- **Scale descent.**  If the dilated family `m·U` is `n`-lonely at `t`, then `U` is `n`-lonely
  at `m·t` (`m ≠ 0`).  Because `‖(m u) t‖ = ‖u (m t)‖`; the packaged form of `lonely_scale`. -/
theorem even_part_descends {κ : Type*} (n : ℕ) (m : ℤ) (hm : m ≠ 0) (U : κ → ℤ) (t : ℝ)
    (h : Lonely n (fun j => m * U j) t) : Lonely n U ((m : ℝ) * t) := by
  have hmR : (m : ℝ) ≠ 0 := by exact_mod_cast hm
  have he : ((m : ℝ) * t) / (m : ℝ) = t := by field_simp
  have key := lonely_scale n U ((m : ℝ) * t) m hm
  rw [he] at key
  exact key.mp h

/-- **The tower step (THM-612 Lemma B, loneliness form).**  Let the full family `v` be `n`-lonely
  at `t` (e.g. `t = t*` of a tight family).  Let `e : κ → ι` pick out an `m`-divisible sub-family
  with `v (e j) = m · U j`.  Then the descended family `U` is `n`-lonely at `m·t`.

  For a tight `S` at `q* = 14 m`, taking `v = S`, `t = t*`, and `e` the `m`-divisible runners
  `E = m·U`, this is exactly "`U` is `14`-lonely at `a/14`" — the reduction one scale down that
  drives the confinement tower.  (By LRC≤13, when `S` is primitive and `m≥2`, `U` then has `≤12`
  runners and is *strictly* loose; the tighteners `F = S∖E` must suppress its super-`1/14` region.
  Bounding `v_max(U)` = the open, LRC-equivalent covering-min piece — NOT proved here.) -/
theorem tower_step {ι κ : Type*} (n : ℕ) (m : ℤ) (hm : m ≠ 0)
    (v : ι → ℤ) (t : ℝ) (h : Lonely n v t)
    (U : κ → ℤ) (e : κ → ι) (huv : ∀ j, v (e j) = m * U j) :
    Lonely n U ((m : ℝ) * t) := by
  have hsub : Lonely n (fun j => m * U j) t := by
    intro j mm
    have hj := h (e j) mm
    rw [huv j] at hj
    exact hj
  exact even_part_descends n m hm U t hsub

#print axioms lonely_subfamily
#print axioms even_part_descends
#print axioms tower_step

end LonelyRunner.EvenDescent
