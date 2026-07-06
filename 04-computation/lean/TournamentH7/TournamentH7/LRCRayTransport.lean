/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S98)
-/
import TournamentH7.LRCLiftRowsL7

/-!
# Ray transport: mod-q witnesses are CRT-ray invariant (HYP-4266)

The arithmetic heart of the periodicity reduction (mac-mini S55), stated and proved:
a witness at `t = a/q` has margins depending only on the speeds MOD `q` — so shifting
any speed by a multiple of `q` (in particular, moving any distance along a CRT-frozen
ray whose period `q` divides, INCLUDING the dilation rays `S·p` with `S ≡ 1 mod lcm`)
changes nothing. One single-period verification therefore covers the whole ray: this is
why witness spectroscopy at `q ≤ Q₀` (with `q` dividing the ray period) is
height-independent, and it is the companion of the S97 two-band transport (which
manufactures scale-tracking witnesses when no small-`q` witness exists).

Three levels: the real margin (`margin_ray_invariant`), the kernel gate
(`speedOK13_ray`, via S82's `speedOK13_congr`), and the whole family
(`family_margin_ray_invariant`).
-/

namespace LonelyRunner
namespace RayTransport

open LiftRowsL7 KernelGate13

/-- **Mod-q margin invariance**: at `t = a/q`, shifting the speed by `M·k` with `q ∣ M`
preserves every integer-distance margin exactly. -/
theorem margin_ray_invariant (v M k a : ℤ) (q : ℕ) (hq : 0 < q) (hM : (q : ℤ) ∣ M)
    (c : ℝ) (h : ∀ m : ℤ, c ≤ |(v : ℝ) * ((a : ℝ) / q) - m|) :
    ∀ m : ℤ, c ≤ |((v + M * k : ℤ) : ℝ) * ((a : ℝ) / q) - m| := by
  intro m
  obtain ⟨M', rfl⟩ := hM
  have hqR : ((q : ℕ) : ℝ) ≠ 0 := by
    exact_mod_cast hq.ne'
  have key : ((v + (q : ℤ) * M' * k : ℤ) : ℝ) * ((a : ℝ) / q)
      = (v : ℝ) * ((a : ℝ) / q) + ((M' * k * a : ℤ) : ℝ) := by
    push_cast
    field_simp
  rw [key]
  have h2 := h (m - M' * k * a)
  convert h2 using 2
  push_cast
  ring

/-- **Gate-level ray invariance**: the strict kernel gate `speedOK13` transports along
rays — a checked row at one point of a CRT ray is a row at every point. -/
theorem speedOK13_ray (s M k num : ℤ) (den : ℕ) (hM : (den : ℤ) ∣ M)
    (h : speedOK13 s num den) : speedOK13 (s + M * k) num den := by
  apply speedOK13_congr ?_ h
  obtain ⟨M', rfl⟩ := hM
  conv_lhs => rw [show s + (den : ℤ) * M' * k = s + (den : ℤ) * (M' * k) by ring]
  first
    | exact Int.add_mul_emod_self_left
    | simp [Int.add_mul_emod_self_left]

/-- **Family-level ray invariance**: a whole-family margin at `a/q` survives shifting
each speed by its own `q`-multiple — one period certifies the entire affine ray. -/
theorem family_margin_ray_invariant {ι : Type*} (v Mv : ι → ℤ) (k a : ℤ) (q : ℕ)
    (hq : 0 < q) (hM : ∀ i, (q : ℤ) ∣ Mv i) (c : ℝ)
    (h : ∀ i, ∀ m : ℤ, c ≤ |(v i : ℝ) * ((a : ℝ) / q) - m|) :
    ∀ i, ∀ m : ℤ, c ≤ |((v i + Mv i * k : ℤ) : ℝ) * ((a : ℝ) / q) - m| :=
  fun i => margin_ray_invariant (v i) (Mv i) k a q hq (hM i) c (h i)

/-- **The census bridge**: a family's margin at `a/q` equals its RESIDUE family's margin
there.  Given a witness for the reduced family `r` (bounded: speeds `< q`) with each
`v i ≡ r i (mod q)`, the ORIGINAL family — at ANY height — inherits the margin.  This is
the periodicity reduction made consumable: census the bounded residue families (mac-mini
/ kps lane), transport to every height for free. -/
theorem margin_of_residue_witness {ι : Type*} (v r : ι → ℤ) (a : ℤ) (q : ℕ) (hq : 0 < q)
    (hres : ∀ i, (q : ℤ) ∣ (v i - r i)) (c : ℝ)
    (h : ∀ i, ∀ m : ℤ, c ≤ |(r i : ℝ) * ((a : ℝ) / q) - m|) :
    ∀ i, ∀ m : ℤ, c ≤ |(v i : ℝ) * ((a : ℝ) / q) - m| := by
  intro i m
  have hshift := margin_ray_invariant (r i) (v i - r i) 1 a q hq (hres i) c (h i)
  have heq : r i + (v i - r i) * 1 = v i := by ring
  simpa [heq] using hshift m

#print axioms margin_ray_invariant
#print axioms speedOK13_ray
#print axioms family_margin_ray_invariant
#print axioms margin_of_residue_witness

end RayTransport
end LonelyRunner
