# Mathlib PR drafts: PowerSeries gap lemmas from the GMC(2) formalization

*kind-pasteur-2026-07-22-S128c154. Three general `PowerSeries` lemmas the GMC(2) proof needed that are
absent from Mathlib. Each is extracted from the repo, **generalized** past the repo's `[Field F]` to the
minimal typeclass, and **verified to compile** standalone (`lake env lean`, exit 0). Ready to submit.*

## Gap 1 — char-0 converse of `derivativeFun_C` (cleanest; `Mathlib/RingTheory/PowerSeries/Derivative.lean`)

Repo source `GMC2DvdKCharZeroClosing` (Mathlib-only). Generalized `[Field F][CharZero F]` →
`[CommRing R][NoZeroDivisors R][CharZero R]`.

```lean
namespace PowerSeries
variable {R : Type*} [CommRing R] [NoZeroDivisors R] [CharZero R]

theorem coeff_eq_zero_of_derivativeFun_eq_zero {f : R⟦X⟧}
    (hf : derivativeFun f = 0) {n : ℕ} (hn : 1 ≤ n) : coeff n f = 0 := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  have hk := coeff_derivativeFun f k
  rw [hf] at hk; simp only [map_zero] at hk
  have hne : (k + 1 : R) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero k
  exact (mul_eq_zero.mp hk.symm).resolve_right hne

theorem eq_C_of_derivativeFun_eq_zero {f : R⟦X⟧} (hf : derivativeFun f = 0) :
    f = C (constantCoeff f) := by
  ext n
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn; rw [coeff_zero_eq_constantCoeff, constantCoeff_C]
  · rw [coeff_eq_zero_of_derivativeFun_eq_zero hf hn, coeff_C, if_neg (by omega)]

theorem derivativeFun_eq_zero_iff {f : R⟦X⟧} :
    derivativeFun f = 0 ↔ f = C (constantCoeff f) :=
  ⟨eq_C_of_derivativeFun_eq_zero, fun h => by rw [h, derivativeFun_C]⟩
end PowerSeries
```

## Gap 2 — geometric inverse `(1 − C w · X)⁻¹ = ∑ wⁿ Xⁿ`

Repo source `GMC2DvdKFrameExtraction`. General over `[CommRing R]`.

```lean
namespace PowerSeries
variable {R : Type*} [CommRing R]
theorem oneSubCX_mul_geometric (w : R) :
    (1 - C w * X) * mk (fun n => w ^ n) = 1 := by
  rw [sub_mul, one_mul, mul_assoc]
  ext n
  rw [map_sub, coeff_C_mul, coeff_mk]
  cases n with
  | zero => simp [coeff_zero_X_mul]
  | succ k =>
    rw [coeff_succ_X_mul, coeff_mk, coeff_one, if_neg (Nat.succ_ne_zero k), pow_succ]; ring
end PowerSeries
```
(and the unit / `Ring.inverse` corollaries in `GMC2DvdKFrameExtraction`).

## Gap 3 — formal `PowerSeries.logDeriv` + it commutes with ring homs

Repo source `GMC2DvdKFrameHSide` (death-star's `logDeriv` def + my `logDeriv_map`). No `PowerSeries.logDeriv`
in Mathlib. General over `[CommRing R]`, `ψ : R →+* S`.

```lean
namespace PowerSeries
variable {R : Type*} [CommRing R]
noncomputable def logDeriv (φ : R⟦X⟧) : R⟦X⟧ := derivativeFun φ * Ring.inverse φ

theorem derivativeFun_map {S} [CommRing S] (ψ : R →+* S) (f : R⟦X⟧) :
    map ψ (derivativeFun f) = derivativeFun (map ψ f) := by
  ext n; rw [coeff_map, coeff_derivativeFun, coeff_derivativeFun, coeff_map, map_mul, map_add,
    map_natCast, map_one]

theorem map_ringInverse_unit {S} [CommRing S] (ψ : R →+* S) {u : R⟦X⟧} (hu : IsUnit u) :
    map ψ (Ring.inverse u) = Ring.inverse (map ψ u) := by
  have hu' : IsUnit (map ψ u) := hu.map (map ψ)
  have h1 : map ψ u * map ψ (Ring.inverse u) = 1 := by rw [← map_mul, Ring.mul_inverse_cancel u hu, map_one]
  calc map ψ (Ring.inverse u)
      = (Ring.inverse (map ψ u) * map ψ u) * map ψ (Ring.inverse u) := by rw [Ring.inverse_mul_cancel _ hu', one_mul]
    _ = Ring.inverse (map ψ u) * 1 := by rw [mul_assoc, h1]
    _ = Ring.inverse (map ψ u) := mul_one _

theorem logDeriv_map {S} [CommRing S] (ψ : R →+* S) {u : R⟦X⟧} (hu : IsUnit u) :
    map ψ (logDeriv u) = logDeriv (map ψ u) := by
  rw [logDeriv, logDeriv, map_mul, derivativeFun_map, map_ringInverse_unit ψ hu]
end PowerSeries
```

## Status
All three verified to compile (`lake env lean`, exit 0) against the pinned Mathlib (v4.30.0). Gap 1 is the
cleanest first PR (Mathlib-only, one file). Gaps 2–3 follow. The full GMC(2) proof is a separate, larger
effort (see `gmc2-verification-and-mathlib-pr-readiness-kps-S128c154`).
