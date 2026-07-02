/-
Copyright (c) 2026 the math-research swarm. All rights reserved.
Released under Apache 2.0 license.
Authors: klein-2026-07-01-S89 (HYP-3845)

# Regular-polygon partitions and the Davenport–Mirsky–Newman–Rado argument

If the vertices of a regular `N`-gon are colored so that every color class is itself the
vertex set of a regular polygon, then two color classes are congruent (they have the same
number of vertices).

Encoding: the vertices are `0, 1, ..., N-1` (as `Finset.range N`). A regular `m`-gon
inscribed in the `N`-gon with vertices among the `N`-gon's vertices is exactly a residue
class `{k < N : k ≡ a [MOD q]}` where `q = N / m` divides `N`. So the statement is the
**Davenport–Mirsky–Newman–Rado theorem on `Z/N`**: a partition of `{0,...,N-1}` into at
least two residue classes with pairwise distinct moduli is impossible.

The proof is the classical generating-function argument: evaluate the partition identity
`∑_i ∑_{k ∈ class i} ω^k = ∑_{k<N} ω^k = 0` at a primitive root of unity `ω` of order the
**largest** modulus `q i₀`. Every class with a strictly smaller modulus sums to `0`
(a full geometric cycle), while the largest-modulus class contributes
`ω^{a} · (N / q i₀) ≠ 0`. Contradiction.

This is the `Z/N` shadow of the "continuous Mirsky–Newman" rigidity (THM-594(C)) used in
the Lonely Runner tower floor: exact tilings by distinct-speed danger arcs are forbidden
for the same reason distinct-modulus exact covers are.

LRC-repo cross-references: HYP-3845, THM-594(C), opus F1 (HYP-3902), mac-mini HYP-3850(a).
-/
import Mathlib

open Finset

namespace PolygonPartition

/-- The residue class `{k < N : k % q = a % q}`, i.e. the vertex set of the inscribed
regular `(N/q)`-gon through vertex `a` with step `q`. -/
def residueClass (N q a : ℕ) : Finset ℕ :=
  (range N).filter (fun k => k % q = a % q)

/-- A residue class with modulus `q ∣ N` is the image of `range (N / q)` under
`j ↦ (a % q) + j * q`. -/
lemma residueClass_eq_image {N q a : ℕ} (hq : 0 < q) (hdvd : q ∣ N) :
    residueClass N q a = (range (N / q)).image (fun j => a % q + j * q) := by
  have hr : a % q < q := Nat.mod_lt _ hq
  ext k
  simp only [residueClass, mem_filter, mem_range, mem_image]
  constructor
  · rintro ⟨hkN, hkr⟩
    refine ⟨k / q, Nat.div_lt_div_of_lt_of_dvd hdvd hkN, ?_⟩
    rw [← hkr]
    exact Nat.mod_add_div' k q
  · rintro ⟨j, hj, rfl⟩
    have hmul : (j + 1) * q ≤ N := by
      calc (j + 1) * q ≤ (N / q) * q := Nat.mul_le_mul_right q hj
        _ = N := Nat.div_mul_cancel hdvd
    have hexp : j * q + q ≤ N := by rw [add_mul, one_mul] at hmul; omega
    refine ⟨by omega, ?_⟩
    rw [Nat.add_mul_mod_self_right]
    exact Nat.mod_mod_of_dvd a dvd_rfl

/-- Geometric evaluation of a residue-class exponential sum at `ω`:
`∑_{k ∈ class} ω^k = ω^(a % q) * ∑_{j < N/q} (ω^q)^j`. -/
lemma residueClass_expSum {N q a : ℕ} (hq : 0 < q) (hdvd : q ∣ N) (ω : ℂ) :
    ∑ k ∈ residueClass N q a, ω ^ k
      = ω ^ (a % q) * ∑ j ∈ range (N / q), (ω ^ q) ^ j := by
  rw [residueClass_eq_image hq hdvd,
    sum_image fun x _ y _ h => Nat.eq_of_mul_eq_mul_right hq (by omega), mul_sum]
  refine sum_congr rfl fun j _ => ?_
  rw [pow_add, pow_mul']

/-- **The polygon-partition theorem (Davenport–Mirsky–Newman–Rado on `Z/N`).**
If residue classes with moduli `q i ∣ N` partition `{0, ..., N-1}` (every vertex lies in
exactly one class) and there are at least two classes, then two classes share a modulus:
the coloring contains two congruent inscribed regular polygons. -/
theorem exists_eq_modulus
    {N : ℕ} (hN : 0 < N) {ι : Type*} [Fintype ι] [DecidableEq ι]
    (a q : ι → ℕ) (hqpos : ∀ i, 0 < q i) (hdvd : ∀ i, q i ∣ N)
    (hcard : 2 ≤ Fintype.card ι)
    (hpart : ∀ k ∈ range N, ∃! i : ι, k % q i = a i % q i) :
    ∃ i j : ι, i ≠ j ∧ q i = q j := by
  by_contra hcon
  have hinj : Function.Injective q := by
    intro i j hij
    by_contra hne
    exact hcon ⟨i, j, hne, hij⟩
  -- the largest modulus
  have hne : (Finset.univ : Finset ι).Nonempty := Finset.univ_nonempty_iff.mpr
    (Fintype.card_pos_iff.mp (by omega))
  obtain ⟨i₀, -, hmax⟩ := Finset.exists_max_image Finset.univ q hne
  -- all other moduli are strictly smaller
  have hlt : ∀ i : ι, i ≠ i₀ → q i < q i₀ := fun i hi =>
    lt_of_le_of_ne (hmax i (mem_univ i)) (fun h => hi (hinj h))
  -- the largest modulus is at least 2
  have h2 : 2 ≤ q i₀ := by
    by_contra hq2
    obtain ⟨i₁, i₂, h12⟩ :=
      Fintype.exists_pair_of_one_lt_card (by omega : 1 < Fintype.card ι)
    have e1 : q i₁ = 1 := by
      have := hmax i₁ (mem_univ i₁); have := hqpos i₁; omega
    have e2 : q i₂ = 1 := by
      have := hmax i₂ (mem_univ i₂); have := hqpos i₂; omega
    exact h12 (hinj (e1.trans e2.symm))
  -- the primitive root of unity of order q i₀
  set ω : ℂ := Complex.exp (2 * Real.pi * Complex.I / (q i₀ : ℂ)) with hω
  have hprim : IsPrimitiveRoot ω (q i₀) :=
    Complex.isPrimitiveRoot_exp _ (by omega)
  have hω1 : ω ≠ 1 := hprim.ne_one (by omega)
  have hωN : ω ^ N = 1 := by
    obtain ⟨M, hM⟩ := hdvd i₀
    rw [hM, pow_mul, hprim.pow_eq_one, one_pow]
  -- the full vertex sum vanishes
  have htotal : ∑ k ∈ range N, ω ^ k = 0 := by
    rw [geom_sum_eq hω1, hωN, sub_self, zero_div]
  -- the classes partition range N
  set A : ι → Finset ℕ := fun i => residueClass N (q i) (a i) with hA
  have hcover : range N = Finset.univ.biUnion A := by
    ext k
    simp only [mem_biUnion, mem_univ, true_and, hA, residueClass, mem_filter]
    constructor
    · intro hk
      obtain ⟨i, hi, -⟩ := hpart k hk
      exact ⟨i, hk, hi⟩
    · rintro ⟨i, hk, -⟩
      exact hk
  have hdisj : Set.PairwiseDisjoint (↑(Finset.univ : Finset ι)) A := by
    intro i _ j _ hij
    refine Finset.disjoint_left.mpr fun k hki hkj => hij ?_
    simp only [hA, residueClass, mem_filter] at hki hkj
    obtain ⟨w, -, huniq⟩ := hpart k hki.1
    exact (huniq i hki.2).trans (huniq j hkj.2).symm
  -- split the vanishing sum over the classes
  have hsplit : ∑ i : ι, ∑ k ∈ A i, ω ^ k = 0 := by
    rw [← Finset.sum_biUnion hdisj, ← hcover, htotal]
  -- classes with smaller modulus vanish
  have hzero : ∀ i : ι, i ≠ i₀ → ∑ k ∈ A i, ω ^ k = 0 := by
    intro i hi
    rw [hA]
    rw [residueClass_expSum (hqpos i) (hdvd i)]
    have hnd : ¬ (q i₀ ∣ q i) := by
      intro hd
      have h1 := Nat.le_of_dvd (hqpos i) hd
      have h2 := hlt i hi
      omega
    have hq1 : ω ^ q i ≠ 1 := fun h => hnd ((hprim.pow_eq_one_iff_dvd _).mp h)
    rw [geom_sum_eq hq1]
    have hcyc : (ω ^ q i) ^ (N / q i) = 1 := by
      rw [← pow_mul, Nat.mul_div_cancel' (hdvd i), hωN]
    rw [hcyc, sub_self, zero_div, mul_zero]
  -- the largest-modulus class does not vanish
  have hbig : ∑ k ∈ A i₀, ω ^ k = ω ^ (a i₀ % q i₀) * ((N / q i₀ : ℕ) : ℂ) := by
    rw [hA]
    rw [residueClass_expSum (hqpos i₀) (hdvd i₀), hprim.pow_eq_one]
    simp
  have hMpos : 0 < N / q i₀ :=
    Nat.div_pos (Nat.le_of_dvd hN (hdvd i₀)) (hqpos i₀)
  have hbig_ne : ∑ k ∈ A i₀, ω ^ k ≠ 0 := by
    rw [hbig]
    exact mul_ne_zero (pow_ne_zero _ (Complex.exp_ne_zero _))
      (Nat.cast_ne_zero.mpr (by omega))
  -- contradiction
  have hsingle : ∑ i : ι, ∑ k ∈ A i, ω ^ k = ∑ k ∈ A i₀, ω ^ k :=
    Fintype.sum_eq_single i₀ hzero
  rw [hsingle] at hsplit
  exact hbig_ne hsplit

end PolygonPartition
