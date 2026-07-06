/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S108)
-/
import Mathlib

/-!
# The coverer dichotomy: the sum-product seam as one lemma (HYP-4406)

The LRC crux is a SUM-PRODUCT rigidity (S107): at the prime 13, `{1,…,12}` is the additive
interval AND the multiplicative group `(ℤ/13)*`, and gap-emptiness is the rigidity of that
coincidence.  Its cleanest concrete face — where the MULTIPLICATIVE (covering) structure
forces the ADDITIVE (height) structure — is this CRT dichotomy:

  a runner `v ≡ r (mod 13)` with `r ∈ {1,…,12}` that COVERS the modulus `q = r` (`r ∣ v`)
  is EITHER `v = r` (the unlifted base runner) OR `v ≥ 14r`.

Because `v = r·m` with `m ≡ 1 (mod 13)` (from `r·m ≡ r` and `gcd(r,13)=1`), so `m = 1` or
`m ≥ 14`.  For the UNIQUE-coverer moduli `r ∈ {7,…,12}` (the only multiple of `q = r` in
`{1..12}` is `r` itself), this says: a pinned covering family either keeps `r` unlifted, or
its `q = r` coverer sits at height `≥ 14r ≥ 98`.  With `LRCDivisorProtection` (no multiple of
`r` ⟹ loose) this is exactly the S78 pigeonhole / S80 height-forcing seam in one line:
covering ⟹ unlifted-or-high.
-/

namespace LonelyRunner
namespace CovererDichotomy

/-- **The coverer dichotomy.**  A positive `v ≡ r (mod 13)` with `1 ≤ r ≤ 12` and `r ∣ v` is
either the base runner `v = r` or height-forced `14·r ≤ v`.  (CRT: `v = r·m`,
`m ≡ 1 (mod 13)`, so `m = 1` or `m ≥ 14`.) -/
theorem coverer_height (r v : ℤ) (hr1 : 1 ≤ r) (hr12 : r ≤ 12) (hv : 0 < v)
    (hres : v % 13 = r % 13) (hdvd : r ∣ v) : v = r ∨ 14 * r ≤ v := by
  obtain ⟨m, rfl⟩ := hdvd
  have hr0 : 0 < r := by omega
  -- m > 0 since r·m > 0 and r > 0
  have hm0 : 0 < m := by
    by_contra h
    push_neg at h
    nlinarith [hv, mul_nonneg hr0.le (neg_nonneg.mpr h)]
  -- 13 ∣ r·(m − 1) from r·m ≡ r (mod 13)
  have h13 : (13 : ℤ) ∣ r * (m - 1) := by
    have he : r * (m - 1) = r * m - r := by ring
    rw [he]; omega
  -- 13 is prime and 13 ∤ r (since 1 ≤ r ≤ 12), so 13 ∣ (m − 1)
  have hdm : (13 : ℤ) ∣ (m - 1) := by
    have hp : Prime (13 : ℤ) := by norm_num
    rcases hp.dvd_mul.mp h13 with hr' | hm'
    · exact absurd (Int.le_of_dvd hr0 hr') (by omega)
    · exact hm'
  -- m ≡ 1 (mod 13), m > 0  ⟹  m = 1 (base) or m ≥ 14 (height)
  obtain ⟨j, hj⟩ := hdm
  have hjnn : 0 ≤ j := by omega
  rcases eq_or_lt_of_le hjnn with hj0 | hj0
  · left
    have hm1 : m = 1 := by omega
    rw [hm1]; ring
  · right
    nlinarith [mul_nonneg hr0.le (show (0 : ℤ) ≤ j - 1 by omega), hj]

#print axioms coverer_height

end CovererDichotomy
end LonelyRunner
