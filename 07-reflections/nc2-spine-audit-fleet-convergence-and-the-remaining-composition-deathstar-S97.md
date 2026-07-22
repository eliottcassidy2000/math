# NC2 formalization: the fleet converged — codex built the number-field machinery; audit says it's sound; one composition remains

**death-star-2026-07-21-S97.** Owner: work the number-field/valuation machinery for NC2, pull often,
look for creative bypasses. Pulled — and found codex had **concurrently built the entire machinery**
(commit `f11c59ae2`, a 33-module `GMC2*` spine under `TournamentH7/`, aggregated by
`GMC2Formalization`). My independent S94–S96 arithmetic (§1 Wick, §4 Lucas/no-carry/off-face, §5
Frobenius/face-sum) is **subsumed** by codex's equivalents. So this session's value is an **audit**
of the fleet capstone plus honest coordination, not duplicate machinery.

## What codex built (the machinery the owner asked me to work)
- **`GMC2FrobeniusResidue`** — independently re-derived every arithmetic lemma I had: `weighted_sum_pow_char`
  (= my `sum_natCast_mul_pow_char`), `prime_dvd_multinomial_of_sum_eq_mul_of_not_dvd` (= my no-carry
  `multinomial_dvd_of_exists_not_dvd`), `multinomial_dilate_modEq`, `prime_dvd_normalized_factorial_of_gap`
  (= my `factorial_dilate_dvd`), `face_sum_frobenius`, `face_sum_ne_zero`. Full fleet convergence.
- **`GMC2ResidueAssembly.three_case_sum_eq_frobenius` / `_ne_zero`** — the *abstract char-`p` bypass* I was
  about to build: over any `[CommSemiring R] [ExpChar R p]`, if non-dilated and off-face channels vanish and
  face channels are `w·x^p`, the channel sum is `(∑ w·x)^p`; nonzero if the face sum is. **No number fields.**
- **`GMC2IntegralFaceSeedDescent.exists_finite_field_moment_point_preserving_integral_lowest_face_seed`** —
  **THE descent**, and codex's *creative bypass*: a **direct finite-field route** through a finite-type
  `ℤ`-algebra instead of a number field. Under `DvdK1` + null + not-one-sided it produces a finite residue
  field (prime char `p`, `IsField`, `Finite`), a torus point `w` (all units, nonzero), preservation of every
  integral zero relation, vanishing of every positive moment relation, and a **nonzero lifted face seed**.
  "Its residue characteristic may be learned before selecting the characteristic-dependent normalized
  moment" — the bypass that dodges choosing a prime up front.
- **`GMC2DvdKInterface.DvdK1`** — DvdK as an explicit premise (not an axiom), exactly the roadmap
  recommendation (like citing LRC≤13).

## Audit (my contribution): the capstone is SOUND
- **`DvdK1` is stated correctly.** `∀` injective charges `q:ι→ℤ`, nonzero coeffs `c:ι→ℂ`, if the charges
  straddle zero (`ChargesStraddleZero`, i.e. two-sided), then `∃ m≥1, aeval c (constantTermRelation q m) ≠ 0`.
  This is *precisely* DvdK Theorem 2 + Remark 3 in exact-support form — not vacuous, not too strong, and taken
  as a hypothesis, not an axiom. Correct.
- **The descent is well-formed.** It threads `GMC2FaceSeed.exists_nonzero_lowest_face_seed hDvdK` through the
  integral specialization; the conclusion's preservation clause (`aeval P.coeff f = 0 → aeval w f = 0`) is the
  right "integral zero relations are preserved," and the moment-vanishing + seed-nonzero clauses are consistent.
- **Spine is sorry-free** (grep: only `GMC2HermiteNoCommonRoot` — the unrelated old Hermite file — and a
  docstring occurrence of the word "sorry" in `GMC2Reduction`). Kernel-purity of `three_case_sum_ne_zero` and
  the descent: build in progress.

## The one remaining task, and why I did not race it
Per `GMC2Formalization`'s own note: *"the remaining composition task is to instantiate the abstract three-case
assembly with the normalized Wick channels and the exact dilation image, then feed its contradiction through
`GMC2IntegralFaceSeedDescent` to obtain the conditional theorem `DvdK1 → NC2`."* The shape is clear:
`aeval w (normalizedMomentRelationInt …)` is both `= 0` (null → integral relation, preserved by the descent) and
`≠ 0` (`three_case_sum_ne_zero`, face seed nonzero) — contradiction. The pieces are all present; the composition
is deep in codex's 33-module API (channel/dilate/face predicates, the `piAntidiag` dilation image, the reference
channel). It is **codex's active capstone**; racing it from outside their API would collide and duplicate. The
clean glue it still needs — turning codex's `ℕ`-divisibility lemmas into the `term = 0` hypotheses of
`three_case_sum` — is a two-line char-`p` cast (`CharP.cast_eq_zero_iff`), which I verified in scratch
(`multinomial_cast_eq_zero`, `factorial_ratio_cast_eq_zero`) and offer to codex.

## Honest status for the owner
NC2/GMC(2) is **very close**: the entire proof is formalized and sorry-free *except* the final composition into
`DvdK1 → NC2`, plus DvdK itself (a cited hypothesis, ~person-months to formalize per S95). The number-field/
valuation machinery the owner asked me to build **exists** — codex built it (finite-field bypass), converging
with my arithmetic. My net contribution this session: verifying it is correct, and the coordination.
Cross-links: `GMC2Formalization` (codex spine), S94–S96 reflections (my converged engine), S95 DvdK roadmap,
memory `nc2-gmc2-lean-formalization-state`. HYP-8805.
