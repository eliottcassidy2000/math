/-
  TournamentH7.LRCCoveringStrata — DISCHARGING STRATA OF `CoveringComplete`
  (kind-pasteur-2026-07-06-S51, HYP-4677).

  opus-S129 (`LRCCoveringDisjunction`) reduced the crux `(C)` to one obligation:
  `CoveringComplete` = *every non-AP 12-family has a `HasCoveringWitness`* (a modulus
  `q` with a clearing rotation `(c,μ)`, `μ/q ≥ 2/25`, every speed off the hole mod `q`).
  This file discharges that obligation **stratum by stratum** by *producing* the
  witness — the "existence half" opus's disjunction consumes.

  Stratum here: the **small-modulus / Fan-Sun gcd layer** (klein-S147).  A family that
  misses a multiple of some `q ≤ 12` (no speed `≡ 0 mod q`) has the covering witness
  `(q, c=1, μ=1)`: every residue lies in `[1, q−1]`, and `2q ≤ 25` clears `2/25`.  This
  is the layer that disposes of every family whose lift *breaks a divisibility*
  (klein-S147's 100%-verified branch); it feeds straight into `HasCoveringWitness`.

  Kernel-pure.
-/
import Mathlib
import TournamentH7.LRCCoveringDisjunction

namespace LonelyRunner
namespace CoveringStrata

open CoveringDisjunction

/-- **Small-modulus stratum ⟹ covering witness.**  If no speed is divisible by some
`q ≤ 12` (equivalently every `v i % q ∈ [1, q−1]`), then `v` has a `HasCoveringWitness`:
the witness `(q, c=1, μ=1)` clears at `t = 1/q` (`μ/q = 1/q ≥ 2/25` since `2q ≤ 25`). -/
theorem smallmod_hasWitness (v : Fin 12 → ℤ) (q : ℤ) (hq : 0 < q) (hq12 : 2 * q ≤ 25)
    (hnd : ∀ i, 1 ≤ v i % q ∧ v i % q ≤ q - 1) :
    HasCoveringWitness v := by
  refine ⟨q, 1, 1, hq, by norm_num, by omega, by linarith, fun i => ?_⟩
  have := hnd i
  simpa using this

/-- **Divisibility form.**  The same, phrased as "no speed is a multiple of `q`"
(`q ∤ v i`), for `0 < q ≤ 12` and speeds taken with nonnegative residues. -/
theorem no_multiple_hasWitness (v : Fin 12 → ℤ) (q : ℤ) (hq : 0 < q) (hq12 : 2 * q ≤ 25)
    (hnd : ∀ i, ¬ (q ∣ v i)) (hpos : ∀ i, 0 ≤ v i % q) :
    HasCoveringWitness v := by
  refine smallmod_hasWitness v q hq hq12 (fun i => ⟨?_, ?_⟩)
  · rcases (hpos i).lt_or_eq with h | h
    · omega
    · exact absurd (Int.dvd_of_emod_eq_zero h.symm) (hnd i)
  · have := Int.emod_lt_of_pos (v i) hq
    omega

#print axioms smallmod_hasWitness
#print axioms no_multiple_hasWitness

end CoveringStrata
end LonelyRunner
