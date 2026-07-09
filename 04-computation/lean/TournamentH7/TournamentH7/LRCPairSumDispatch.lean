/-
  TournamentH7.LRCPairSumDispatch — the grid-free THM-668 pair-sum dispatch leg
  (kind-pasteur-2026-07-09-S114).

  mac-mini's THM-668 (pair-sum ruler theorem): every local max of `m(t) = minᵢ ‖vᵢ t‖` sits at a
  rational `t = p/q` whose denominator is a **pair-sum modulus** `q = vᵢ + vⱼ`.  So the LRC witness is
  always such a rational, and loneliness can be certified by **integer residue arithmetic on `q`** — no
  grid.  This file supplies the Lean primitive mac-mini asked me to wire in: a witness `t = p/q` clears
  every runner iff each residue `(vₗ·p) mod q` lands in the band `[q/14, 13q/14]`, and then
  `Mreach ≥ 1/14` (via `Mreach_ge_of_lonely_instant`, kps-S99b).

  The **Schur-triple count `E₃`** is the *kill budget* of this dispatch (mac-mini THM-668 + opus LEM-015):
  each height-3 relation `vᵢ+vⱼ ∣ vₗ` forces `‖vₗ·t‖ = 0` at the `(i,j)` ruler, deleting a live ruler.
  The AP (`E₃ = C(13,2) = 78`, the maximiser) kills 12 rulers and keeps exactly one live band-ruler
  `q = 14`, sitting at tangency `σ = 0` — the extremal.  So `E₃` is the scalar the dispatcher reads to
  size the live-ruler budget.  Self-contained apart from `LRCReachWitness`.
-/
import Mathlib
import TournamentH7.LRCReachWitness

namespace LonelyRunner
namespace LRC14Concrete

/-- **Residue-band ⟹ cleared (the dispatch primitive).**  For integers `a` and `q > 0`, if the residue
`a mod q` lies in the band `[q/14, 13q/14]` (written as the integer inequalities `q ≤ 14·(a % q)` and
`14·(a % q) ≤ 13·q`), then the rational `a/q` is at least `1/14` from every integer:
`nearInt(a/q) ≥ 1/14`.  This is the integer-exact clearance test at a fixed modulus — no grid. -/
theorem nearInt_intDiv_ge_of_band (a q : ℤ) (hq : 0 < q)
    (hlo : q ≤ 14 * (a % q)) (hhi : 14 * (a % q) ≤ 13 * q) :
    (1 : ℝ) / 14 ≤ nearInt ((a : ℝ) / (q : ℝ)) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  set r : ℤ := a % q with hr
  have hr0 : (0 : ℝ) ≤ (r : ℝ) := by exact_mod_cast Int.emod_nonneg a hq.ne'
  have hrq : (r : ℝ) < (q : ℝ) := by exact_mod_cast Int.emod_lt_of_pos a hq
  -- decompose a/q = ⌊⌋ + r/q, hence fract(a/q) = r/q
  have key : (a : ℝ) = ((a / q : ℤ) : ℝ) * (q : ℝ) + (r : ℝ) := by
    have h : q * (a / q) + a % q = a := Int.ediv_add_emod a q
    have hc := congrArg (fun z : ℤ => (z : ℝ)) h
    push_cast at hc
    rw [hr]; linarith
  have hdiv : (a : ℝ) / (q : ℝ) = ((a / q : ℤ) : ℝ) + (r : ℝ) / (q : ℝ) := by
    rw [key]; field_simp
  have hfract : Int.fract ((a : ℝ) / (q : ℝ)) = (r : ℝ) / (q : ℝ) := by
    rw [hdiv, Int.fract_intCast_add]
    exact Int.fract_eq_self.mpr ⟨by positivity, by rw [div_lt_one hqR]; exact hrq⟩
  unfold nearInt
  rw [hfract, le_min_iff]
  refine ⟨?_, ?_⟩
  · rw [le_div_iff₀ hqR]
    have : (q : ℝ) ≤ 14 * (r : ℝ) := by exact_mod_cast hlo
    linarith
  · have h14 : 14 * (r : ℝ) ≤ 13 * (q : ℝ) := by exact_mod_cast hhi
    have hle : (r : ℝ) / (q : ℝ) ≤ 13 / 14 := by rw [div_le_iff₀ hqR]; linarith
    linarith

/-- **The pair-sum dispatch leg (THM-668, grid-free).**  Given a pair-sum modulus `q > 0` and a
multiplier `p`, if every runner's residue `(vₗ·p) mod q` lies in the band `[q/14, 13q/14]`, then the
rational witness `t = p/q` makes the runner set lonely: `Mreach v ≥ 1/14`.  This is the exact
integer-arithmetic dispatch primitive — enumerate pair-sum events (provably complete by THM-668), test
the band, done.  Composes directly with `Mreach_ge_of_lonely_instant`. -/
theorem mreach_ge_of_pairsum_band (v : Fin 13 → ℤ) (p q : ℤ) (hq : 0 < q)
    (hband : ∀ i, q ≤ 14 * ((v i * p) % q) ∧ 14 * ((v i * p) % q) ≤ 13 * q) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  apply Mreach_ge_of_lonely_instant
  refine ⟨(p : ℝ) / (q : ℝ), fun i => ?_⟩
  have hi := hband i
  have hcl := nearInt_intDiv_ge_of_band (v i * p) q hq hi.1 hi.2
  rw [show ((v i : ℝ) * ((p : ℝ) / (q : ℝ))) = (((v i * p : ℤ)) : ℝ) / (q : ℝ) by push_cast; ring]
  exact hcl

end LRC14Concrete
end LonelyRunner
