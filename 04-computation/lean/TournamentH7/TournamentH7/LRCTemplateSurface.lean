/-
  ⚠️ CORRECTION (kind-pasteur-2026-07-05-S11, HYP-4137 / MISTAKE-110):
  `TemplateDichotomy` BELOW IS FALSE as stated (the `s ≤ 50` bound).  An explicit,
  independently-verified single-scale Fin-13 family (height ~1e22) satisfies every
  hypothesis, is not tight-shaped, and has NO 2/25-witness at any `s ≤ 50` (min
  witness s = 53).  Reason: the profile pins residues only mod `q ≤ 25`; witnesses
  at FREE moduli (27,32,49, primes 29..47,53,…) depend on residues the profile does
  not control and are killed by a CRT lift.  Only witnesses at PINNED-ONLY moduli
  (`q | lcm(2..25)`) are height-independent.  The theorem `lrc14_of_template_and_corner`
  is kernel-pure and TRUE as an implication, but `htempl` is UNPROVABLE — a DEAD
  reduction.  CORRECTED target: witness at a pinned-only `q ≤ Q₀` (Q₀ ≥ 69, census max;
  height-independent).  The real surface `lrc14_of_dichotomy_and_corner` (loose = ∃ real
  tstar) is UNAFFECTED — the counterexample has a real witness at t = 13/53.  See
  05-knowledge/results/lrc14_q50_refutation_kps_S11.md.

  TournamentH7.LRCTemplateSurface — THE Q50 CRUX, PINNED INTO THE LEAN SURFACE
  (kind-pasteur-2026-07-05-S10, HYP-4127).

  mac-mini-S55 crystallized the spectral-gap crux: the profile filters are
  CRT-ray-periodic (their floating 7-cluster passes every filter at every
  frozen scale), so NO residue filter can bound heights — the absolute-height
  mechanism must be WITNESS TEMPLATES: every one of their 511,947 full-profile
  census survivors (height ≤ 48) carries a `2/25`-margin witness at a REDUCED
  denominator `q ≤ 45`.  THE Q50 CONJECTURE: this is uniform.

  This file makes that crux FORMAL:

  * `TemplateWitness s k` — the fixed integer conditions
    `μ(s) ≤ (vᵢ·k) mod s ≤ s − μ(s)` with `μ(s) = ⌈2s/25⌉ = (2s+24)/25`,
    DECIDABLE per family (a `Decidable` instance + a `decide` demonstrator on
    the spectrum's second value).
  * `loose_of_template` — a template yields the dichotomy's loose branch, via
    the HYP-4102 rational-point atom.
  * `TemplateDichotomy` — THE REFINED OPEN PREDICATE: every primitive compressed
    covering family at the argmax peel is tight-shaped OR satisfies a template
    with modulus `≤ 50`.  The loose side is now a FIXED FINITE template family.
  * `lrc14_of_template_and_corner` — THE NEW PINNED SURFACE:

        LRC(14)  ≤  citation  +  TemplateDichotomy  +  CornerLonely.

  Proving `TemplateDichotomy` = the Q50 conjecture + tight classification; its
  loose half is a finite CRT statement (finitely many `(s ≤ 50, k)` templates,
  each a residue condition), exactly the shape the census/pyramid machinery
  sweeps.  Kernel-pure; no `sorry`; `decide` only in the demonstrator example.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace TemplateSurface

open LRC14 HcompSurface HarmonicGate

/-- The integer margin floor `⌈2s/25⌉`, in integer arithmetic. -/
def muInt (s : ℤ) : ℤ := (2 * s + 24) / 25

/-- `μ(s)/s ≥ 2/25` — the template margin clears the gap threshold. -/
lemma muInt_ge (s : ℤ) (hs : 0 < s) : (2 : ℝ)/25 ≤ (muInt s : ℝ) / s := by
  have h25 : 2 * s ≤ 25 * muInt s := by
    unfold muInt
    omega
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  rw [le_div_iff₀ hsR]
  have hR : (2 : ℝ) * s ≤ 25 * (muInt s : ℝ) := by exact_mod_cast h25
  linarith

/-- **The template**: fixed integer residue conditions on the base making every
base runner `2/25`-far at `t = k/s`.  Pure integer arithmetic — decidable. -/
def TemplateWitness (s k : ℤ) (v : Fin 13 → ℤ) (istar : Fin 13) : Prop :=
  ∀ i, i ≠ istar → muInt s ≤ (v i * k) % s ∧ (v i * k) % s ≤ s - muInt s

instance (s k : ℤ) (v : Fin 13 → ℤ) (istar : Fin 13) :
    Decidable (TemplateWitness s k v istar) := by
  unfold TemplateWitness
  infer_instance

/-- The spectrum's second value carries the `(s, k) = (25, 2)` template — the
decidable pipeline end-to-end, by `decide`. -/
example : TemplateWitness 25 2 ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 24, 999] 12 := by
  decide

/-- **Template ⟹ the dichotomy's loose branch** (via the HYP-4102 atom). -/
theorem loose_of_template (v : Fin 13 → ℤ) (istar : Fin 13) (s k : ℤ) (hs : 0 < s)
    (ht : TemplateWitness s k v istar) :
    ∃ tstar : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, 2/25 ≤ |(v i : ℝ) * tstar - m| :=
  loose_branch_of_cert v istar (2/25) k s (muInt s) hs (muInt_ge s hs) ht

/-- **THE TEMPLATE DICHOTOMY** — the Q50-refined open predicate: every primitive
compressed covering family at the argmax peel is tight-shaped, or satisfies a
witness template with modulus `≤ 50`.  Empirical basis: 511,947/511,947 census
survivors to height 48 witness at `q ≤ 45` (mac-mini S54/S55). -/
def TemplateDichotomy : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → tupleGcd v = 1 →
    ∀ istar : Fin 13, (∀ i, |v i| ≤ |v istar|) →
      (∃ c : ℤ, 2 ≤ c ∧ ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j) ∨
      (∃ s k : ℤ, 0 < s ∧ s ≤ 50 ∧ TemplateWitness s k v istar)

/-- The template dichotomy refines the surface's dichotomy. -/
theorem tightLooseDichotomy_of_template (h : TemplateDichotomy) :
    TightLooseDichotomy := by
  intro v hv hcov hcomp hprim istar hmax
  rcases h v hv hcov hcomp hprim istar hmax with htight | ⟨s, k, hs, _, ht⟩
  · exact Or.inl htight
  · exact Or.inr (loose_of_template v istar s k hs ht)

/-- **THE NEW PINNED SURFACE**: LRC(14) from the citation node, the template
dichotomy (finite loose side), and the corner. -/
theorem lrc14_of_template_and_corner (cite : LonelyRunner.LRCUpTo13)
    (htempl : TemplateDichotomy) (hcorner : CornerLonely) :
    LRC14Statement :=
  lrc14_of_dichotomy_and_corner cite (tightLooseDichotomy_of_template htempl) hcorner

#print axioms loose_of_template
#print axioms tightLooseDichotomy_of_template
#print axioms lrc14_of_template_and_corner

end TemplateSurface
end LonelyRunner
