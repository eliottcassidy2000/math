/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S129)
-/
import TournamentH7.LRCTorusReduction
import TournamentH7.LRCRankRigidity
import TournamentH7.LonelyRunner

/-!
# Route 2 assembly: `(A) ⟸ (C)` wired end-to-end, and the J-K top-level obligation (HYP-4652)

The LRC(14) "Route 2" (owner's intent, `00-navigation/LRC14-PROOF-MAP.md`) is

```
LRC14Statement
  ⟸ [J-K reduction]   CITATION (Jain–Kravitz / Giri–Kravitz 2024: the LRC-spectrum accumulation
                       points are governed by the relative spectra of rank-2 subtori)
  ⟸ (A)   no coupled rank-2 subtorus U has M(U) ∈ (1/13, 2/25)
  ⟸ (A)⟸(C)   projection floor (GREEN) + pigeonhole rigidity (core+wrapper GREEN) + the C-bridge
  ⟸ (C)   the 1-D 12-speed Farey gap: the dilated AP is the unique 12-family with M < 2/25.
```

This file **wires the middle** — `(A) ⟸ (C)` — into a single conditional theorem, composing the two
already-GREEN kernels

* the **projection floor** `TorusReduction.torus_loose_of_loose_direction` (opus-S99/S101): one
  loose integer direction `a·u + b·v` ⟹ the coupled 2-torus has a `2/25`-loose point; and
* the **pigeonhole rigidity wrapper** `RankRigidity.dep_of_infinite_common_proportional` (opus-S102):
  if every `(1,N)` direction is proportional to a *finite*-classifier vector, then `u, v` are
  dependent (rank ≤ 1) — the infinite-pigeonhole `Sym 12` argument, already formalized (the proof
  map's "wrapper OPEN" note is **stale**: it is GREEN).

The only remaining input is the **C-bridge** `hbridge`: *a `(1,N)` direction that is not loose is
(after centering) proportional to its finite dilated-AP ordering vector.*  This is exactly the crux
`(C)` in wrapper-ready form — it folds (i) `(C)` proper [not-loose ⟹ dilated AP] and (ii) the
centering normalization that makes the ordering classifier **finite** (`Sym 12`, not `Sym 12 × ℚ`:
centering removes the AP's additive offset, so the residual data is only the permutation).  It is
carried as a `Prop` hypothesis, not asserted — the honest boundary.

Given `hbridge`, the composition is unconditional: **a genuinely rank-2 integer torus is loose**
(`M(U) ≥ 2/25`, so `M(U) ∉ (1/13, 2/25)`), which is `(A)`.  The `by_cases` is on whether *some*
`(1,N)` direction is loose:
* **yes** → projection floor → torus loose;
* **no**  → `hbridge` makes *every* direction proportional → rigidity wrapper → `u,v` dependent →
  contradicts rank-2.

Finally we record the J-K top link as an obligation `Prop` and the trivial modus-ponens wiring to
`LRC14Statement`, so the fleet's `(C)` work registers against the top-level target (parallel to
Route 1's `LRC14.lrc14_from_witness_floor`).  The J-K dimension bookkeeping (13-speed LRC(14) →
rank-2 in `(ℝ/ℤ)¹²` → the 12-speed gap) is the citation to pin against the paper.
-/

namespace LonelyRunner
namespace Route2Assembly

/-- A coupled rank-2 integer torus with generators `u, v : Fin 12 → ℤ` is **loose** if some point
`(t,s)` keeps every speed `≥ 2/25` from `ℤ`: this is `M(U) ≥ 2/25`, hence `M(U) ∉ (1/13, 2/25)`. -/
def TorusLoose (u v : Fin 12 → ℤ) : Prop :=
  ∃ t s : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(u i : ℝ) * t + (v i : ℝ) * s - m|

/-- The 1-parameter integer direction `u + N·v` is **loose** at radius `2/25`: some `τ` keeps every
`(u i + N·v i)` at `≥ 2/25` from `ℤ`. -/
def DirLoose (u v : Fin 12 → ℤ) (N : ℤ) : Prop :=
  ∃ τ : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |((u i + N * v i : ℤ) : ℝ) * τ - m|

/-- The integer pair `(u,v)` spans a genuine **rank-2** torus: they are `ℝ`-linearly independent
(no nonzero real combination vanishes pointwise). -/
def Rank2 (u v : Fin 12 → ℤ) : Prop :=
  ∀ p q : ℝ, (∀ i, p * (u i : ℝ) + q * (v i : ℝ) = 0) → p = 0 ∧ q = 0

/-- **`(A) ⟸ (C)`, wired.**  A genuinely rank-2 integer torus is loose (`M(U) ≥ 2/25`), given the
C-bridge `hbridge` (a not-loose `(1,N)` direction is proportional to its finite dilated-AP classifier
vector).  Composes the GREEN projection floor and the GREEN rigidity wrapper; `hbridge` is the sole
`(C)`-obligation. -/
theorem torus_loose_of_rank2 (u v : Fin 12 → ℤ) (hindep : Rank2 u v)
    {S : Type*} [Finite S] (cls : ℤ → S) (W : S → Fin 12 → ℝ) (lamf : ℤ → ℝ)
    (hbridge : ∀ N : ℤ, ¬ DirLoose u v N →
        ∀ i, (u i : ℝ) + (N : ℝ) * (v i : ℝ) = lamf N * W (cls N) i) :
    TorusLoose u v := by
  by_cases hsome : ∃ N : ℤ, DirLoose u v N
  · -- some direction is loose → the projection floor makes the torus loose
    obtain ⟨N, τ, hτ⟩ := hsome
    apply TorusReduction.torus_loose_of_loose_direction u v
    refine ⟨1, N, τ, fun i m => ?_⟩
    have h1 : (1 * u i + N * v i : ℤ) = u i + N * v i := by ring
    rw [h1]; exact hτ i m
  · -- no direction is loose → hbridge makes every direction proportional → rigidity → dependent
    have hnone : ∀ N : ℤ, ¬ DirLoose u v N := fun N hN => hsome ⟨N, hN⟩
    have hprop : ∀ N : ℤ, ∀ i, (u i : ℝ) + (N : ℝ) * (v i : ℝ) = lamf N * W (cls N) i :=
      fun N => hbridge N (hnone N)
    obtain ⟨p, q, hpq, hzero⟩ :=
      RankRigidity.dep_of_infinite_common_proportional
        (fun i => (u i : ℝ)) (fun i => (v i : ℝ)) cls W lamf hprop
    exact absurd (Prod.ext (hindep p q hzero).1 (hindep p q hzero).2) hpq

/-- `(u,v)` is **affinely rank-2**: no nonzero real combination `p·u + q·v` is a constant vector.
Equivalently `{u, v, 𝟙}` are `ℝ`-independent — the natural non-degeneracy of a *coupled proper*
rank-2 subtorus (it is not contained in a coset of the diagonal `𝟙`-direction).  Strictly stronger
than `Rank2`. -/
def Rank2Affine (u v : Fin 12 → ℤ) : Prop :=
  ∀ p q : ℝ, (∀ i, p * (u i : ℝ) + q * (v i : ℝ) = p * (u 0 : ℝ) + q * (v 0 : ℝ)) → p = 0 ∧ q = 0

/-- **`(A) ⟸ (C)`, the honest form.**  Same conclusion as `torus_loose_of_rank2`, but the bridge
`hbridge` is now *literally* the crux `(C)`'s conclusion — "a not-loose `(1,N)` direction is a
dilated AP `a + d·(ordering)`" — with **no centering pre-baked**.  The affine rigidity wrapper
`dep_of_infinite_common_affine` handles the AP's additive offset internally (anchoring at
coordinate 0), and the offset-direction non-degeneracy is the explicit hypothesis `Rank2Affine`.
So the three concerns are cleanly separated: `(C)` = `hbridge`, centering = DONE (in the wrapper),
`𝟙`-exclusion = `Rank2Affine` (tied to J-K's "coupled proper" subtorus). -/
theorem torus_loose_of_rank2_affine (u v : Fin 12 → ℤ) (hindep : Rank2Affine u v)
    (hbridge : ∀ N : ℤ, ¬ DirLoose u v N →
        ∃ a d : ℤ, ∃ σ : Equiv.Perm (Fin 12), ∀ i, u i + N * v i = a + d * ((σ i).val : ℤ)) :
    TorusLoose u v := by
  by_cases hsome : ∃ N : ℤ, DirLoose u v N
  · obtain ⟨N, τ, hτ⟩ := hsome
    apply TorusReduction.torus_loose_of_loose_direction u v
    refine ⟨1, N, τ, fun i m => ?_⟩
    have h1 : (1 * u i + N * v i : ℤ) = u i + N * v i := by ring
    rw [h1]; exact hτ i m
  · have hall : ∀ N : ℤ, ¬ DirLoose u v N := fun N hN => hsome ⟨N, hN⟩
    choose a d σ hAP using fun N => hbridge N (hall N)
    obtain ⟨p, q, hpq, hconst⟩ :=
      RankRigidity.dep_of_infinite_common_affine
        (fun i => (u i : ℝ)) (fun i => (v i : ℝ))
        σ (fun τ i => ((τ i).val : ℝ)) (fun N => (d N : ℝ)) (fun N => (a N : ℝ))
        (by
          intro N i
          have h := hAP N i
          have h' : (u i : ℝ) + (N : ℝ) * (v i : ℝ)
              = (a N : ℝ) + (d N : ℝ) * ((σ N i).val : ℝ) := by exact_mod_cast h
          simp only []
          linear_combination h')
    have hconst' : ∀ i, p * (u i : ℝ) + q * (v i : ℝ) = p * (u 0 : ℝ) + q * (v 0 : ℝ) := hconst
    exact absurd (Prod.ext (hindep p q hconst').1 (hindep p q hconst').2) hpq

/-! ## The top-level obligations and the J-K wiring

We state `(A)` (no rank-2 torus in the gap — here in the stronger *all rank-2 tori loose* form that
the reduction actually delivers) and the C-bridge as global `Prop`s, then package the reduction and
the J-K citation. -/

/-- The top-level LRC(14) statement (mirror of `LRC14.LRC14Statement`), stated locally so this file
depends only on `Lonely`.  Every family of 13 nonzero integer speeds has a 14-lonely time. -/
def LRC14Target : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∃ t : ℝ, Lonely 14 v t

/-- **`(A)` (the rank-2 gap statement), delivered form.**  Every coupled rank-2 integer torus is
loose, so none has `M(U) ∈ (1/13, 2/25)`. -/
def AStatement : Prop := ∀ u v : Fin 12 → ℤ, Rank2 u v → TorusLoose u v

/-- **The global C-bridge obligation.**  For every rank-2 integer pair there is a finite classifier
witnessing that each not-loose `(1,N)` direction is proportional to its dilated-AP ordering vector.
This is `(C)` in wrapper-ready (centered) form — the sole open analytic/combinatorial node of
`(A)⟸(C)`. -/
def CBridge : Prop :=
  ∀ u v : Fin 12 → ℤ, Rank2 u v →
    ∃ (S : Type) (_ : Finite S) (cls : ℤ → S) (W : S → Fin 12 → ℝ) (lamf : ℤ → ℝ),
      ∀ N : ℤ, ¬ DirLoose u v N →
        ∀ i, (u i : ℝ) + (N : ℝ) * (v i : ℝ) = lamf N * W (cls N) i

/-- **`(A) ⟸ (C)` at the statement level.**  The global C-bridge implies `(A)`.  Pure composition of
`torus_loose_of_rank2` over all rank-2 pairs. -/
theorem AStatement_of_CBridge (h : CBridge) : AStatement := by
  intro u v hrank
  obtain ⟨S, hS, cls, W, lamf, hbridge⟩ := h u v hrank
  exact torus_loose_of_rank2 u v hrank cls W lamf hbridge

/-- **The J-K reduction obligation — ⚠ LIKELY FALSE, NOT a valid citation (MISTAKE-117, opus-S130).**
`AStatement → LRC14Target`.  This was believed to be the cited Jain–Kravitz / Giri–Kravitz 2024
"rank-2 governs the spectrum" reduction, but the S130 audit (verified against arXiv:2304.01462)
found that Giri–Kravitz study the **accumulation points** of the LR spectrum (`acc(S(n)) = S(n-1)`),
NOT the **supremum** that the LRC bounds — the abstract says verbatim *"Rather than attack this
conjecture, we study the structure of the sets S(n)."*  Controlling rank-2 subtori (accumulation-point
data) does NOT bound the sup, so this implication is UNWARRANTED and probably FALSE.  It is retained
only as the formal statement of the (broken) top link; **do not treat it as dischargeable**.  The
theorem below is a valid implication but rests on this likely-false hypothesis. -/
def JKReduction : Prop := AStatement → LRC14Target

/-- **Route 2, top to bottom (conditional).**  The J-K citation applied to `(A)`, with `(A)` supplied
by the C-bridge: `[J-K] + (C-bridge) ⟹ LRC(14)`.  Every hypothesis is an explicit obligation; the
interior (`(A)⟸(C)`, projection floor, rigidity) is GREEN.  Parallel to Route 1's
`lrc14_from_witness_floor`. -/
theorem lrc14_via_route2 (hjk : JKReduction) (hbridge : CBridge) : LRC14Target :=
  hjk (AStatement_of_CBridge hbridge)

#print axioms torus_loose_of_rank2
#print axioms torus_loose_of_rank2_affine
#print axioms AStatement_of_CBridge
#print axioms lrc14_via_route2

end Route2Assembly
end LonelyRunner
