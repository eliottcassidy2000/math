# HYP-3956: the section formula + the x-side gap-surplus engine — the Lean-facing restatement of THM-599

**Status:** VERIFIED/DELIVERED — kind-pasteur-2026-07-01-S33. kps block (3950+).
**Canonized as:** the "Lean-facing restatement" section of THM-599 (01-canon).

## What was restated, and why it translates

THM-599's original proof rides compact-group Haar pushforward + Pontryagin duality — correct but
heavy for Lean. The restatement does the c-integral in CLOSED FORM per x-section:

    (SF)  A(U) = ∫₀¹ F_U(x) dx,   F_U(x) = Σ_gaps (g − w)⁺   (the gap surplus),

and every d-fold overlap becomes a one-circle integral of an explicit piecewise-linear integrand:
pair = ∫(w − ‖mx‖)⁺ = w² (one line), AP-triple = ∫(w − 2‖mx‖)⁺ = 2h². Mathlib inputs: Fubini on
T², AddCircle ball volume, x ↦ m·x measure-preserving. No duality, no subtorus Haar, no annihilator
lattices in the mechanized path.

## The enumeration improvement

F_U's breakpoints lie ONLY at { j/m } ∪ { (7j±1)/(7m) } over pairwise differences m — collisions
and gap-through-w crossings. Two consequences:
1. **Computation:** breakpoint counts 17–549 across the k = 2..13 argmins (the S32 c-engine's sets
   were Σ 4·u·u′ ≈ 10³–10⁵); per-pattern time drops to 0.00–0.04 s (~100× at k = 13).
2. **Lean form:** all breakpoints ∈ (1/M)ℤ, M = 7·lcm(diffs), so
   **(GT): A(U) = (1/M)·Σ_{r<M} F_U((2r+1)/(2M))** — a finite rational sum (`Finset.sum` of an
   explicit function). The ONE analytic lemma needed (generic, once): the integral of a
   piecewise-affine function is its midpoint sum on a refining grid. Everything else is rational
   arithmetic and sorting — exactly what Lean mechanizes well. (For evaluation inside Lean use the
   breakpoint sum, not the possibly-huge uniform grid; same lemma covers both.)

## Verification (lrc14_section_formula_xengine_kps.py + .out)

- X1: closed-form identities exact (36/49 at three pairs; 61/98 at two AP-triples).
- X2: **the x-engine reproduces the entire S32 c-engine ledger k = 2..13 in exact rational
  equality** — twelve rationals, two INDEPENDENT enumerations (c-side interval intersections vs
  x-side gap sums). This is the strongest correctness certificate the ledger can have short of Lean.
- X3: scale demo — the COMPLETE exact census of all 924 canonical 7-patterns of {1..13} in 5 s:
  min = 173/588 at (0,2,3,4,5,6,8) (= the sampled argmin's translate), runner-ups 521/1764,
  869/2940, 73/245; ALL ≥ witnessMP (min margin ×5.21). The k = 7 admissible pattern universe is
  now EXACTLY certified, not spot-checked.

## Lean plan (for whoever picks it up; pairs with mac-mini HYP-3856's normal form)

1. `gapSurplus (U : Finset ℤ) (x : ℚ) : ℚ` — sort, fold, (·−w)⁺: computable, no reals.
2. The piecewise-affine midpoint lemma (intervalIntegral, generic).
3. `A_eq_breakpointSum` : ties ∫F to the finite sum via the breakpoint enumeration
   (the completeness of the breakpoint set = "no collision ⟹ order constant ⟹ affine": a
   sorting-stability argument, mechanizable).
4. The ledger rows become `decide`/`norm_num` facts: A(pattern) = exact rational ≥ witnessMP.
Natural neighbors in-repo: LonelyRunnerMathlib.lean (defs), LRCThreeGapSampling (the count),
mac-mini's polygon/two_congruent_classes file.

## Artifacts
- 04-computation/lrc14_section_formula_xengine_kps.py (+ .out)
- THM-599 canon addendum (the restatement section)

## Depends on / relates to
THM-599, HYP-3953 (F), HYP-3955 (c-engine — now cross-verified), mac-mini HYP-3856, THM-565,
LonelyRunnerMathlib.lean, OPEN-Q-108.
