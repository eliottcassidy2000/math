# LRC(14) Lean corpus — Mathlib-readiness report

**mac-mini-2026-07-06-S11 (HYP-4362).**  An honest assessment of what the
TournamentH7 corpus needs to be Mathlib-submittable (per the owner's directive:
get everything built and submittable, don't split into PRs).

## Executive summary

- **Sorry-free: YES.**  Zero `sorry` in code across all 213 modules (every
  `sorry` occurrence is in a docstring documenting "no sorry").  Machine-checked
  by `AxiomAudit.lean` / `LRC14AxiomAudit.lean` (`#print axioms`, no `sorryAx`).
- **The main blocker: `native_decide`.**  204 `native_decide` tactic uses across
  ~40 modules add the `Lean.ofReduceBool` axiom (trusts the compiler, not the
  kernel).  **Mathlib does not accept `native_decide` in merged code.**  These
  must become kernel `decide`, `Decidable`-instance proofs, or genuine proofs —
  and many are used *because* `decide` is too slow (heavy finite census gates).
- **A kernel-pure conditional CORE exists.**  The canonical conditional theorem
  and its entire dependency chain are `native_decide`-free (verified: 0 in
  LRCHcompSurface, LRCHarmonicGate, LRCOneWindowPeel, LRCTightAPFreeRider,
  LRCBandMargin, LRCMergeExclusion, LRCClusterGcdSharp, LRCPeelCompression).
  The 204 `native_decide` uses are confined to the ALTERNATE data-census route
  (LRCWindowData22 — a 31,750-line generated `winData22` census — and
  exploratory modules), which the conditional core does not need.

## The Mathlib-submittable core (kernel-pure, conditional)

    lrc14_of_dichotomy_and_corner (cite : LRCUpTo13)
        (hdich : TightLooseDichotomy) (hcorner : CornerLonely) : LRC14Statement
    lrc14_of_spread_dichotomy_and_corner  -- the same, dichotomy on spread bases only

- `cite : LRCUpTo13` — the owner-sanctioned LRC(≤13) citation, a NAMED HYPOTHESIS
  (LRC13Citation.lean), not an axiom, not a sorry.  For a real Mathlib PR this
  would carry the standard preprint caveat or the ≤13 results as cited lemmas.
- `hdich`, `hcorner` — the two remaining REAL-ANALYTIC hypotheses.
- **Axiom footprint CONFIRMED (machine-checked):** both
  `lrc14_of_dichotomy_and_corner` and `lrc14_of_spread_dichotomy_and_corner`
  depend on `[propext, Classical.choice, Quot.sound]` ONLY — kernel-pure, NO
  `Lean.ofReduceBool`, NO `sorryAx`.  The conditional core is Mathlib-axiom-clean.

**CORRECTION integrated (kps-S11 / MISTAKE-110):** the bounded-modulus finite
census (`TemplateDichotomy`, Q50 at s ≤ 50) is FALSE — free-modulus witnesses
(s = 53, and pinned-only bound → ∞) exist, so NO fixed-modulus template closes
LRC(14).  `lrc14_of_template_and_corner` is kernel-pure and true-as-an-implication
but its hypothesis is unprovable — a DEAD reduction, correctly marked in
LRCTemplateSurface.lean.  **The remaining hypotheses are irreducibly
real-analytic** (`TightLooseDichotomy` = ∃ a real t*), NOT a finite check.  This
also tempers my S8–S10 "the tail is one finite census" framing: the residue
bridge finitizes the *pinned* part, but the free-modulus part is real-analytic.

## Genuinely Mathlib-general, standalone-ready

- **LRCTorusRate.lean** (mac-mini-S10, kernel-pure): the Lipschitz-net max
  transfer — `exists_net_ge` (an L-Lipschitz f's max transfers to any ε-dense
  set with loss ≤ L·ε), `lipschitz_ge_at_near`, `two_point_rate`.  This is
  general metric-space content, LRC-independent, and is the closest thing in the
  corpus to a directly-submittable Mathlib lemma (modulo checking it isn't
  already derivable from `LipschitzWith` + density API).
- The integer-arithmetic atoms (LRCKStratification `dInt_scale`/`binder_dvd`,
  the `rational_point_margin` atom) are clean, kernel-pure, `decide`/`omega`/
  `linarith`-proved — Mathlib-style, though LRC-specific.

## What a full Mathlib submission needs (effort estimate)

1. **Eliminate `native_decide` from any submitted theorem** (204 uses).  Two
   sub-cases: (a) small finite checks → replace `native_decide` with `decide`
   (works where the kernel can reduce in time); (b) heavy census gates (winData22,
   factorial atoms) → these will NOT reduce in the kernel; they need either
   genuine proofs or exclusion from the submission.  **Largest single task.**
2. **Pick ONE canonical statement** — `lrc14_of_dichotomy_and_corner` — and prune
   the superseded surfaces (covering-far/22, template) from the submission set.
3. **The `LRCUpTo13` citation** — package as cited results with the standard caveat.
4. **Style pass** — Mathlib naming (`lowerCamelCase` defs, `snake_case` theorems
   per section conventions), docstrings on every public decl, `variable`
   discipline, no `open` leakage.  Mechanical but broad (213 modules → the
   submitted subset).
5. **The remaining hypotheses are OPEN** (real-analytic dichotomy + corner), so
   the submission is a CONDITIONAL result: "LRC(14) reduces to
   {LRC(≤13)} + {TightLooseDichotomy} + {CornerLonely}", machine-checked.  This
   is a legitimate, publishable conditional theorem — the honest deliverable.

## This session's concrete build actions (results)

- Full-corpus `lake build`: the kernel-pure CORE (≈150 modules incl. the whole
  dichotomy chain) built green; the slow census/exploratory tail (winData22 etc.)
  compiles but is the native_decide route.
- **Axiom audit of the canonical core: CONFIRMED kernel-pure** — `[propext,
  Classical.choice, Quot.sound]` only, no `Lean.ofReduceBool`, no `sorryAx`.
- My modules (LRCTorusRate, LRCKStratification, LRCLiftFloorRows/Assembly,
  LRCResiduePinning, LRCWindowMargin13, LRCKernelGate13, LRCLiftRigidityRows)
  verified native_decide-free and kernel-pure.

## Bottom line

The **conditional theorem is Mathlib-axiom-clean today**: LRC(14) ⟸ {LRC(≤13)}
+ {TightLooseDichotomy} + {CornerLonely}, sorry-free, `[propext, Classical.choice,
Quot.sound]` only.  What stands between this and an actual Mathlib PR is (1)
pruning to this one surface + its native_decide-free chain (dropping the
winData22/exploratory route), and (2) a mechanical style/naming/docstring pass on
that subset.  The two remaining hypotheses are genuinely open (real-analytic), so
the submission is a legitimate CONDITIONAL result — not the full theorem.
