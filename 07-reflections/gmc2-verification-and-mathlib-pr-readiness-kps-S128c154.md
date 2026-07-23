# GMC(2): verification audit + Mathlib-PR-readiness verdict

*kind-pasteur-2026-07-22-S128c154. Owner: double-check the GMC(2) formalization is in the best state
possible and ready for a Mathlib PR. This is the audit result.*

## Correctness — as strong as a formalization gets

| Check | Result |
|---|---|
| `gmc2_unconditional` builds green | **yes** — 8542-job dependency chain (`GMC2DvdKOmegaWiring`) |
| `#print axioms gmc2_unconditional` | `[propext, Classical.choice, Quot.sound]` — no `sorryAx`, no `native_decide` |
| `sorry`/`admit`/`native_decide`/`axiom` in the chain | **none** (grep hits were docstrings + `#print axioms` lines) |
| Independent adversarial referee (full chain) | **SOUND** — statement genuine, non-vacuous, non-circular, no gap |
| GMC2 chain vs the failing LRC modules | **fully independent** — no GMC2 file imports any LRC module |

**The statement is the real GMC(2)** (referee-confirmed): `∀ P Q : MvPolynomial (Fin 2) ℂ, (∀ m ≥ 1,
E(Pᵐ) = 0) → ∃ N, ∀ m ≥ N, E(Q·Pᵐ) = 0`, with `E` the genuine Wick/Gaussian expectation
(`∑ P.coeff s · wt s`, `wt s = if s 0 = s 1 then (s 0)! else 0`). `SinglePolyCrux` is a genuine ∀→∃
(empty packet gives product `1 ≠ algebraMap(cc·X)`), and the vanishing hypothesis is genuinely consumed
by `hderiv_final`. The Ω-instance (`Ω = AlgebraicClosure (LaurentSeries ℂ)`, base algebra via
`rfToL.toAlgebra`) is sound: the two `rfl` guards (`IsScalarTower.of_algebraMap_eq`, `halg`) force the
intended faithful composite, and a single `hS_of_dvd_value` application makes an instance mismatch a
type error, not a silent bug.

## The one real caveat about the *repo* (not GMC(2))

`lake build TournamentH7` (the whole root aggregator) **fails**, but exclusively on in-progress **LRC**
modules (`LRCCoherentBlockerChronology`, `LRCPairTowerValuation`, `LRCTwoCircleII`, …) — a different
research area. **No GMC2 module imports any LRC module**, so this does not touch GMC(2)'s correctness or
a GMC(2) PR. For a fully-green root, the LRC WIP needs fixing (fleet's lane).

## Mathlib-PR readiness — honest scoping

**The full proof is not a single Mathlib PR.** It is a large research formalization: 87 GMC2 files,
47 in the `gmc2_unconditional` chain, 32 `import Mathlib` wholesale, TournamentH7-specific naming, and
many exploratory/dead modules. Upstreaming it whole would be a multi-week refactor (specific imports,
prune dead modules, Mathlib naming, consolidate, and the analytic heart line-audited for readability —
the kernel already verifies it, but Mathlib review wants human-legible proofs).

**What IS near-term PR-able** — reusable general lemmas, confirmed genuine Mathlib gaps:
1. **`eq_C_of_derivativeFun_eq_zero`** (char-0 converse of `derivativeFun_C`) — the cleanest. Extracted,
   **generalized to `[CommRing R] [NoZeroDivisors R] [CharZero R]`**, and **verified to compile**
   standalone (draft in scratch `PR-draft-charzero-derivative.lean`). Ready as the first PR.
2. **`(1 − C w · X)⁻¹ = ∑ wⁿ Xⁿ`** (`oneSubCX_mul_mkGeom` / `inverse_oneSubCX`) — general `section`,
   extractable; no `PowerSeries` geometric-inverse currently in Mathlib.
3. **`derivativeFun_map`** (formal derivative commutes with `map ψ`) and a formal
   **`PowerSeries.logDeriv` + `logDeriv_map`** — clean and general; no `PowerSeries.logDeriv` in Mathlib.

## Verdict

GMC(2) is **rigorously verified** — kernel-checked end-to-end *and* independently adversarially reviewed
SOUND, with clean axioms and no sorries, and isolated from the repo's unrelated LRC build failures. It is
in the best state a formalization can be for correctness. For Mathlib, the path is to upstream the 3
general-lemma gaps first (the char-0 converse is drafted and compiles); the full GMC(2) proof is a
separate, larger effort.

## Cross-links
`GMC2DvdKOmegaWiring` (gmc2_unconditional) · reflection `gmc2-proven-unconditional-omega-wiring-closed`
(the proof) · `GMC2DvdKCharZeroClosing` (char-0, Mathlib-only) · `GMC2DvdKFrameExtraction` (geometric
inverse) · `GMC2DvdKFrameHSide` (logDeriv_map) · HYP-9020.
