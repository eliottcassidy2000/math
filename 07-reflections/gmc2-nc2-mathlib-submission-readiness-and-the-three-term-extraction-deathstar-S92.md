# GMC(2)/NC2 formalization: Mathlib-submission readiness, and the general three-term no-common-root extraction

**death-star-2026-07-21-S92** (HYP-8805). Owner: get the NC2/GMC(2) results formalized to Mathlib-submission
quality. Worked the Lean corpus (`04-computation/lean/TournamentH7`), verified kernel-purity, and extracted the
one genuinely Mathlib-new, fully-general result into a standalone PR-ready file.

## The submittable results (all verified kernel-pure)
`#print axioms` on each returns **only** `[propext, Classical.choice, Quot.sound]` — the Mathlib-standard set;
**no `sorryAx`, no `Lean.ofReduceBool`** (native_decide). So they clear Mathlib's hardest gate.

| result | file | Mathlib-readiness |
|---|---|---|
| **`ThreeTerm.no_common_root` / `exists_nonvanishing`** — a monic three-term recurrence with `b n ≠ 0` has no two consecutive members sharing a root, over any integral domain | **`ThreeTermRecurrence.lean` (NEW, this session)** | **PR-ready.** General, `autoImplicit false`, docstringed, kernel-pure; **Mathlib has no such lemma** (checked). |
| `GMC2.mathieuZhao_of_charge_pos` — if every monomial of `P ∈ ℂ[Z,W]` has charge ≥ 1 then `E(Q·P^m)=0` for `m ≫ 0` (the NC2 ⇒ GMC(2) charge-arithmetic step) | `GMC2Reduction.lean` | Kernel-pure, clean `MvPolynomial (Fin 2) ℂ` proof; the *statement* is Mathlib-expressible (no Mathieu-Zhao name needed). PR-viable after light packaging. |
| `GMC2Hermite.*` (Hermite/trExp instances) | `GMC2HermiteNoCommonRoot.lean` | Kernel-pure; superseded for Mathlib by the general `ThreeTermRecurrence` (Mathlib already has `Polynomial.hermite`). |

## What I did this session
**Extracted and generalized the flagship result.** The `GMC2Hermite.ThreeTerm.no_common_root` proof (over `ℝ`,
via `linarith`) → new file `TournamentH7/ThreeTermRecurrence.lean`, **generalized to any `[CommRing R] [IsDomain R]`**
(replacing `linarith` by `linear_combination hrec`, `[Zero R]` on the structure, `mul_eq_zero` for the domain
step), with `set_option autoImplicit false`, full docstrings, and an `ℝ`-Hermite instance. **Builds clean (40s),
kernel-pure.** This is the abstract core of Favard/orthogonal-polynomial no-common-roots, proved from the
recurrence alone — a natural, self-contained Mathlib contribution.

## Honest scope — what is and isn't submittable
- **Post-integration correction (THM-2022):** NC2 and GMC(2) are proved in
  the repository by Frobenius amplification of the lowest balanced Wick face,
  but that proof is not yet formalized in Lean. The distinction here is
  proved-on-paper versus formalized/submission-ready, not proved versus open.
- **NC2 and GMC(2) as *full Lean theorems* are not yet available.** The
  current submission payload is the proved reductions and recurrence lemmas,
  not yet a formal proof of THM-2022. That is exactly what should go to Mathlib.
- **The Lean-proved DvdEZ/NC2 ⇒ GMC(2) implication** (`mathieuZhao_of_charge_pos` and the fleet's
  `GMC2Reduction`) is real and kernel-pure — a genuine, citable formalized implication.
- Project-internal computational lemmas (`GMC2MomentBasics`: `E_Pspan_sq`, `E_Pfake_*` by kernel `decide`) are
  correct and kernel-pure but **specific instances**, not general theorems — they stay in-repo, not for Mathlib.

## PR-packaging progress (S93 update)
1. **`ThreeTermRecurrence` recast — DONE.** Added the actual polynomial sequence `ThreeTerm.poly : ℕ → R[X]`
   (`noncomputable`, since `Polynomial` over a general `CommRing` is), the evaluation bridge
   `eval_poly : (T.poly n).eval x = T.p n x`, and the Mathlib-idiomatic **`no_common_root_poly`** stated in
   `Polynomial.IsRoot` form. Kept the function-level `p`/`no_common_root` as the computational core. All 5 public
   results verified kernel-pure (`[propext, Classical.choice, Quot.sound]`).
2. **Minimal imports — DONE.** Stripped `import Mathlib` to the exact four modules `#min_imports` reports and a
   test-build confirms: `Mathlib.Algebra.Polynomial.Eval.Defs`, `Mathlib.Data.Real.Basic`,
   `Mathlib.Tactic.LinearCombination`, `Mathlib.Tactic.Positivity`. Build dropped **8475 → 1202 jobs** — a real,
   large reduction, and the standard Mathlib-PR requirement (`shake`-clean) met by construction.
3. **Hermite bridge — DONE, and it is itself Mathlib-new.** New file `TournamentH7/HermiteThreeTerm.lean`
   (kernel-pure, warning-free) proves what Mathlib lacked:
   - **`Polynomial.derivative_hermite_succ`**: the ladder relation `(hermite (n+1))' = (n+1) • hermite n`.
     Single-step induction from Mathlib's `hermite_succ` — the `H_n'` terms cancel, closed by the `module` tactic.
   - **`Polynomial.hermite_recurrence`**: the classical **three-term recurrence**
     `hermite (n+2) = X·hermite (n+1) − (n+1)·hermite n` (a one-liner from the ladder relation).
   - **`ThreeTermRecurrence.hermite_no_common_root`**: *no two consecutive Hermite polynomials share a root* —
     `Polynomial.hermite` exhibited as the `ThreeTerm ℤ` instance `hermiteZ` (`a≡0`, `b m = m`), so the general
     `no_common_root_poly` applies verbatim. This is the compelling **application** that motivates the abstract
     lemma for Mathlib: a concrete, named orthogonal family that Mathlib already ships, now known root-coprime.
   To make `hermiteZ` land cleanly the structure hypothesis was **weakened to `hb : ∀ n, b (n+1) ≠ 0`** (only
   `b 1, b 2, …` occur in the recurrence; Hermite has `b 0 = 0`) — a strict generalization, rebuilt kernel-pure.
   (`module`, `derivative_smul`, and `smul_eq_C_mul` all need care: the Module-ℤ smul in the statements does not
   syntactically match the `SMulZeroClass` smul those lemmas are stated with, so they fire via `map_smul` /
   `show … from`, not bare `rw`.)
4. **`mathieuZhao_of_charge_pos`**: still on `import Mathlib`; minimal-import trim pending (statement stands alone).

## The owner's "divide by (pA₀)!" is exactly THM-2022's crux (vindication + correction)
THM-2022 (codex) proves NC2/GMC(2) in full, and its **§4 normalization is literally the owner's S91 directive**:
*divide the moment of order `p·m0` by the common factorial `(p·A0)!`*. That instinct was right and is the hinge of
the actual proof. What my S91 execution got **wrong** (now MISTAKE-215): the residue after dividing is **not** a
Vandermonde of channel degrees — it is the Frobenius power **`Q̄^p`** of the lowest-balanced-face constant term
`Q = CT_u(f_F^{m0})`, isolated by Kummer's carry formula (`p > m0` ⟹ no-carry ⟹ only `p`-dilated channels
survive) and Lucas' congruence (`binom(pm0; ps) ≡ binom(m0; s)`). Likewise the S89/S90 "equal Vandermonde nodes =
regular/Paley scores" identification is withdrawn (MISTAKE-214); the free-probability central-trinomial identity
survives only as analogy. **Honest net:** the owner's normalization idea is vindicated as the mechanism of the
proof; the tournament/Vandermonde/Paley *dressing* I put on it was the error. The correct object at the wall is
`Q̄^p`, not a discriminant.

## Status
`ThreeTermRecurrence.lean` is **Mathlib-PR-ready**: general (`[CommRing R] [IsDomain R]`), kernel-pure, docstringed,
`autoImplicit false`, minimal imports, both function- and `Polynomial`/`.IsRoot`-level statements, new to Mathlib
(checked: Mathlib lacks even the three-term no-common-root). Wired into the project root and building. The GMC(2)
charge-arithmetic reduction (`mathieuZhao_of_charge_pos`) is kernel-pure and PR-viable after a minimal-import trim.
The full NC2/GMC(2) proof is **THM-2022** (proved on paper; Frobenius/Kummer/Lucas over a lowest balanced face,
with the one-variable Duistermaat–van der Kallen constant-term theorem THM-1630 as its deep analytic input — a
citation, not in Mathlib); its formalization is the real multi-session target and is correctly **excluded** from
the present submission claim.
A **second** PR-ready file now accompanies it: `HermiteThreeTerm.lean` supplies the three-term Hermite recurrence
(Mathlib-missing) and the "consecutive Hermite share no root" corollary, kernel-pure — the concrete Mathlib-object
application of the abstract lemma.
Cross-links: GMC2Reduction/GMC2HermiteNoCommonRoot/GMC2MomentBasics (existing kernel-pure corpus), S62 (Hermite
no-common-root origin), THM-2022 (full proof; corrects S89–S91 via MISTAKE-214/215), THM-1630 (DvdK constant-term,
the citation input), THM-1540 (NC2⇒GMC(2)), memory `lean-mathlib-cast-pitfalls`. Files
`04-computation/lean/TournamentH7/TournamentH7/ThreeTermRecurrence.lean` (recast + minimal imports, weakened `hb`,
S93) and `…/HermiteThreeTerm.lean` (three-term Hermite recurrence + no-common-root, S93). HYP-8805.
