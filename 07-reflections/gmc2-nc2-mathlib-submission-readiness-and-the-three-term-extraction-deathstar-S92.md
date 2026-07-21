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
- **NC2 and GMC(2) as *full theorems* are OPEN** (the fleet's live residual is radial-channel noncancellation at
  the regular/Paley wall — my S87–S91). So "our NC2/GMC(2) results" = the **proved reductions and lemmas**, not a
  proof of the conjecture. That is exactly what should go to Mathlib.
- **The Lean-proved DvdEZ/NC2 ⇒ GMC(2) implication** (`mathieuZhao_of_charge_pos` and the fleet's
  `GMC2Reduction`) is real and kernel-pure — a genuine, citable formalized implication.
- Project-internal computational lemmas (`GMC2MomentBasics`: `E_Pspan_sq`, `E_Pfake_*` by kernel `decide`) are
  correct and kernel-pure but **specific instances**, not general theorems — they stay in-repo, not for Mathlib.

## Remaining work toward an actual PR (scoped, not blocking)
1. **`ThreeTermRecurrence`** (closest to ready): optionally recast `p : ℕ → R → R` as `Polynomial R` with
   `.IsRoot`, and connect `hermiteReal` to Mathlib's `Polynomial.hermite` (show they satisfy the same recurrence).
   Decide placement (`Mathlib/RingTheory/Polynomial/…` near orthogonal polynomials). The mathematics is done and
   verified; this is packaging/idiom.
2. **`mathieuZhao_of_charge_pos`**: strip the project `import Mathlib` to minimal imports; the statement stands
   alone (a charge-graded expectation-vanishing lemma) and needs no Mathieu-Zhao scaffolding.
3. A Mathlib PR wants module docstrings (present), no `set_option` hacks (none), and CI-green — the extracted file
   already satisfies these.

## Status
`ThreeTermRecurrence.lean` is **Mathlib-PR-ready** (general, kernel-pure, docstringed, `autoImplicit false`, new
to Mathlib), wired into the project root and building. The GMC(2) charge-arithmetic reduction is kernel-pure and
PR-viable with light packaging. NC2/GMC(2) themselves remain open (correctly excluded from any submission claim).
Cross-links: GMC2Reduction/GMC2HermiteNoCommonRoot/GMC2MomentBasics (existing kernel-pure corpus), S62 (Hermite
no-common-root origin), S87–S91 (the open NC2 residual), memory `lean-mathlib-cast-pitfalls`. New file
`04-computation/lean/TournamentH7/TournamentH7/ThreeTermRecurrence.lean`. HYP-8805.
