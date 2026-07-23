# GMC(2) formalization: post-closure assessment and what a Mathlib PR still needs

*klein-2026-07-22-S401. Owner: "finish the GMC(2) formalization and get it into the best possible
state for a mathlib PR"; then "pull and assess, no builds." An independent read-level assessment
after the fleet closed GMC(2) mid-session. I did not build (per directive); every claim below is
either read off the source or attributed to the agent who verified it.*

## What changed mid-session

I began this session mapping GMC(2) as **conditional** — reduced kernel-purely to the one-variable
Duistermaat–van der Kallen input (`DvdK1`), itself reduced to `SinglePolyCrux`, with the residual
being transpose bookkeeping. **That is now obsolete.** While I was working, the last hypothesis was
discharged:

- **kind-pasteur-S128c153** (`0784107bd`): `GMC2DvdKOmegaWiring.singlePolyCrux_holds` — the Ω-wiring
  over `Ω = AlgebraicClosure (LaurentSeries ℂ)` proves `SinglePolyCrux` outright; `gmc2_unconditional
  = gmc2_of_crux singlePolyCrux_holds`.
- **kind-pasteur-S128c154/b**: verification audit (green 8542-job chain, axioms clean, independent
  adversarial referee returned SOUND), the canonical capstone `GMC2Main.GMC2.gmc2`, 3 Mathlib
  PR-draft gap lemmas, and the minimal-chain footprint (69 essential / 19 off-chain of 88 modules).
- **mac-mini-S167**: built the same Ω-lift *independently and blind*, then **ceded** it as a
  duplicate; contributed the analytic-heart **line-audit** (no defect) and a full-chain axiom audit.
- **boxeph-S243**: bridge + divisibility load-bearing in the capstone.

**Correction on record:** my earlier framing in this session ("no file proves GMC(2)
unconditionally"; "one glue-composition away") was accurate when written and is now superseded. The
premise ledger is empty.

## My own read-level assessment (what I checked, not what I was told)

**The capstone statement is the genuine conjecture.** `GMC2Main.GMC2.gmc2 (P Q : MvPolynomial (Fin 2)
ℂ) (hnull : ∀ m ≥ 1, E (P^m) = 0) : ∃ N, ∀ m ≥ N, E (Q * P^m) = 0` — no residual hypothesis, and `E`
is the real Wick expectation (`E P = ∑ₛ P.coeff s · wt s`, `wt s = (s 0)!` when `s 0 = s 1`, else
`0`, i.e. `E[Zᵃ Z̄ᵇ] = a!·δ_{ab}`). This is the Mathieu–Zhao/moment form, not a weakened proxy.

**The premise ledger is empty.** `singlePolyCrux_holds : SinglePolyCrux` is a closed theorem taking
no hypothesis; `dvdK1_of_crux : SinglePolyCrux → DvdK1` then *removes the DvdK citation* rather than
assuming it. No `sorry` / `admit` / `native_decide` appears in code anywhere in the GMC2 tree (the
sole grep hit is a docstring sentence in `GMC2Main`).

**The Ω-wiring reads correctly.** Three points I specifically checked:
1. *The vanishing hypothesis is genuinely consumed* — `hvanish` flows into
   `smallRootFactor_coeff0_of_vanish` to produce `Pω.coeff 0 = algebraMap v`. The proof is not
   vacuous.
2. *The instance diamond is `rfl`-guarded* — the non-synthesizable `Algebra (RatFunc ℂ)
   (LaurentSeries ℂ)` is supplied locally via `rfToL.toAlgebra`, and the compatibility
   `(algebraMap (LaurentSeries ℂ) Ω).comp rfToL = algebraMap (RatFunc ℂ) Ω` is proved by `rfl`. This
   is exactly where a silent scalar-tower mismatch would hide, and it is closed definitionally.
3. *`Ω` is the algebraic closure, deliberately* — not `LaurentSeries ℂ` itself, because the small
   roots are ramified (`t^{1/M} ∉ ℂ((t))`). Both kps and mac-mini reached this independently, and
   agreed even on the constant (`cc = (−1)^{M+1} r₀`). For a claim this size, blind agreement on
   normalization is worth more than either derivation alone.
4. *`CharZero` is load-bearing and correctly scoped* — required by
   `smallRootFactor_coeff0_of_vanish` and false in char `p` (`t^p` has zero derivative); instantiated
   at `F = ℂ`. mac-mini confirmed `GMC2DvdKCharZeroClosing` carries an explicit `omit [CharZero F]`
   on the one step that does not need it.

**Verdict: it is looking good.** Structurally this is as strong as a formalization gets —
kernel-pure, sorry-free, green, statement-faithful, independently corroborated, adversarially
reviewed.

## Where the residual risk actually sits (the honest part)

Kernel-purity proves "no `sorry`"; it does **not** prove the statements are the ones the mathematics
needs. The shallowest link (the Ω-wiring) has now had the *most* eyes. The risk is concentrated
upstream, in decreasing order of scrutiny:

1. **The analytic heart** — `hderiv_final` / `smallRootFactor_coeff0_of_vanish` (`d_t h(0,t) = 0`
   under vanishing, in the `(F⸨x⸩)⟦t⟧` frame). This is the deepest link, and its line-audit was
   performed **by its own author** (mac-mini, S167). Self-audit of one's own deepest lemma is the
   weakest position in the review chain — mac-mini says so explicitly and recommends adversarial
   review target exactly here.
2. **`irreducible_Phi`, the Weierstrass factorization** (which leans on Mathlib's
   `PowerSeries.exists_isWeierstrassFactorization`), and **the height package** — mac-mini states
   plainly they did *not* re-verify these.

So the right next referee pass is upstream, not on the Ω-wiring. That is a recommendation, not a
defect claim: **I found no defect, and I did not audit those proofs line-by-line either.**

## What a Mathlib PR still needs (packaging, not mathematics)

- **Import hygiene — the biggest concrete blocker.** **32** GMC2 modules carry a blanket
  `import Mathlib` (including `GMC2DvdKOmegaWiring` itself, line 1). Mathlib requires precise
  imports. This is mechanical, and mac-mini flagged it **UNCLAIMED**.
- **Pruning the footprint** — 19 of 88 modules are off-chain (superseded additive route, early
  exploratory). Keep as research records; exclude from a PR footprint. Also **retire or clearly mark
  the superseded conditional surfaces** (`GMC2NC2Capstone`'s `..._of_heightWitnessSupplier`,
  `GMC2DvdKInterface`) so the "remaining inputs" story is not misread as still-open.
- **Exclude the two self-flagged defective files**: `GMC2DvdKUniqueChannelBypass` (kernel-checked but
  **vacuous** — inconsistent hypothesis, MISTAKE-240) and `GMC2HermiteNoCommonRoot` (chain not fully
  formal).
- **Ship the extracted gap lemmas first** — kps has already drafted and standalone-verified exactly
  the three I independently arrived at: the char-0 converse of `derivativeFun_C` (generalized to
  `[CommRing][NoZeroDivisors][CharZero]`), the geometric inverse `(1 − C w·X)⁻¹ = mk(wⁿ)`, and
  `logDeriv_map` (log-derivative commutes with ring homs). Gap 1 is the cleanest first PR. These are
  genuine Mathlib gaps and are independent of whether the GMC(2) result itself is ever submitted.
- **Generalize `ℂ` → a general char-0 field** where cheap (flagged diminishing-returns by kps).

## Fleet hygiene — a mistake of mine

mac-mini's S167 log records that "a sibling agent's broad `pkill -f "lake build"` killed several of
my runs (exit 144)". **That was me**, clearing what I thought were my own stuck builds; it killed
mac-mini's full-chain axiom audit mid-run and contributed to the `plausible` dep transiently failing
to resolve. In a shared working tree, `pkill` must be scoped to a specific PID, never to a pattern
another agent's jobs match. Recorded so the fleet does not repeat it.

*No new mathematics from me here. This is an independent read-level assessment: the proof is the
fleet's (kind-pasteur's Ω-wiring closing it, on death-star / mac-mini / boxeph / codex inputs).
My contribution is the corrected status, the confirmation that the premise ledger is empty and the
capstone statement is faithful, the localization of residual risk upstream of the Ω-wiring, and the
PR blocker list — of which import hygiene is the largest and is unclaimed.*
