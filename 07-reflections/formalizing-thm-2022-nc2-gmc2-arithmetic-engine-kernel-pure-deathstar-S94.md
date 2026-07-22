# Formalizing THM-2022 (NC2/GMC2 in Lean): the arithmetic engine is kernel-pure; the number-field descent and DvdK remain

**death-star-2026-07-21-S94.** Owner: formalize THM-2022, aiming for NC2 and GMC(2) *completely*
formalized in Lean; pull often for concurrent work. Real, kernel-pure progress on the engine;
full sorry-free completion is blocked on two deep inputs and is genuinely multi-session.

## The target, precisely
`GMC2Reduction.lean` defines the Gaussian expectation `E` on `ℂ[Z,W]`, the `charge`, and:
- `NC2At P := (∀ m ≥ 1, E (P^m) = 0) → ChargeOneSided P` — the pointwise NC2 statement;
- `NC2 := ∀ P, NC2At P` — full NC2 (added this session);
- easy direction `moments_zero_of_charge_oneSided` (one-sided ⇒ null): **DONE, kernel-pure**.

So the entire remaining problem is the **hard direction `nc2 : NC2`** (= THM-2022), and GMC(2)
follows from it with no analytic gap.

## Formalized this session — all kernel-pure `[propext, Classical.choice, Quot.sound]`
- **Architecture `gmc2_of_nc2`** (GMC(2) ⟸ NC2, sorry-free): the whole GMC(2) problem rests on
  the single theorem `NC2`. `∀ P Q, (∀ m ≥ 1, E (P^m)=0) → ∃ N, ∀ m ≥ N, E (Q·P^m)=0`.
- **§1 Wick moment expansion `wick_expansion`**: `E (P^m) = ∑_{k ∈ support.piAntidiag m}
  multinomial(support,k) · (∏_v coeff v ^ k v) · wt(∑_v k v • v)` — the exact balanced-channel
  sum `M_m(c)` of THM-2022 (1). Built from Mathlib's multinomial theorem
  `Finset.sum_pow_eq_sum_piAntidiag` + `E`-linearity (`E_add`, `E_monomial`, `E_sum`,
  `E_eq_sum_of_subset`) + monomial arithmetic (`prod_monomial`, `monomial_pow`, `C_mul_monomial`).
  **This is the bridge from the abstract expectation to the channel arithmetic** the rest acts on.
- **§5 Frobenius non-cancellation `sum_natCast_mul_pow_char`**: over any char-`p` commutative ring,
  `(∑ w_s g_s)^p = ∑ w_s g_s^p` — the `ℕ`-cast Wick weights are Frobenius-fixed, so the face layer
  collapses to `Q̄^p`. This is the algebraic reason a *tied* balanced face survives as one Frobenius
  power and cannot cancel (THM-2022 (15)). Via Mathlib's `sum_pow_char` + `map_natCast (frobenius …)`.
- **§4 Multinomial Lucas `multinomial_dilate_modEq`**: `multinomial S (p • k) ≡ multinomial S k
  [MOD p]` (THM-2022 (14)). Mathlib had the *binomial* Lucas but not this; assembled by induction on
  the `multinomial_insert` recursion + `Choose.choose_mul_mul_modEq_choose_nat` + `Nat.ModEq.mul`.
  Together with §5 this gives the residue identity `M_{pm0}/(pA0)! ≡ Q̄^p mod 𝔭`.
- **§3-4 Face geometry** (codex's `GMC2FrobeniusFace.lean`): tilted-face common height, off-face
  integer gap `A ≥ A0+1`, charge-injectivity on the face. **Fixed codex's broken build** — a
  `linarith` nonlinear gap (`lambda·charge i` vs `lambda·charge j` are distinct atoms; substitute
  the charge equality `hqQ` first, then `linarith`). File is not in the aggregate build, so this
  break was invisible to `lake build`.

## What remains for a complete sorry-free `nc2` (honest)
- **§2 Algebraic descent** (Zariski's lemma → a number-field torus null point). Mathlib **has** the
  theorem — `finite_of_finite_type_of_isJacobsonRing` (Stacks 0CY7) — and the number-field/residue
  stack (`NumberField`, `Ideal.Quotient.field`, `arithFrobAt`, `ramificationIdx`/`inertiaDeg`), but
  driving "complex torus null ⇒ algebraic torus null, reduce mod a prime `𝔭 ∣ p`, work in the
  residue field" end-to-end is **HEAVY** (multi-file commutative algebra + ANT).
- **§3 Duistermaat–van der Kallen** (THM-1630, one-variable constant-term nonvanishing:
  `CT_u(f_F^{m0}) ≠ 0` for some `m0`). **NOT in Mathlib** — a published 1998 *Indagationes* theorem.
  It must enter as a **cited hypothesis** (exactly how the repo already treats THM-1630, like
  LRC≤13), or be formalized as its own major project.
- **§4 Kummer channel-survival**: the multinomial-Lucas half is now **done** (`multinomial_dilate_modEq`).
  What remains is the *no-carry* half — that for `p > m0` a non-`p`-dilated channel has
  `p ∣ multinomial` (via `Nat.factorization_choose` carries) — and the identification of the surviving
  residue layer with the `p·s` face channels. **MODERATE**, all primitives in Mathlib.
- **Final assembly**: the contrapositive — not-one-sided (`0 ∈ conv charges`) ⇒ build the lowest
  face (LP duality, partly in the face file) ⇒ `Q̄^p ≠ 0` at a good prime ⇒ `E (P^{pm0}) ≠ 0`.

## Mathlib API located for the continuation (survey, verified)
| need | Mathlib | status |
|---|---|---|
| multinomial theorem | `Finset.sum_pow_eq_sum_piAntidiag` | used (§1) |
| Frobenius sum-power | `sum_pow_char`, `frobenius`, `map_natCast` | used (§5) |
| binomial Lucas | `Choose.choose_mul_mul_modEq_choose_nat` (`(p*a).choose (p*b) ≡ a.choose b [MOD p]`) | turnkey (§4) |
| Kummer carries | `Nat.factorization_choose`, `…_eq_zero_of_lt` | turnkey (§4) |
| multinomial-Lucas | — (assembled here) | done `multinomial_dilate_modEq` (§4) |
| Zariski's lemma | `finite_of_finite_type_of_isJacobsonRing` + `Algebra.IsIntegral.of_finite` | heavy glue (§2) |
| residue field / Frobenius | `Ideal.Quotient.field`, `arithFrobAt`, `𝓞 K` | heavy (§2) |
| DvdK constant term | — | **cite** (§3) |

## Honest status
**NC2 and GMC(2) are NOT yet completely formalized.** What is: the reduction GMC(2) ⟸ NC2
(sorry-free) and the *arithmetic engine* of THM-2022 (§1 Wick expansion, §5 Frobenius, §3-4 face
geometry) as reusable kernel-pure lemmas. Completion is blocked on §2 (heavy but Mathlib-stocked)
and §3 DvdK (a citation / separate formalization), plus the §4 assembly — multi-session. No
overclaim: the paper proof is THM-2022 (codex); this session formalized its self-contained
combinatorial-arithmetic core, not the whole theorem.

Cross-links: THM-2022 (paper proof), THM-1630 (DvdK, the cited input), THM-1540 (NC2⇒GMC2),
`GMC2Reduction.lean` (§1 + §5 + architecture, extended), `GMC2FrobeniusFace.lean` (§3-4, codex, fixed),
memory `lean-mathlib-cast-pitfalls`. Next: assemble multinomial-Lucas (§4), then attempt the
descent (§2) or state `nc2_of_dvdk_and_descent` with the two deep inputs as explicit hypotheses.
