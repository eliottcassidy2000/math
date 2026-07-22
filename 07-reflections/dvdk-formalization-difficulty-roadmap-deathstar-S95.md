# Duistermaat–van der Kallen (one variable): Lean-formalization difficulty roadmap

**death-star-2026-07-21-S95.** Owner asked, alongside the NC2/GMC2 continuation, to roadmap how hard
it would be to formalize DvdK — the one input to THM-2022 that is a *citation*, not yet formal.
Researched against the primary sources (DvdK's `powers.pdf`, Monsky arXiv:0906.1836, van den
Essen–Schoone arXiv:2305.10062 §5, ESV arXiv:0908.2609) and the local Mathlib checkout (v4.30.0).

## What THM-2022 actually needs
Exactly the **one-variable** DvdK "Theorem 2 + Remark 3" (THM-1630 in canon): for the face Laurent
polynomial `f_F(u) = Σ_{i∈F} c_i u^{q_i}` with `0 ∈ conv{q_i}` (two-sided), some power has a nonzero
constant term — `∃ m0 ≥ 1, CT(f_F^{m0}) ≠ 0`.  Char 0 is **essential** (over 𝔽_p, `u+u⁻¹` has
`CT(fⁿ) ≡ 0 ∀n`), which is why THM-2022 applies DvdK over a number field, *before* reducing mod `p`.

## Two proof routes, and the verdict
**(A) DvdK's original analytic proof** (residues + Liouville, ~1 page): `F(t)=Σ CT(fⁿ)t^{n-1}` is a
contour integral; a residue `−1/t` at `u=0` (two-sidedness) cannot be cancelled by the branch
residues `O(t^{−1−1/m})`; if `F` were single-valued around every nonzero critical value it would be
entire and vanish at ∞, so Liouville forces `F≡0`, contradiction.
- **Mathlib status: multi-person-YEAR.** Liouville, Cauchy integrals, circle integrals, removable
  singularities, the meromorphic library are all present — but the **residue theorem / argument
  principle is ABSENT**, and, fatally, the proof needs **local Puiseux expansions of the algebraic
  function `ζ(r)` and monodromy** (non-single-valuedness) around branch points. Multivalued-function /
  monodromy theory over ℂ is a large unbuilt area of Mathlib. **Not recommended.**

**(B) Monsky's purely algebraic proof** (cleanest exposition: van den Essen–Schoone §5). Reformulates
as: `CT(fⁿ)=0` for all large `n` ⟹ `f` one-sided. Steps: form `U(X)=Xˢ(1−z·f) ∈ ℂ((z))[X]`; split
its roots (Newton–Puiseux over `ℂ((z))`) by a sign-valuation into `S⁺`/`S⁻` (inside/outside, none of
valuation 0); char-0 separability gives simple roots and a partial-fraction identity for
`1/(1−zf)`; expand each `(X−aᵢ)⁻¹` as a **two-sided** formal Laurent series, apply the band-limited
functional `L=CT`, get a closed form `W(z)`; show `W(z)≠1`, contradicting `W(z)=L(1)=1`.
- **Mathlib status: a few MONTHS (≈4–9 person-months for an expert).** Present and reusable:
  `LaurentPolynomial`, `LaurentSeries` with X-adic `Valued` instance, `AlgebraicClosure`,
  `Polynomial.PartialFractions` (existence+uniqueness), `FieldTheory.Separable`, `Polynomial.roots`
  /`Splits`, `RatFunc.laurent`, Krasner's lemma + `SpectralNorm` + `Henselian`. **Two genuinely new
  chunks:** (1) *the single hardest prerequisite* — extending the `ℂ((z))` valuation to
  `AlgebraicClosure ℂ((z))` (a Puiseux surrogate; buildable from `SpectralNorm`/`Krasner` on the
  complete field `ℂ((z))`, giving the `<1`/`>1` root dichotomy directly, but not off-the-shelf);
  (2) the bespoke **bi-infinite `ℂ[[X,X⁻¹]]` with the band-limited `L`** (Mathlib's Hahn/Laurent
  series are one-sided-bounded — the two-sided decay object must be built; elementary algebra, no
  analysis). **Recommended route** if DvdK is to be formalized at all.

## The shortcut, and why it does NOT rescue THM-2022
For a **concrete** `f` (fixed rational coefficients), `∃ m0, CT(f^{m0})≠0` is a *finite computation*
(days): clear denominators to `g=Xˢf ∈ ℤ[X]`, then `CT(fᵐ)=Polynomial.coeff (gᵐ) (s·m)`, exhibit the
first nonzero `m0` — the general theorem guarantees the search terminates; ESV's effective bound
(`≤ m+n` for doubly-monic, conjectural) says `m0` is small.
**But THM-2022's face Laurent polynomial has SYMBOLIC coefficients** — after the §2 descent the `c_i`
are the arbitrary algebraic coordinates of the hypothetical torus null point, universally quantified.
So `CT(f_F^{m0})` is a *polynomial in the `c_i`* and must be shown `≢ 0` **as a polynomial**, which is
precisely the general theorem. The finite-computation shortcut applies only to checking a *specific*
counterexample candidate, never to proving NC2. **Confirmed: THM-2022 needs general DvdK.**
No elementary general existence proof is known — both proofs route through the same core (`W(z)` is a
nontrivial algebraic function); dropping `limsup=|v|` to bare existence does not unlock a cheaper
proof. (One narrow easy case: if `f = g(z)·g(1/z)`-type, then `CT(fⁿ)=(1/2π)∫|g|^{2n} > 0` trivially —
but that does not cover general face Laurent polynomials.)

## Bottom line for the NC2/GMC2 programme
- **Formalizing DvdK is a standalone ≈4–9 person-month project** (Monsky/§5 route), with the valuation
  extension to `AlgebraicClosure ℂ((z))` as the one hard prerequisite. The analytic route is years.
- **Recommended for the NC2 formalization: keep DvdK as a cited hypothesis** (exactly as THM-1630 is
  cited on paper — like LRC≤13), i.e. state `nc2` with the DvdK constant-term nonvanishing as an
  explicit input, and formalize everything else. This is the honest, tractable division of labour: the
  arithmetic engine (§1 Wick / §4 multinomial-Lucas / §5 Frobenius + `face_sum_frobenius`, already
  kernel-pure this session) + face geometry + §2 descent are formalizable; DvdK is imported.
- If someone *does* want DvdK in Lean, this note is the blueprint: Monsky, vdE–Schoone §5, start with
  the valuation-on-`AlgebraicClosure ℂ((z))` module.

Cross-links: THM-1630 (DvdK citation, TNC identification), THM-2022 (the proof needing it), reflection
`formalizing-thm-2022-...-S94` (the engine), `GMC2Reduction.lean` (`face_sum_frobenius`/`face_sum_ne_zero`).
Sources: DvdK Indag. Math. 9 (1998) 221–231; Monsky arXiv:0906.1836; van den Essen–Schoone
arXiv:2305.10062 §5; Erman–Smith–Várilly-Alvarado arXiv:0908.2609.
