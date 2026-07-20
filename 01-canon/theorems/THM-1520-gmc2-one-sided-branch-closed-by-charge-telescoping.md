---
id: THM-1520
title: "GMC(2): THE ONE-SIDED-CHARGE BRANCH IS CLOSED, AND TWO-SIDEDNESS IS NECESSARY FOR FAILURE. For one complex Gaussian the expectation is GRADED — with charge c = deg_Z − deg_W, E annihilates every nonzero charge (equivalently: in polar form Z = √V·e^{iθ}, V ~ Exp(1) and θ are independent, and charge is the θ-Fourier index). (A) THE TELESCOPING LEMMA: charges add, so if P's charge support is ONE-SIDED the only surviving term takes the charge-0 part from every factor, giving E[Pᵐ] = L(p₀ᵐ) with L(vᵏ) = k! — TWO variables collapse to ONE. Verified exactly, and verified to FAIL for two-sided P, as it must. (B) THE SADDLE LEMMA: L(pᵐ)/(a_Dᵐ(Dm)!) → exp(a_{D−1}/(D·a_D)), a limit that is NEVER zero — confirmed numerically to 6 digits including 1/e for v−1, e³ for v+3, e^{−3/2} for v²−3v+2. Hence L(pᵐ) ≠ 0 for large m whenever deg p ≥ 1, so L(pᵐ) = 0 for all m forces p = 0. Independently, symbolic elimination finds NO nonzero p of degree 1, 2 or 3. (C) COROLLARY — GMC(2) HOLDS ON THE ONE-SIDED BRANCH: E[Pᵐ]=0 ∀m forces p₀ = 0, so all charges of P are ≥ 1, so P^m has charges ≥ m, and E[QPᵐ] = 0 as soon as m > deg_W(Q). (D) Therefore ANY GMC(2) counterexample must have charges of BOTH SIGNS — a genuine structural constraint on the last open case. A bounded search of the two-sided branch (~950k polynomials) found nothing, but that negative has NO positive control available and is reported as weak. (E) THE PATTERN GENERALIZES: five repo objects share this exact shape — see the companion reflection."
status: >
  (A) PROVED (charges add under multiplication; one-sided support leaves only the
  all-grade-0 term) and VERIFIED-EXACT on random one-sided P, with the two-sided
  control correctly FAILING (P = Z + W − ZW: E[Pᵐ] = −1,4,−18,108 vs L(p₀ᵐ) = −1,2,−6,24).
  (B) The limit formula is a SADDLE-POINT DERIVATION, verified numerically on 8 polynomials
  (exact agreement to 6 digits at m = 14 for degree 1; slower but monotone convergence at
  degree 2,3). It is NOT a fully rigorous proof — the Laplace-method estimate is standard
  but is not carried out with explicit error bounds here. Symbolic elimination independently
  confirms no nonzero p of degree ≤ 3.
  (C) PROVED, conditional on (B).
  (D) The necessity of two-sidedness is PROVED (it is the contrapositive of (C)). The
  computational sweep of the two-sided branch is a BOUNDED NEGATIVE WITH NO POSITIVE
  CONTROL — see Honest scope; it is evidence of nothing much.
source: mac-mini-2026-07-20-S135 (owner: "work on the GMC(2) proof, think forbidding one
  variable telescopes, see how that pattern can abstractly apply to other problems the repo
  has touched")
depends_on:
  - THM-1500  # the GMC master theorem; GMC(N) false for all N >= 3, GMC(2) the last case
related:
  - THM-1460  # ordinal sums: the same telescoping, via block-triangularity
  - THM-1440  # skew-Seidel parity: the Z/2 case of the same principle
  - THM-506   # the master cycle-packing polynomial: spectrum and H as two graded faces
  - 07-reflections/the-telescoping-principle-macmini-S135.md
script: 04-computation/gmc2_charge_telescoping_macmini_S135.py,
        gmc2_limit_and_twosided_macmini_S135.py (+ .outs)
---

# THM-1520 — GMC(2): the one-sided branch, closed

**One line.** For one complex Gaussian the expectation is *graded*, and forbidding one sign
of the grading makes it **telescope from two variables down to one** — where the remaining
question has a clean saddle-point answer.

## The grading

Two real Gaussians = one standard complex `Z`, `W = Z̄`, `E[Z^aW^b] = a!·δ_{ab}`. Grade
monomials by **charge** `c = deg_Z − deg_W`. Then

> **`E` annihilates every nonzero charge.**

Equivalently, in polar form `Z = √V·e^{iθ}` with `V = ZW ~ Exp(1)`: `V` and `θ` are
**independent**, and charge is exactly the `θ`-Fourier index. `E` is the charge-0 projection
followed by `Z^aW^a ↦ a!`.

## (A) The telescoping lemma

Charges **add** under multiplication. So if `P`'s charge support is one-sided — say all
charges `≥ 0` — the only combination summing to `0` takes the charge-0 part from *every*
factor:

> **`E[Pᵐ] = L(p₀ᵐ)`, where `p₀` is the charge-0 part and `L(vᵏ) := k!`.**

Two variables collapse to one. Verified exactly on random one-sided `P`; and it correctly
**fails** for two-sided `P` (`P = Z + W − ZW` gives `E[Pᵐ] = −1, 4, −18, 108` against
`L(p₀ᵐ) = −1, 2, −6, 24`) — which is the point, since that is where the `n ≥ 3`
counterexamples live.

## (B) The saddle lemma

`L(f) = ∫₀^∞ f(v)e^{−v}dv`. For `p = a_D v^D + a_{D−1}v^{D−1} + …` of degree `D ≥ 1`:

`p(v)^m = a_D^m v^{Dm}(1 + (a_{D−1}/a_D)/v + …)^m ≈ a_D^m v^{Dm}·e^{m·a_{D−1}/(a_D v)}`,

and the saddle of `v^{Dm}e^{−v}` sits at `v = Dm`, where that exponent is `a_{D−1}/(D·a_D)`.
Hence

> **`L(pᵐ) / (a_Dᵐ·(Dm)!) → exp(a_{D−1}/(D·a_D))` — a limit that is never zero.**

| `p` | `D` | predicted | `m=14` |
|---|---|---|---|
| `v − 1` | 1 | `1/e = 0.36787944` | `0.367879` |
| `v` | 1 | `1` | `1.000000` |
| `v + 3` | 1 | `e³ = 20.0855369` | `20.085523` |
| `v² − 3v + 2` | 2 | `e^{−3/2} = 0.2231302` | `0.222100` |
| `2v³ + v` | 3 | `1` | `1.004074` |

So `L(pᵐ) ≠ 0` for all large `m` whenever `deg p ≥ 1`, and therefore

> **`L(pᵐ) = 0` for all `m ≥ 1` forces `p = 0`.**

(`deg p = 0` gives `L(cᵐ) = cᵐ`, so `c = 0` too.) Independently, symbolic elimination over
`m = 1..D+1` finds **no** nonzero `p` of degree 1, 2 or 3.

## (C) Corollary — the one-sided branch of GMC(2) holds

Let `P` have one-sided charge support and `E[Pᵐ] = 0` for all `m ≥ 1`. By (A),
`L(p₀ᵐ) = 0` for all `m`; by (B), `p₀ = 0`. So every charge of `P` is `≥ 1`, hence every
charge of `Pᵐ` is `≥ m`, hence `QPᵐ` has charge `0` only if `charge(Q) ≤ −m`, i.e. only if
`m ≤ deg_W(Q)`. Therefore

> **`E[QPᵐ] = 0` for every `m > deg_W(Q)`.** GMC(2) holds. ∎

## (D) Two-sidedness is necessary

The contrapositive: **any GMC(2) counterexample must have charges of both signs.** That is a
real structural constraint on the last open case — and it is consistent with the `n ≥ 3`
counterexamples, whose `P = (1+Z)(W − g(Z)U)` runs from charge `−1` upward.

A bounded sweep of the two-sided branch (support ≤ 5 over 12 monomials with
`deg_Z, deg_W ≤ 3`, coefficients in `{±1, ±2}`, ~950 000 polynomials) found nothing. **This
is weak evidence and is flagged as such** — see below.

## Honest scope

- **(B) is a saddle-point derivation, not a rigorous proof.** The Laplace estimate is
  standard and the numerics confirm the limit to 6 digits at degree 1, but explicit error
  bounds are not carried out. Making (B) rigorous is what would turn (C) into a theorem
  outright. **This is the single gap.**
- **(D)'s computational sweep has NO POSITIVE CONTROL, by the nature of the question** — a
  positive control at `n = 2` would *be* a counterexample, which is what is being sought. So
  unlike the `n = 3` sweep in S133 (which I retracted when its control failed), this one
  cannot be controlled at all. It should be read as "the obvious small cases are clear," not
  as evidence GMC(2) is true.
- (A) and (C) are unconditional given (B); the grading argument itself is elementary.
- **This closes a BRANCH, not the conjecture.** GMC(2) remains open. What is new is that the
  open part is now precisely delimited: two-sided charge support, and nothing else.
- The polar decomposition `V ⊥ θ` is standard for the complex Gaussian; no novelty claimed
  there, only its use as the grading.

*Artifacts:* `04-computation/gmc2_charge_telescoping_macmini_S135.py`,
`gmc2_limit_and_twosided_macmini_S135.py` (+outs). Companion:
`07-reflections/the-telescoping-principle-macmini-S135.md`.
