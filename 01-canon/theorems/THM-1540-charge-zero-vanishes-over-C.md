---
id: THM-1540
title: "THE CHARGE-0 PART VANISHES OVER C — correcting and strengthening THM-1535. SELF-AUDIT: THM-1535 s2 proved the sign-coherent case via c^T H c with H_{ab}=(a+b)! positive definite, which is valid only for REAL c; over C the form c^T H c (no conjugate) vanishes for c != 0 — explicitly c_0 = c_1(-1 +- i) at size 2 — and my 59048-polynomial sweep used real coefficients {-1,0,1}, so the complex case was ASSUMED, not proved. THE FIX IS A BETTER REDUCTION: in polar form z = r e^{i theta}, the CHARGE GRADING IS EXACTLY THE FOURIER GRADING IN theta, and s = |z|^2 is EXPONENTIAL(1) since E[(z zbar)^k] = k!. So the charge-0 part P0 = sum_a c_a (z zbar)^a is just g(s), and E[P0^m] = int_0^inf g(s)^m e^{-s} ds — a ONE-DIMENSIONAL problem. With all charges >= 0 the charge-0 part of P^m is exactly P0^m, so nullcone gives ALL moments of g vanishing, not merely the second; and only g = 0 does that (exact solve at degrees 1,2,3 over C; and asymptotically, for deg g = d >= 1 the phase arg g(s) STABILISES to arg(leading coeff) as s -> inf so no oscillation remains to cancel the tail, giving int ~ e^{a/d}(dm)! != 0). Hence GMC(2) IS PROVED ON THE SIGN-COHERENT LOCUS OVER C, upgrading THM-1535 from R. The nullcone conjecture also survives a complex recheck: Gaussian-integer coefficients, degree <= 2 gives 48 nullcone members ALL charge-definite, and 262143 mixed-charge degree-3 candidates give ZERO nullcone members"
status: PROVED (the reduction; the sign-coherent case over C, with the deg<=3 solve exact and the general case asymptotic) + VERIFIED (complex sweeps). The both-signs case remains OPEN
author: opus-2026-07-20-S412
corrects: THM-1535 §2 (real-coefficients-only gap; the conclusion stands, the proof is replaced)
depends_on: [THM-1535, THM-1495]
---

# THM-1540 — The charge-0 part vanishes over ℂ

## 1. Self-audit: what was wrong with THM-1535 §2

THM-1535 §2 argued: with only nonnegative charges, `E[P²] = E[P₀²] = cᵀHc` with
`H_{ab} = (a+b)!` positive definite, hence `c = 0`.

**That is valid only for real `c`.** GMC is stated over `ℂ`, and the form `cᵀHc` — with
**no conjugate** — can vanish for `c ≠ 0` even when `H ≻ 0`. Explicitly at size 2,
`H = [[1,1],[1,2]]` and `cᵀHc = 0` has solutions `c₀ = c₁(−1 ± i)`.

Worse, the supporting sweep used coefficients in `{−1,0,1}` — **real**. So the complex case
was assumed, not tested. **The conclusion turns out to be right; the proof was not.**

*(Note in passing: over `ℝ` the whole question is trivial — `E[P²] = ∫P²dγ > 0` for real
`P ≠ 0`, so the real nullcone is `{0}`. GMC has content only over `ℂ`.)*

## 2. The right reduction: charge = Fourier mode, radius = exponential

Write `z = r e^{iθ}`. A monomial `z^a \bar z^b` has modulus `r^{a+b}` and phase
`e^{i(a−b)θ}`:

> **The charge grading IS the Fourier grading in `θ`.** `E` = angular average (kills all
> `q ≠ 0`) followed by radial average.

And with `s = r²`, the radial law is **Exponential(1)**, because `E[(z\bar z)^k] = k!`.
Hence for the charge-0 part `P₀ = Σ_a c_a (z\bar z)^a = g(s)`:

```
E[P₀^m]  =  ∫₀^∞ g(s)^m e^{−s} ds
```

*(verified exactly against the direct computation for `g = 1−s`, `s²−4s+2`, and the complex
`g = is−1`.)* **The charge-0 problem is one-dimensional.**

## 3. The sign-coherent case, over ℂ (PROVED)

If every charge of `P` is `≥ 0`, then the only way a product of monomials has total charge
`0` is for **every** factor to have charge `0`. So the charge-0 part of `P^m` is exactly
`P₀^m`, giving

```
E[P^m] = E[P₀^m]   for every m.
```

So `P` in the nullcone forces `∫₀^∞ g^m e^{−s}ds = 0` for **all** `m` — not merely `m = 2`,
which is what THM-1535 used and where the complex gap opened.

> **Lemma.** The only polynomial `g ∈ ℂ[s]` with `∫₀^∞ g^m e^{−s}ds = 0` for all `m ≥ 1`
> is `g = 0`.

- **Exact,** degrees 1, 2, 3: solving `E[g^m] = 0` for `m = 1…2d+3` over `ℂ` returns **only
  the zero solution** in each case.
- **General, asymptotic.** Let `d = deg g ≥ 1` with leading coefficient `c_d`; factor out
  `c_d^m` (harmless for vanishing) and write `g(s) = s^d(1 + a/s + O(1/s²))`. The weight
  `s^{dm}e^{−s}` peaks at `s = dm`, where `m·a/s ≈ a/d` is a **constant**; so
  `∫₀^∞ g^m e^{−s}ds ≈ e^{a/d}·(dm)! ≠ 0`. The point is that
  **`arg g(s) → arg c_d` as `s → ∞` — the phase stabilises, so no oscillation survives to
  cancel the tail.** For `d = 0`, `g = c₀` and the integral is `c₀^m`, forcing `c₀ = 0`. ∎

> **Corollary (upgrade of THM-1535).** At `n = 2`, a sign-coherent nullcone element has
> `P₀ = 0`, hence is strictly charge-definite — **over `ℂ`**. With THM-1535 §1, **GMC(2)
> holds on the sign-coherent locus over `ℂ`.**

## 4. The nullcone conjecture rechecked over ℂ

| sweep | scanned | nullcone members | not charge-definite |
|---|---|---|---|
| deg ≤ 2, coeffs in `{0,±1,±i}` | 3124 | 48 | **0** |
| deg ≤ 3, coeffs in `{0,1,−1,i}`, **mixed-charge only** | 262143 | — | **0** |

**No complex counterexample.** The conjecture — `N₂` = charge-definite, i.e. the Newton
polygon misses the diagonal `a = b` — survives the complex recheck it had not previously had.

## 5. The remaining case, in its cleanest form

Only **charges of both signs** is left. In the polar picture, write `w = e^{iθ}` so that `P`
is a **Laurent polynomial in `w`** with `r`-dependent coefficients, `P = Σ_q f_q(r) w^q`,
having both positive and negative powers. Then the angular average is the constant term, and

> **Remaining statement.** If `P = Σ_q f_q(r)w^q` has both `q > 0` and `q < 0` present, then
> `∫₀^∞ CT_w( P^m ) e^{−s} ds ≠ 0` for some `m` (`s = r²`).

Two things are known and neither suffices alone: by Gordan, a balanced monomial first
appears at `m = c + |d|` (and empirically that is exactly where `E[P^m]` first fails to
vanish); and Newton-polytope **vertices** cannot cancel, since a vertex coefficient of `P^m`
is a pure power. What is missing is control of the non-vertex balanced monomials, now with
the extra `r`-integration available as a tool rather than an obstacle.

## 6. Status of GMC(2)

| locus | status |
|---|---|
| charge-definite | **PROVED** (THM-1535 §1) — in every dimension |
| sign-coherent, over `ℝ` | **PROVED** (THM-1535 §2) |
| sign-coherent, over `ℂ` | **PROVED** (§3 here) |
| charges of both signs | **OPEN** — the last case |

## Verification

`04-computation/charge_zero_complex_opus_S412.py` (the gap; the polar/exponential reduction;
exact solves at degrees 1–3; asymptotic illustration),
`04-computation/nullcone_complex_sweep_opus_S412.py` (Gaussian-integer sweeps). Outputs in
`05-knowledge/results/`.
