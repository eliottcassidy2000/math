---
id: THM-1530
title: "THE TORAL NULLCONE, AND AN EXACT LAGRANGE PROOF AT EXTREME WEIGHT ±1. (A) NC2's leading order is a purely algebraic question: writing Z = √r·e^{iθ} with r ⊥ θ, E[P^m] ~ Γ(Dm+1)·CT_u(Λ(u)^m) where Λ is the LEADING SYMBOL, a Laurent polynomial. So the nullcone question becomes the TORAL NULLCONE (TNC): which Λ have CT(Λ^m) = 0 for every m ≥ 1? (B) TNC IS PROVED WHEN THE MINIMUM EXPONENT IS −1 (or the maximum is +1), exactly and without asymptotics: Λ = u^{−1}R(u) gives CT(Λ^m) = [u^m]R(u)^m, and Lagrange–Bürmann turns Σ_m t^m[u^m]R^m = 1/(1 − tR'(w(t))) with w = tR(w); demanding vanishing forces R'(w(t)) ≡ 0, and w is nonconstant since w'(0) = R(0) ≠ 0, so R' vanishes on an open set, so R is CONSTANT and Λ is a single weight. (C) THE CONTRAST AT MINIMUM EXPONENT −M, M ≥ 2 IS STRUCTURAL, NOT MERELY HARDER: branching over the M-th roots of unity with H = R^{1/M} gives the condition [u^m]H^m = 0 only for m divisible by M — a strictly weaker demand — so the M ≥ 2 case is a different problem. Tested: 3,456 three-term Laurent polynomials with min ≤ −2 and max ≥ 2, zero hits. (D) THE EXACT GENERATING FUNCTION for the {−1,0,1} stratum, derived by residues and verified: Σ_m L_m(α,β)t^m = ((1−βt)² − 4αt²)^{−1/2}, so Σ_m E[P^m]t^m = E_r[((1−β(r)t)² − 4α(r)t²)^{−1/2}] with α = r·a·c and β = b. At constant α, β this forces β = 0 and α = 0 by hand — TNC at M = N = 1. (E) The r-dependent stratum: 15,576 triples searched, ZERO witnesses; every case has E_r[L_m] ≠ 0 at large m, the sole exception being the parity-forced odd-m vanishing at β ≡ 0, which is the two-weight {−1,+1} case already closed by THM-1510."
status: >
  (A) DERIVATION (leading-order Laplace); the reduction is heuristic in the sense that the
      subleading terms are not controlled here.  Stated as the frame, not as a theorem.
  (B) PROVED, exactly, by Lagrange-Buermann.  No asymptotics, no genericity.  Machine-checked
      exhaustively over small R: zero counterexamples.
  (C) PROVED (the root-of-unity branching is a computation); the M >= 2 case is OPEN and is
      tested, not settled -- 3,456 cases, zero hits.
  (D) PROVED (residue computation) + VERIFIED against direct expansion on four (α,β) pairs.
  (E) BOUNDED SEARCH, NOT A PROOF.  Box stated explicitly.
  GMC(2) REMAINS OPEN.  This adds one exactly-proved case and a cleaner frame; it does not
  close the conjecture.
source: klein-2026-07-20-S345 (owner: work to prove GMC(2), think nullcone, keep up with concurrent work and extend not duplicate)
concurrency: >
  THIS FILE EXTENDS AND DOES NOT DUPLICATE.  Explicitly CEDED to concurrent work:
  - mac-mini THM-1500 (pushed BEFORE my THM-1510) already had the master theorem AND the
    stronger statement that f = f(0)/(1+s) is UNIQUELY FORCED by Lagrange inversion.  My
    THM-1510 §A/§B duplicated it with only sufficiency; theirs is better and has priority.
  - mac-mini THM-1520 independently derived the same Laplace constant exp(a_{D-1}/(D a_D))
    that my THM-1510 EMP uses, and applied it to the one-sided-charge branch.  My EMP was
    pushed first; their extension is theirs.
  - boxeph HYP-8320/8350: the radial-family no-go, and the Fock/vacuum-Mathieu bridge.
  - death-star HYP-8330: GMC(2) evidence (532 kernel elements, zero counterexamples).
  NEW HERE: the toral-nullcone framing (A), the Lagrange proof at extreme weight ±1 (B), the
  M >= 2 structural contrast (C), and the square-root generating function (D).
depends_on:
  - THM-1510  # EMP, NC2, and the two-weight theorem -- this continues that line
related:
  - THM-1500  # mac-mini's master theorem (ceded)
  - THM-1520  # mac-mini's one-sided-charge branch
  - THM-1490  # the weight-polytope obstruction
script: 04-computation/toral_nullcone_klein_S345.py (+ .out)
---

# THM-1530 — the toral nullcone

Setting as in THM-1510: `n = 2`, one complex Gaussian `Z`, `E[Z^aZ̄^b] = δ_{ab}a!`, weight of
`Z^aZ̄^b` is `a−b`, and `NC2` is the conjecture that the nullcone
`N_2 = {P : E[P^m] = 0 ∀m ≥ 1}` contains only the one-signed-weight polynomials. Recall
**NC2 ⟹ GMC(2)** (THM-1510 §C).

## A. The frame: NC2's leading order is toral

Polar-decompose `Z = √r·e^{iθ}` with `r ~ Exp(1)` **independent** of `θ ~ Unif`. Then `P` is
literally a function of `(r,θ)` and `E[P^m] = E_{r,θ}[P^m]` — the `θ`-average *is* the weight-0
projection. Writing `u = e^{iθ}` and `P(r,u) = Σ_q g_q(r)u^q`, let `D = max_q deg g_q` and

```text
Λ(u) := Σ_{q : deg g_q = D} c_q u^q        (the LEADING SYMBOL, a Laurent polynomial)
```

To leading order `E[P^m] ~ Γ(Dm+1)·CT_u(Λ(u)^m)`. So NC2 demands, at leading order:

> **Toral nullcone (TNC).** For which Laurent polynomials `Λ` is `CT(Λ^m) = 0` for every
> `m ≥ 1`? Conjecturally: exactly those with all exponents of one strict sign.

This strips the analysis down to an algebraic question with no Gaussian in it.

## B. TNC is PROVED at extreme weight ±1 — exactly

**Theorem.** If the minimum exponent of `Λ` is exactly `−1` and `CT(Λ^m) = 0` for all `m ≥ 1`,
then `Λ = c·u^{−1}` (a single weight).

*Proof.* Write `Λ(u) = u^{−1}R(u)` with `R` a polynomial, `R(0) ≠ 0`. Then `Λ^m = u^{−m}R^m`, so

```text
CT(Λ^m) = [u^m] R(u)^m.
```

Lagrange–Bürmann: for `w = tR(w)` (unique small root, analytic, `w'(0) = R(0) ≠ 0`),

```text
Σ_{m≥0} t^m [u^m]R(u)^m  =  1/(1 − t·R'(w(t))).
```

Vanishing for every `m ≥ 1` makes the left side identically `1`, so `t·R'(w(t)) = 0`, i.e.
`R'(w(t)) = 0` for `t ≠ 0`. Since `w'(0) = R(0) ≠ 0`, `w` is nonconstant and sweeps a
neighbourhood of `0`; a polynomial vanishing on an open set is zero, so `R' ≡ 0`, `R` is
constant, and `Λ = c·u^{−1}`. ∎

By `u ↦ 1/u` the same argument settles *maximum exponent exactly `+1`*. **No asymptotics and
no genericity** — this is the one piece of the GMC(2) programme that is exact.

*Machine check:* exhaustive over `R` of degree `≤ 3` with coefficients in `{0,±1,±2,±i}` —
**zero** Laurent polynomials with min exponent `−1`, more than one term, and `CT(Λ^m) = 0` for
`m = 1..8`.

## C. `M ≥ 2` is a different problem, not a harder one

For minimum exponent `−M`, `Λ = u^{−M}R(u)` and `CT(Λ^m) = [u^{Mm}]R^m`. Substituting
`t = τ^M` and branching over the `M`-th roots of unity, with `H = R^{1/M}`, the same Lagrange
computation gives the condition

```text
[u^m] H(u)^m = 0   only for   M | m,
```

a **strictly weaker** demand than the `M = 1` case, where it is required for *every* `m`. So
the obstruction that closes `M = 1` genuinely does not transfer, and `M ≥ 2` must be attacked
differently. This is worth recording precisely because it looks like a routine generalisation
and is not.

*Bounded evidence:* 3,456 three-term Laurent polynomials with `min ≤ −2` and `max ≥ 2`, none
with `CT(Λ^m) = 0` for `m = 1..10`. (Note the trivial escape: if every exponent is divisible
by `g > 1` then `Λ(u) = μ(u^g)` and nothing is new — such hits would be rescalings.)

## D. The square-root generating function on the `{−1,0,1}` stratum

For `P = Z̄a(r) + b(r) + Zc(r)`, put `α = r·a(r)c(r)` and `β = b(r)`. A residue computation on
`CT_w(β + x(w + w^{−1}))^m` with `x² = α` gives the exact ordinary generating function

```text
Σ_m L_m(α,β) t^m  =  ((1 − βt)² − 4αt²)^{−1/2},      L_m = Σ_k m!/(k!²(m−2k)!) α^k β^{m−2k}
```

hence `Σ_m E[P^m]t^m = E_r[((1 − β(r)t)² − 4α(r)t²)^{−1/2}]`. *Verified* against direct
expansion at `(α,β) = (2,3), (1,0), (0,5), (−1,2)`.

At **constant** `α, β` setting this to `1` gives `(1−βt)² − 4αt² = 1`, i.e. `β = 0` and
`α = 0` — TNC at `M = N = 1`, by hand. Note the exponent `−1/2`: the same square root that
boxeph identified as the repo's Wallis/fiber-fraction family, arriving here from residues
rather than from the `U`-side.

## E. The `r`-dependent stratum — searched, not proved

15,576 triples `(a,b,c)` linear with integer coefficients in `[−2,2]`: **zero witnesses** with
`E_r[L_m] = 0` for `m = 1..14`. Exact arithmetic (`L_m(α(r),β(r))` is a polynomial in `r` and
`E_r[r^k] = k!`).

**Self-correction.** My run's summary line reported "no case has `E_r[L_m]` vanishing at large
`m`: **False**". That was a **parity artifact of my own flag**, not a witness: the offending
case is `β ≡ 0`, where `L_m` contains only `β^{m−2k}` and therefore vanishes identically for
*odd* `m`. Its even moments are `2, 12, 120, 1680, 30240, 665280 = (2k)!/k!`, all nonzero — and
`β ≡ 0` is the two-weight `{−1,+1}` case already closed by THM-1510 via EMP. So the correct
statement is: **every tested case has `E_r[L_m] ≠ 0` at large `m`**, with the only zeros being
parity-forced.

## F. Where GMC(2) now stands

| stratum | status |
|---|---|
| one-signed weights | trivially in the nullcone; NC2 permits |
| at most two weights, any degree | **proved** (THM-1510, EMP) |
| leading symbol with extreme weight `±1` | **proved** (§B, exact) |
| leading symbol with `min ≤ −2` and `max ≥ 2` | **open** (§C explains why) |
| full `{−1,0,1}` with `r`-dependence | searched clean; Laplace route sketched, not carried |

The honest remaining obstacles are (i) the `M,N ≥ 2` toral case, and (ii) promoting the
leading-symbol reduction (§A) from a leading-order heuristic to a controlled induction. Both
are now sharply stated.

*Files: `04-computation/toral_nullcone_klein_S345.py` (+ `.out`).*
