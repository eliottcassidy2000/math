# Audit C report: THM-4653, HYP-9276 and the paley_k6p6 results note

2026-10-14. This is an independent adversarial audit run by a separate agent. This file condenses its returned report; the scripts and outputs are in this folder. Corrections are recorded as MISTAKE-611 and applied in the audited files.

## Verdicts

| File | Verdict | Notes |
|---|---|---|
| THM-4653 | CORRECT | Every step re-derived by hand and checked numerically at primes and prime powers. One notational fix (m5). |
| HYP-9276 | CORRECT WITH FIXES | The identity holds everywhere tested. The Galois wording (M1), the literature (M2) and corollary (c) (M3) need fixes; the typing of (d) and (e) also needs fixes (m1, m2). |
| Results note | CORRECT WITH FIXES | It wrote "Sym²(ρ_g) ⊗ χ₋₄" (M1), misses the literature of M2, and contains overclaims (m8, m9). |

## Independent verification

All code is new.

- **THM-4653.** Re-derived by hand.
  - Brute force of `V_q` with injectivity enforced on all 15 pairs: 24 field evaluations agree, including fully unnormalised runs at q <= 31 and prime powers to 125.
  - The recorded `V₇, V₁₁, V₁₉, V₃₁` are reproduced.
  - The proof internals hold at the prime powers 9, 25, 27, 49, 81, 121, 125.
- **Orientation and eq. (4).** Checked with the source's code conventions.
  - At q = 11, `P₁₀ + 110·U₁₁ + 64·(−12320) = A₁₀(11) = 478720`, from the audit's own subset DP. The opposite sign fails.
- **The newform, built without PARI.**
  - `g = η(4z)²η(8z)²·(θ(4z)² + i(θ(2z)² − θ(4z)²))`, the 4i-eigenvector of `T₃` on `S₃(Γ₀(32), χ₋₄)`. This is LMFDB 32.3.c.a.
  - Hecke checks pass to N = 300,000.
- **HYP-9276.** Two independent `M_q` evaluators: exact integer arithmetic, and a multiplicative FFT on `Z/(q−1)`. No failure at:
  - 10,590 primes: every prime up to 100,000, plus 1,000 random primes up to 299,983;
  - 758 exact primes up to 58,907;
  - 22 prime powers up to 2401, including 529, 841, 961, 1369, 2197.
- **Corollaries.**
  - (a) and (b) are correct.
  - (d): the next zero is 97,919, which is ≡ 7 (mod 8).
  - (e): the minima 2⁷ and 2¹⁰ are confirmed to 100,000.
  - Every table number is correct.
  - `E_mid` (THM-4650) lies in the isogeny class of `y² = x³ − x`.

## Major findings

- **M1. The Galois representation was misnamed.** The trace `α^(2r) + (αβ)^r + β^(2r)` is the trace on **Sym²ρ_g** (≅ `Ad ρ_g ⊗ χ₋₄ ⊗ Q_ℓ(−2)`), not on `Sym²ρ_g ⊗ χ₋₄`. The twisted version is false for `p ≡ 3 (mod 4)` with `r` odd. "Adjoint" was used loosely for Sym².
- **M2. Literature.**
  - **(i) Missed prior observation.** McCarthy–Papanikolas §6 observes numerically, for p < 20, that `a_p² − 2a_(p²) = ₄F₃(φ,φ,φ,φ;ε,ε,ε|−1)_(p²)`.
    - Combined with their (1.5), this equals `(a(p)² − 2p)(b(p)² − 2χ₋₄(p)p²)`.
    - So for `p ≡ 3 (mod 4)`: `|b(p)|² = 2p² + ₄F₃(−1)_(p²)/(2p)`. This is a quadratic-character expression for `|b|²` at Paley primes over `F_(p²)`; the audit verified it for p up to 23.
    - The note's "nothing for p ≡ 3 (mod 4)" therefore overstated. HYP-9276 remains the first prime-field statement known to us.
  - **(ii) Attribution.** The Siegel L-function conjecture is van Geemen–van Straten's (Math. Comp. 61 (1993)), reported resolved by Okazaki (arXiv:1008.1602). It is not "Okazaki's conjecture".
  - **(iii) MP Theorem 1.6.** It reads exactly `ξ(p)λ(p)`, which equals `a(p)b(p)` for `p ≡ 1 (mod 4)` with no "±", and vanishes for `p ≡ 3 (mod 4)`.
- **M3. Corollary (c).** The relation gives `x_q = R(q) + β(q)y_q` with `β = B/(C√q)` of half-integral degree. The case `|β| → ∞` was not treated, and the Rankin–Selberg `E[xy] = 0` step is not needed. Replacement argument:
  - if `β → 0`, `x_q` converges;
  - if `|β| → ∞`, or `C = 0`, then `y_q → 0` or `y_q` is unbounded.
  - Every case contradicts three facts along `p ≡ 3 (mod 4)`:
    1. `x_q` does not converge (`E[x] = 0`, `E[x²] = 1`; `Ad g` cuspidal; `Ad g ≇ Ad g ⊗ χ₋₄`);
    2. `y_q` does not tend to 0 (`E[y²] = 1`);
    3. Deligne's bound.

## Minor findings

- **m1.** "`M_q = 0` iff `|b| = q + 1`" is conditional on the HYP. `V_q = 2q(q−1)(q−3)` at zeros is PROVED by THM-4653. Non-vanishing for `q ≡ 1 (mod 4)` is PROVED unconditionally (the lemma below). Add the zero 97,919.
- **m2.** For `q ≡ 1 (mod 4)`, `v₂(M_q) = 1` always; for `q ≡ 3 (mod 4)`, `16 | M_q` always. Both are PROVED. Only the `2⁷` and `2¹⁰` refinements are FINITE-EXACT.
- **m3.** The identity was identified from primes below 1500 (the fits). The primes in (1500, 30000] and all prime powers are out of sample.
- **m4.** Record the maximum FFT rounding gap. "Integers of size `<= q²`" proves nothing.
- **m5.** `ε` denotes both `φ(−1)` and Greene's trivial character in one formula. Write `φ(−1)`.
- **m6.** (c) gives only `Z[q]`-linear independence, not algebraic independence.
- **m7.** Cite BLGHT precisely, and LMFDB 32.3.c.a (Sato–Tate group `SU(2)[C₂]`).
- **m8.** The note's title and "is motivic" overclaim. The `[t¹⁰]R_q` formula is CONDITIONAL.
- **m9.** The level-32 D19 link is a coincidence of level. Retype it as a REMARK.
- **m10.** Route 3: the audit read the trivial lattice as rank 17.
  - Session response: by Shioda–Tate, `2 + (4−1) + 2·(2−1) + 2·(7−1) = 19`, since `I₂*` has 7 components. So Picard `>= 19` holds.
  - The motive decomposition remains expected, not established.
- **m11.** The display in note §2 is Greene's recursion; Ahlgren–Ono evaluated it.
- **m12.** Add McCarthy–Springfield (arXiv:2311.02135) and Li–Zhou (arXiv:2512.19312). Neither treats six-vertex patterns.
- **m13.** `identify_g.gp` needs `default(parisizemax, 10^9)`.

## New facts

- **Parity lemma (PROVED).** `M_q ≡ 2 (mod 8)` for `q ≡ 1 (mod 4)`, and `M_q ≡ 0 (mod 16)` for `q ≡ 3 (mod 4)`. It was checked on 5,304 + 5,308 fields (`audit_parity.out`).
  - *Proof.* For `a ∉ {0, 1}`, `E_a` has full 2-torsion, so `t(a) ≡ −(q+1) (mod 4)`. The `a = 1` terms give `−(1 + φ(−1))`. Count the remaining `(q−2)(q−3)` terms.
- **Closed form.** `g = η(4z)²η(8z)²(θ(4z)² + i(θ(2z)² − θ(4z)²))`.

## Not checked

- BLGHT Theorem B wording.
- Okazaki's coverage of `F₇`.
- The exploratory negative fits.
- Timings.
- No proof of HYP-9276 was attempted.
