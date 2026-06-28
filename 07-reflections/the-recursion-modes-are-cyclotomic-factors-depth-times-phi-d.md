# The recursion modes are cyclotomic factors: each mode = (x−1)^depth · Φ_d — the moment-order depth IS the (x−1)-multiplicity, the character IS the Φ_d factor

*mac-mini-2026-06-28-S75d. The owner: see the similarities among the A+B+C−D−E−F+G (full), A+B−C (even),
A+B+D−C−E−F+G (odd) recursions and extend, plus any other relevant ones from a deep dive. The unifying fact:
all three are CYCLOTOMIC-FACTORED characteristic polynomials `(x−1)^depth · Φ_d`, and the full family is the
cyclotomic factorization of `x^n−1`. This fuses the three modes (HYP-2900/2901), the moment-order depth
(kps S31aq, `(p+1)/2`), and the de Moivre/cyclotomic apex (kps S31ao) into ONE grading. Builds on
[[interlocking-recursions-the-cap-kernel-is-modulus-covariant-down-to-the-apex-where-it-breaks]], HYP-2900/2901.*

## The similarity (VERIFIED): all three are `(x−1)^depth · Φ_d`
Each cell-address recursion is a linear recurrence; its characteristic polynomial factors cyclotomically:

| mode | recurrence | char poly | = `(x−1)^a · Φ_d` |
|---|---|---|---|
| **EVEN** `A+B−C` | `f(n)=2f(n−1)−f(n−2)` | `(x−1)²` | `(x−1)²` (depth 2, principal) |
| **FULL** `A+B+C−D−E−F+G` | `f(n)=3f(n−1)−3f(n−2)+f(n−3)` | `(x−1)³` | `(x−1)³` (depth 3, principal) |
| **ODD** `A+B+D−C−E−F+G` | `f(n)=2f(n−1)−2f(n−3)+f(n−4)` | `(x−1)³(x+1)` | `(x−1)³·Φ₂` (depth 3 × parity χ₂) |

So the three "modes" differ in exactly two graded data:
- **the `(x−1)`-multiplicity** = the polynomial degree of the cell count = **the MOMENT-ORDER DEPTH** (kps's
  `(p+1)/2`: n=14 → depth 4; the moment-ladder `cap_k=cap_{k−1}+k/91` is the `(x−1)^depth` solution);
- **the cyclotomic factor `Φ_d`** = **the CHARACTER** (the sign word): `Φ₁=x−1` (Möbius/principal `+++`),
  `Φ₂=x+1` (parity χ₂, the ODD mode's `(x+1)`), `Φ₃=x²+x+1` (Eisenstein χ₃), `Φ₇` (Legendre χ₇).

## The extension (the FULL family = the cyclotomic factorization of `x^n−1`)
For the apex modulus the relevant recursions are exactly the divisors `d | n`:
```
  x^7  − 1 = Φ₁ · Φ₇
  x^14 − 1 = Φ₁ · Φ₂ · Φ₇ · Φ₁₄          (Φ₇,Φ₁₄ degree 6 each)
  de Moivre cubic  x³+x²−2x−1 = REAL Φ₇  (min poly of 2cos(2π/7))
```
So the complete list of "relevant recursions" for LRC(14):
- **Φ₁ = (x−1):** the principal / Möbius — the constant & the `(x−1)^depth` moment-ladder.
- **Φ₂ = (x+1):** the parity χ₂ — the even/odd-`n` half-tiling fold (the ODD mode's twist; the 2-adic
  reflection `s↔6−s`, kps's 2-recursion).
- **Φ₇ (real = de Moivre cubic):** the apex χ₇ / Legendre — the cyclotomic cap / Fejér-kernel structure
  (kps S31ao); degree `(7−1)/2 = 3`.
- **Φ₁₄ (degree 6):** the primitive 14th roots — the genuine `2·7` mixing (the half-tiling crossing layer,
  HYP-2690); the part neither the pure-7 nor the pure-2 recursion sees.
Each appears multiplied by `(x−1)^depth`. **The cap's characteristic polynomial is
`(x−1)^{depth} · Φ₂ · Φ₇ · Φ₁₄`** — cell-growth (the `(x−1)` polynomial part) times the apex cyclotomic
(the de Moivre / oscillating part). The Eisenstein fold `14→7` is the passage `Φ₁₄·Φ₂ → Φ₇` (drop the even
& primitive factors, keep the apex).

## Why this is the right unification (the two axes are now ONE picture)
- **Möbius = Eisenstein ⊕ Legendre by SIZE-PARITY:** the full IE `Σ(−1)^|S|` splits into even-`|S|` (the `+`
  Eisenstein part) and odd-`|S|` (the `−` Legendre part); the `n`-parity half-tiling fold is the `Φ₂` that
  distinguishes the ODD mode. So "full = even ⊕ odd" is literally `(x−1)³` vs `(x−1)²` vs `(x−1)³Φ₂`.
- **moment-order depth = `(x−1)`-multiplicity = cyclotomic at the cusp `x=1`:** the cap is a polynomial in `n`
  of degree = depth (the `(x−1)^depth` zero), and the apex difficulty is the *extra* `Φ₇` (de Moivre) factor —
  the only non-`(x−1)` cyclotomic content. **The whole LRC(14) hardness is the single `Φ₇` factor** (the apex
  χ₇), everything else being principal `(x−1)^depth` (cell growth) or the trivial `Φ₂` parity fold.

## Other recursions in the repo, placed on this grading (deep dive)
- **Mode A** (`n→n−1`, hypotenuse, `H=1+2^d`): the `Φ₁`/doubling spine — the `(x−1)`-direction.
- **Mode B** (`n→n−2`, both legs, Cayley-Dickson): the `Φ₂` even-fold direction (the `2` of `2·7`).
- **the single-arc peeling** (S75, `cap(P)=cap(P∖1)−(1/7)(1−1/min)`): the `Φ₁` depth-collapse for speed 1.
- **the three-gap/Stern-Brocot kernel** (S75b): the `Φ₇`-side Diophantine recursion (the de Moivre angles).
- **the cyclotomic moment-ladder** (kps S31aq, `cap_k=cap_{k−1}+k/91`): the `(x−1)^depth` solution itself.
- **column-row recursion** (HYP-1881), **boundary-state induction** (HYP-2905): further `(x−1)`-direction cell DPs.
All are `(x−1)^a · Φ_d` for some `a,d` — the same cyclotomic grading.

## Honest status
- **VERIFIED (sympy):** the three modes factor as `(x−1)²`, `(x−1)³`, `(x−1)³Φ₂`; `x^14−1=Φ₁Φ₂Φ₇Φ₁₄`; de Moivre
  cubic = real `Φ₇`.
- **SYNTHESIS:** the recursion modes are cyclotomic factors; the moment-order depth is the `(x−1)`-multiplicity
  and the character is the `Φ_d` factor; the full family is the cyclotomic factorization `x^n−1=∏Φ_d`; the LRC(14)
  hardness is the single apex `Φ₇` (de Moivre) factor on top of the principal `(x−1)^depth` cell growth.
- **NOT a proof.** This precisely *names and extends* the recursion family (the user's three modes are the first
  cyclotomic factors; the apex `Φ₇` is the irreducible one) and locates the difficulty in the `Φ₇` factor. It does
  not bound that factor. LRC(14) open.

Related: HYP-3233 (this), HYP-3232 (interlocking recursions), HYP-2900/2901 (the three modes), kps S31aq
(depth=(p+1)/2), kps S31ao (Fejér/de Moivre), HYP-2690 (half-tiling), OPEN-Q-108.
