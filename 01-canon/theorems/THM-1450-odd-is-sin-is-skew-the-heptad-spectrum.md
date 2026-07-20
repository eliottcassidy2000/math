---
id: THM-1450
title: "ODD = SIN = SKEW: x(x^2+7)(x^4+14x^2+17) IS THE MODAL 7-TOURNAMENT SPECTRUM. All its roots are purely imaginary, the signature of a real skew matrix; and by the switching-invariance of the characteristic polynomial plus THM-474 (switching classes = tilings) only 2^15 = 32768 representatives need scanning. Result: exactly ELEVEN distinct characteristic polynomials exist among all 32768 switching classes of 7-tournaments, and P = x^7+21x^5+115x^3+119x is the MOST COMMON at 10080 classes (30.8%), while Paley's x(x^2+7)^3 is RARE at 240 (0.73%). So P is the GENERIC heptad spectrum and Paley is the extremal one. REFUTED en route: P is NOT the sin(7t) multiple-angle polynomial (coeffs [-64,0,112,0,-56,0,7] vs [1,0,21,0,115,0,119]). CONFIRMED: x^2+7 = 0 gives +-i*sqrt7 = +-g(chi_7), the Gauss sum; the quartic has Galois group D_4; and P = x*(u^3 - 32u) with u = x^2+7, a DEPRESSED cubic with roots {0, +-4sqrt2}. THE UNIFICATION: skew (odd) <-> imaginary spectrum <-> sin/cos rotation; symmetric (even) <-> real spectrum. The skew/symmetric split of THM-1445 IS the sin/cos split, and 'odd-dimensional skew matrices are singular' is the SAME parity argument as THM-1445-A step 3 — an involution on an odd set must fix a point"
status: VERIFIED-BY-EXHAUSTION (all 32768 switching classes, exact integer Faddeev-LeVerrier) + symbolic checks
author: opus-2026-07-20-S407
depends_on: [THM-474 (mac-mini: switching classes = tilings — the 64x reduction), THM-1445 (parity eigenspaces), THM-1440]
---

# THM-1450 — Odd is sin is skew: the heptad spectrum

## 1. What the polynomial is

`P(x) = x(x²+7)(x⁴+14x²+17) = x⁷ + 21x⁵ + 115x³ + 119x` is **odd**, and every root is
**purely imaginary**: `0, ±i√7, ±1.15895i, ±3.55765i`.

Purely imaginary spectrum with a zero is the signature of a **real skew-symmetric matrix**.
A `7×7` skew matrix has spectrum `{0, ±iλ₁, ±iλ₂, ±iλ₃}`, so its characteristic polynomial
is odd of degree 7 with exactly this shape. **`P` is a tournament spectrum.**

Here `λ₁, λ₂, λ₃ = √7, √(7−4√2), √(7+4√2)`, and `λ₁²λ₂²λ₃² = 7·(49−32) = 7·17 = 119 = c₁`.

## 2. The census (exhaustive)

The characteristic polynomial is a **switching invariant**: switching is `S ↦ DSD` with
`D = diag(±1)`, a *similarity*, so the spectrum is unchanged. By **THM-474** the switching
classes of tournaments on `[n]` are exactly the tilings — the tournaments containing the
base path — so only `2^C(6,2) = 32768` representatives need scanning instead of `2²¹`. **A
64× reduction, taken directly from canon.**

| characteristic polynomial | `c₃` | `c₁` | switching classes | share |
|---|---|---|---|---|
| **`P = x(x²+7)(x⁴+14x²+17)`** | **115** | **119** | **10080** | **30.8 %** |
| — | 131 | 231 | 5040 | 15.4 % |
| — | 83 | 23 | 5040 | 15.4 % |
| — | 99 | 71 | 5040 | 15.4 % |
| — | 115 | 183 | 1680 | 5.1 % |
| … | | | | |
| **Paley `x(x²+7)³`** | 147 | 343 | **240** | **0.73 %** |

**Only ELEVEN distinct characteristic polynomials exist across all 32768 switching
classes** — a severe invariant collapse. `c₅ = 21 = C(7,2)` is forced for every 7-tournament
(it counts arcs), so all the information sits in `(c₃, c₁)`.

**So `P` is the GENERIC heptad spectrum — the modal one — and Paley's `x(x²+7)³` is the
RARE, extremal one.** That is the honest headline: this polynomial is what a *typical*
7-tournament sounds like, and the repo's beloved Paley heptagon is the outlier.

## 3. Anatomy, confirmed and refuted

**CONFIRMED.**
- `P(x) = x·(u³ − 32u)` with `u = x²+7`. The inner cubic is **depressed**, roots
  `{0, +4√2, −4√2}` — the **`{0, r, −r}`** shape: one fixed point plus one free 2-orbit of
  the involution `u ↦ −u`. The same shape as the JC fibre over a `τ`-fixed target
  (THM-1440/1445).
- `x² + 7 = 0` gives `±i√7`, and `i√7 = g(χ₇)`, the **Gauss sum** of the quadratic character
  mod 7 (verified numerically to 20 digits). This is the `p ≡ 3 (mod 4)` imaginary Gauss sum
  that kind-pasteur's heptagon reflection identifies.
- `x⁴+14x²+17` has Galois group **`D₄`** (order 8), computed.

**REFUTED.** `P` is **not** the `sin(7θ)` multiple-angle polynomial. Those coefficients are
`[−64, 0, 112, 0, −56, 0, 7]`; `P`'s are `[1, 0, 21, 0, 115, 0, 119]`. The degree and parity
match, the polynomial does not. Recorded so the next session does not re-guess it.

## 4. The unification the polynomial is pointing at

**Odd/even is sin/cos is skew/symmetric.** These are three names for the `±1` eigenspaces of
one involution:

| | ODD (`−1`) | EVEN (`+1`) |
|---|---|---|
| functions | `sin` | `cos` |
| matrices | **skew** `Sᵀ = −S` | **symmetric** `Sᵀ = S` |
| spectrum | purely **imaginary** `±iλ` | purely **real** |
| `exp(tS)` | **rotation** (orthogonal) — `cos λt`, `sin λt` | stretching |
| repo object | tournament `A − Aᵀ` | Seidel `J − I − 2A` |
| census | skew two-graphs `1,2,2,6` (THM-1415) | two-graphs `2,3,7,16,54` (THM-1430) |

A skew matrix is an *infinitesimal rotation*: `exp(tS)` is orthogonal, and it acts on each
eigenplane as `[[cos λᵢt, −sin λᵢt],[sin λᵢt, cos λᵢt]]`. **So the three rotation rates of a
generic 7-tournament are `√7, √(7−4√2), √(7+4√2)`** — that is what `P` encodes. Paley's
tournament rotates at the *single* rate `√7` in all three planes, which is exactly why its
spectrum degenerates to `x(x²+7)³`: it is the **isoclinic** case.

## 5. The parity argument is the same one, again

A skew matrix of **odd** size is singular: `det S = det(−Sᵀ) = (−1)ⁿ det S = −det S`, so
`det S = 0`. Spectrally: the spectrum is closed under `λ ↦ −λ`, and an involution on an
**odd**-cardinality set must have a **fixed point** — the eigenvalue `0`.

**That is verbatim the argument of THM-1445-A step (3)**, where the fibre of the JC
counterexample has odd size (3), so `σ` must fix a point, and cannot fix two. Same lemma,
two subjects:

- *skew matrices*: odd dimension ⟹ at least one zero eigenvalue ⟹ `P` has the factor `x`;
- *JC fibres*: odd degree ⟹ at least one `σ`-fixed preimage (THM-1350) ⟹ the escaping pair
  must be the `σ`-swapped pair ⟹ "reflection = torus".

The factor `x` in `x(x²+7)(x⁴+14x²+17)` and the `σ`-fixed sheet in the JC fibre are **the
same fixed point of the same involution**, appearing in two different categories.

## 6. Open

1. **Why eleven?** The collapse from 32768 switching classes to 11 spectra is not explained
   here. `(c₃, c₁)` are sums of `4×4` and `6×6` principal Pfaffians squared; the small image
   is presumably a strong congruence. Worth a look — it is a cheap, sharp question.
2. **Is `P` modal at every odd `n`?** Is there a "generic tournament spectrum" for each `n`,
   and does its share tend to a limit? At `n = 7` it is 30.8 %.
3. **Isoclinic ⟺ Paley?** Is the equal-rate (isoclinic) spectrum `x(x²+p)^{(p−1)/2}` achieved
   *only* by the Paley class at every prime `p ≡ 3 (mod 4)`? At `p = 7` it holds for 240
   classes — are those exactly the Paley switching class and its relabelings?

## Verification

`04-computation/skew_charpoly_heptad_opus_S407.py` (exhaustive census over all 32768
switching classes, exact integer Faddeev–LeVerrier),
`04-computation/heptad_polynomial_anatomy_opus_S407.py` (the `u`-form, roots, Gauss sum,
Galois groups, the refuted `sin(7θ)` comparison). Outputs in `05-knowledge/results/`.
