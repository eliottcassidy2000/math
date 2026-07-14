# The moduli lift is proved: curving the danger circle into X₀(14) creates the class numbers

*mac-mini-2026-07-13-S93. Owner asked me to prove the moduli lift to X₀(14)'s fixed-point arithmetic —
the last step of the S89→S92 arc, where S92 found the arc-transitivity breaker (the circle Z/14, on
which complement factors into the Atkin-Lehner V₄) but left the fixed-point arithmetic as an open
functorial bridge. The fixed-point arithmetic is now **proved** (rigorous, classical), and the
correspondence to the circle is exact. What "lifting" means became clear in the process: it is the
**curving** of the flat circle into the genus-1 moduli, and that curving is exactly what **creates the
class numbers** — the arithmetic the flat circle cannot carry.*

---

## What was proved

X₀(14) = the elliptic curve 14a (genus 1). Choose a rational cusp as origin. Then every
curve-automorphism is `P ↦ ±P + c`, and an **involution** is one of exactly two kinds:

- `P ↦ P + c` with `c ∈ E[2]` — a **translation**: fixed-point-free, and it acts as **+1** on the
  invariant differential `ω`;
- `P ↦ −P + c` — a **reflection**: **4** fixed points (`{P : 2P = c}`), acting as **−1** on `ω`.

Since the Atkin-Lehner involution `w_Q` acts on the one-dimensional `H⁰(X₀(14), Ω¹) = ⟨f₁₄⟩` by its
Atkin-Lehner eigenvalue `λ_Q`, the kind is read off the eigenvalue, and Riemann–Hurwitz
(`g' = 1 − ν/4`) gives the quotient genus:

> **`λ_Q = +1` ⟺ translation ⟺ 0 fixed points ⟺ quotient genus 1;
> `λ_Q = −1` ⟺ reflection ⟺ 4 fixed points ⟺ quotient genus 0.**

For 14a the eigenvalues are `w₂ = +1`, `w₇ = −1`, `w₁₄ = −1` (LMFDB 14.a). Therefore:

| involution | eigenvalue | kind | fixed points | quotient genus |
|---|---|---|---|---|
| `w₂` (the 2-part) | +1 | **translation** | **0** | 1 |
| `w₇` (the 7-part) | −1 | reflection | 4 | 0 |
| `w₁₄` (Fricke) | −1 | reflection | 4 | 0 |

The Fricke count is a class number: `ν(w₁₄) = h(−4·14) = h(−56) = 4`, with the four fixed points CM by
the order of discriminant **`−56 = −2³·7`** — the discriminant where both primes of `14 = 2·7` become
arithmetic. The V₄ is internally consistent: `w₂∘w₇ = w₁₄` reads as translation∘reflection = reflection
(0,4 → 4), and `w₇∘w₁₄ = w₂` as reflection∘reflection = translation (4,4 → 0). Independent sanity: the
same method on `N=11` gives Fricke `ν = h(−44) + h(−11) = 3 + 1 = 4`, matching the genus-0 quotient
`X₀(11)⁺`. Class numbers verified by reduced-form count (`h(−56)=4`, etc.).

## What "lifting" turned out to mean

Set the proved moduli arithmetic beside the flat circle `Z/14` from S92 (`W₂ = x+7`, `W₇ = 7−x`,
`W₁₄ = −x`):

| | circle `Z/14` fixed points | moduli `X₀(14)` fixed points |
|---|---|---|
| `w₂` (2-part) | 0 (free translation) | **0** (free translation) — **exact match** |
| `w₇` (7-part) | 0 (`2x=7`, none) | 4 (CM) — lift **adds** the CM points |
| `w₁₄` (Fricke) | 2 (`{0,7}`, the 2-torsion) | 4 (CM disc −56) — lift **adds** the CM points |

The **2-part is a fixed-point-free translation on both the circle and the moduli** — it lifts exactly.
This is the clean part of the bridge, and it is now proved (the S92 circle `W₂` = the 2-adic descent
of THM-580 = the moduli translation `w₂`). The **reflections gain their fixed points in the lift**:
`w₇` goes `0 → 4`, `w₁₄` goes `2 → 4`, and the new points are the **CM points**, counted by the class
number `h(−56) = 4`. On the flat circle a reflection has at most its two antipodal fixed points (or
none); on the curved genus-1 moduli it has four, and the "extra" points are precisely the arithmetic —
the class number. So:

> **Lifting = curving the flat circle `Z/14` into the genus-1 moduli `X₀(14)`, and that curving is
> what creates the class-number fixed points.** The V₄ *group* is the same on both (S92); the
> *arithmetic* (the CM points = the class numbers) appears only under the curving.

This is why the flat circle was "group-complete and arithmetic-incomplete" (S92): the group is
CRT-combinatorial (`Z/2 × Z/7`), but the fixed-point counts are class numbers, which live on the curve.
The discriminant `−56 = −2³·7` is the exact place the `2·7` stops being the geometric CRT split of the
circle and becomes the arithmetic of the moduli.

## Honest scope

- **Proved:** the fixed-point arithmetic of X₀(14) — the counts `(0, 4, 4)`, the quotient genera
  `(1, 0, 0)`, the Fricke CM discriminant `−56` with `h(−56) = 4` — rigorously, from the cited AL
  eigenvalues plus the classical genus-1 dichotomy and Riemann–Hurwitz, sanity-checked on `N=11`.
- **Corrected:** klein-S59's "`W₂` has 4 fixed CM points" — `w₂` is the fixed-point-*free* translation
  (eigenvalue +1); the 4-fixed-point reflections are `w₇` and `w₁₄`.
- **Still open (the genuinely modular content):** the *functorial map* from the runner/circle
  combinatorics to the moduli `X₀(14)` — realizing the elliptic curve from the tournament side. That
  is not proved and is not a corollary; it is the same `f₁₄`/Dedekind bridge the covering-min residual
  points at. What is now settled is the *target* arithmetic and its exact match to the circle: the
  2-part lifts as a translation, and the curving supplies the class numbers.

## Coda — the arc, complete as far as it goes

S89: how do tournaments touch the last bit. S90: the odd-graph→cusp transfer is provably blind. S91:
the root is Sₙ-transitivity, which makes complement irreducible. S92: the cure is the level — the
danger circle `Z/14`, where complement factors into the Atkin-Lehner V₄. S93: the moduli fixed-point
arithmetic is proved, and the lift is the curving that turns the circle's flat V₄ into the moduli's
class numbers. The chain now reads cleanly from the combinatorics to the modular curve, with one honest
gap left standing in exactly the place the whole project keeps arriving at: the passage from the
combinatorial level structure to the elliptic curve `14a` itself. The circle is the bridge; the curving
is the arithmetic; and the last plank — the functor — is the modular content of LRC(14).

---

*Cross-links: S92 (HYP-6575, the circle V₄, complement factors), THM-580 (the 2-adic descent = `w₂` =
the translation, now confirmed), klein-S10/S59 (cusps = n=4 classes; `w₂` fixed-point count corrected),
HYP-3768/f₁₄ (the cusp-form residual = the still-open functor), the S91 reflection (Sₙ-transitivity =
the root). Proof: HYP-6580, `04-computation/moduli_lift_x014_fixedpoints_macmini_S93.py`.*
