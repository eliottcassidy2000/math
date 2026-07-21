---
id: THM-1800
title: "BINARY FORMS <-> TOURNAMENTS: (in)transitivity IS representation theory. The dictionary: n players = n roots of a binary form of degree n on P^1; a tournament = an orientation of every edge; TRANSITIVITY (a linear order) = the SIGN character of S_n = the Vandermonde prod(x_i-x_j) (the odd/skew part, THM-1450), while the DISCRIMINANT prod(x_i-x_j)^2 is the even/S_n-invariant part, orientation-blind; INTRANSITIVITY (3-cycles) = the first SL(2)-invariant BEYOND the sign (the cubic covariant/Hessian), whose vanishing = coincident roots = the ramification R of THM-1770. THE CHARACTER/DISCRIMINANT TOURNAMENT: on F_p (p = 3 mod 4), i->j iff chi(j-i)=+1, chi = the quadratic (Legendre) character = the discriminant character of the binary quadratic x^2-a; this is the Paley tournament, and it is MAXIMALLY INTRANSITIVE -- doubly-regular, #3-cycles = (p+1)p(p-1)/24 (verified p=3,7,11,19,23), a Jacobi-sum invariant, exceeding the random count p(p-1)(p-2)/24 by the factor (p+1)/(p-2). THE Sym^3 END: the cyclic 3-tournament <-> discriminant SQUARE <-> A_3 Galois; the transitive <-> non-square <-> S_3; the S_3 resolvent of a binary cubic IS the triangle's transitive/intransitive dichotomy (the repo's 'Redei sign = discriminant character', 'generic fibre = cyclic 3-tournament', THM-1375 IV, THM-1770)"
status: SYNTHESIS + VERIFIED (the Paley cyclic-triple count (p+1)p(p-1)/24 exact for p=3,7,11,19,23; the Vandermonde=sign, discriminant=even identities). The dictionary is exact; the 'intransitivity = cubic covariant' SL(2) statement is a precise classical identification.
author: opus-2026-07-20-S432
depends_on: [THM-1450 (odd=skew=tournament, even=symmetric), THM-1770 (Sym^3 counterexample; ramification R), THM-1375 (Redei sign = discriminant character; cyclic-3 fibre), THM-1200 (the two sevens; Paley)]
---

# THM-1800 — Binary forms ↔ tournaments: (in)transitivity is representation theory

Owner: work the representation theory of binary forms and how it relates to tournaments — which
are in/transitivity itself. Here is the dictionary, a concrete invariant, and the tie to the
repo's Paley / Sym³ / Rédei-sign threads.

## 1. The dictionary

| binary form / `SL(2)` object | tournament object |
|---|---|
| `n` roots of a degree-`n` form on `ℙ¹` | `n` players / vertices |
| a chosen **orientation** of each edge | a tournament ∈ `(ℤ/2)^{C(n,2)}` |
| **Vandermonde** `∏_{i<j}(x_i−x_j)` = the **sign** character of `S_n` | the **transitive** (linear-order) orientation and its parity |
| **discriminant** `∏(x_i−x_j)²` (`S_n`-invariant, even) | orientation-**blind** data (which pairs are edges) |
| the **cubic covariant / Hessian** (first invariant beyond sign) | **intransitivity** — the `3`-cycle content |
| **coincident roots** `x_i = x_j` (discriminant `= 0`) | the **ramification** `R` (THM-1770) |

The odd/even split is THM-1450's **skew/symmetric** split: the tournament (orientation) lives in
the **odd** part (Vandermonde flips under a transposition = one edge reverses); the
discriminant is the **even** part.

## 2. The character tournament is the discriminant construction (verified)

On `F_p` (`p ≡ 3 mod 4`, so `χ(−1) = −1` gives antisymmetry), define `i → j ⟺ χ(j−i) = +1`
where `χ` is the **quadratic (Legendre) character**. `χ(a)` is exactly the **discriminant
character** of the binary quadratic `x² − a` (`a` is a QR iff `x²−a` splits over `F_p`). This is
the **Paley tournament**.

> **The Paley tournament is maximally intransitive.** It is doubly regular, with
> ```
> #(3-cycles) = (p+1)p(p−1)/24 ,
> ```
> verified `p = 3,7,11,19,23` (values `1, 14, 55, 285, 506`). This is a **Jacobi-sum invariant**
> of `χ` (a cubic character sum), and it **exceeds** the random-tournament count
> `p(p−1)(p−2)/24` by the factor `(p+1)/(p−2) → 1`. Doubly-regular = the maximum of the 3-cycle
> count, so **the discriminant/character orientation is the *most intransitive* tournament** —
> intransitivity is not noise, it is what the quadratic character *maximises*.

## 3. The Sym³ end: the resolvent is the triangle's dichotomy

A binary **cubic** has 3 roots; `disc = ((a−b)(b−c)(c−a))²`, and its square root
`(a−b)(b−c)(c−a)` is the `S_3`-**sign** = the cyclic orientation `a→b→c→a`. Hence:

- `disc` a **square** (√ rational) `⟺` Galois `= A_3` (cyclic) `⟺` the **cyclic 3-tournament**
  (a coherent 3-cycle = intransitive);
- `disc` a **non-square** `⟺` `S_3` `⟺` the **transitive** (ordered) triple.

This is precisely the repo's **"Rédei sign = discriminant character"** and **"generic fibre =
cyclic 3-tournament"** (THM-1375 IV, THM-1770): **the `S_3` resolvent of a binary cubic *is* the
triangle's transitive/intransitive dichotomy.** The Sym³ counterexample (THM-1770) is this at
the level of the whole space `Sym³(ℙ¹) = ℙ³`, small diagonal = twisted cubic.

## 4. The `SL(2)` reading

`V_n = Sym^n(ℂ²)` is the `(n+1)`-dimensional irreducible `SL(2)`-representation (binary forms of
degree `n`). The tournament orientation is a **sign-twisted** datum: it transforms by the sign of
`S_n ⊂ ` the monodromy, i.e. it is a section of the sign-twisted line bundle whose square is the
discriminant. **Intransitivity is the first `SL(2)`-covariant beyond the sign** — for cubics the
**Hessian** (a binary quadratic covariant), whose discriminant is the cubic's discriminant.
Coincident roots (Hessian/discriminant vanishing) are exactly the **ramification** where the
tournament degenerates (an edge becomes undirected). So:

> **(In)transitivity is the representation theory of the sign-twisted line on the configuration
> space of `n` points**, with the discriminant as its square and the covariants (Hessian,
> Jacobian, …) as the higher-order intransitivity data.

## 5. Why this matters for the repo

- **Paley's role is explained invariant-theoretically:** it is the *discriminant/quadratic-
  character* orientation, hence maximally intransitive (doubly regular) — which is why it recurs
  as the extremal object (THM-1200 the two sevens, THM-1450 the isoclinic spectrum, the LRC
  heptagon). Its 3-cycle count is a Jacobi sum.
- **The odd=skew=tournament / even=symmetric split (THM-1450)** is the sign/discriminant split
  here — one dictionary entry.
- **The JC fibre (cyclic-3)** and the **Sym³ counterexample (THM-1770)** are the `n=3` end of the
  same dictionary; the tangent-not-osculating hyperplane there is the discriminant-tangency that
  distinguishes `A_3` (cyclic) from `S_3` (transitive).

## 6. Open / tangential

1. **Higher intransitivity invariants.** For `n ≥ 4`, the `SL(2)`-covariants of binary quartics
   (the invariants `I`, `J`, the sextic covariant) should index the higher tournament statistics
   (4-cycles, the score sequence's symmetric functions). Which covariant = which cycle count?
2. **The Jacobi-sum count and `H`.** The tournament determinant `H` (the project's central
   invariant) of the Paley class is a character-sum; expressing `H_Paley(p)` as an explicit
   product of Jacobi/Gauss sums would tie `H` to binary-form invariant theory directly.
3. **The maximally-intransitive extremal.** Doubly-regular = discriminant-character = maximal
   3-cycles; is *every* extremal-`H` tournament a character construction (a binary-form
   discriminant orientation)? The census (Paley maximal at several `p`) suggests a
   representation-theoretic characterisation of the `H`-extremiser.

## Verification

`04-computation/binary_forms_tournaments_opus_S432.py` — the Paley cyclic-triple count
`(p+1)p(p−1)/24` for `p=3,7,11,19,23`; the Vandermonde = sign, discriminant = `S_n`-invariant
identities; the cubic discriminant-square = cyclic-3. Output in `05-knowledge/results/`.
