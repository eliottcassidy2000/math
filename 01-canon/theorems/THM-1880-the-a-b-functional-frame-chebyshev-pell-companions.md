---
id: THM-1880
title: "THE a/b FUNCTIONAL FRAME: a(x)=x+1 AND b(x)=x/2 GENERATE THE TRANSITIVE TOURNAMENT'S CHEBYSHEV–PELL COMPANION POLYNOMIALS, AND b∘a IS THE HALF-DICTIONARY. Writing a(x)=x+1, ā(x)=x−1 (the conjugate shift), b(x)=x/2 (the symmetriser), the transitive tournament's two invariant polynomials and their companion are all a/b-composites: char_A = xⁿ (the GIT nullcone monomial, THM-1810); char_S = E_n := b(aⁿ+āⁿ) = ((x+1)ⁿ+(x−1)ⁿ)/2 (the even/skew form, THM-1875); and the ODD companion O_n := b(aⁿ−āⁿ) = ((x+1)ⁿ−(x−1)ⁿ)/2. (1) THE PELL/CHEBYSHEV IDENTITY, verified exactly n=1..8: E_n² − O_n² = (x²−1)ⁿ, so (E_n,O_n) is the (cosh,sinh) pair of the metric x²−1 — E_n and O_n are Chebyshev polynomials of imaginary/cotangent argument. (2) THE CORRECTED COUPLED cos/sin RECURSION: E_n = x·E_{n−1}+O_{n−1}, O_n = E_{n−1}+x·O_{n−1}, with E_n+O_n = aⁿ, E_n−O_n = āⁿ — i.e. a and b build the pair by one multiplication + one symmetrisation per step. (3) THE SPECTRA are pure cotangent (roots of unity of the odd/even circle): E_n has roots i·cot((2k−1)π/2n), O_n has roots i·cot(kπ/n). (4) b∘a = (x+1)/2 IS THE HALF-DICTIONARY (THM-1555, the {−1,0,1}↔{0,½,1} conjugation), its inverse a⁻¹∘b⁻¹ = 2x−1; so the ½ that runs through the corpus (fiber fraction, Legendre exponent, Re=−½ line, half-angle) is literally the b in this monoid. (5) TRIANGULAR NUMBERS ARE THE COEFFICIENTS: E_n = Σⱼ C(n,2j) x^{n−2j}, so its subleading coefficient is C(n,2) = T_{n−1} = the number of ARCS of the tournament — the triangular number is the second symmetric function of the skew spectrum, and E_n/O_n carry the even/odd binomials C(n,2j)/C(n,2j+1). So 'triangular numbers, and thus tournaments, composed of a and b recursively' is exact: the tournament's skew characteristic form is a(x) and b(x) folded n times, its arc-count triangular number is a coefficient, and its spectrum is trigonometric (cot) — the odd/even (sin/cos/tan) = SL₂ Weyl axis made a two-generator functional monoid"
status: >
  (1),(2),(3) VERIFIED-EXACT in sympy (E_n²−O_n²=(x²−1)ⁿ, the corrected coupled recursion residual = 0,
  and the cot spectra) for n = 1..8; all are classical identities for the Chebyshev/Fibonacci-Pell
  polynomials of the pair (x+1, x−1) and are provable in a line from E_n±O_n = (x±1)ⁿ.
  (4) THM-1555 (b∘a = the half-dictionary) restated; exact.
  (5) The coefficient C(n,2) = T_{n−1} = #arcs is e₂ of the skew spectrum = Σ_{i<j} S_{ij}² =
  C(n,2), exact.
  This is a UNIFYING FRAME (a two-generator functional monoid tying THM-1555/1810/1875 to
  Chebyshev/Pell/trig), not a new open-problem advance; every equation is verified/classical.
source: kind-pasteur-2026-07-21-S128c138 (owner: think of triangular numbers/tournaments as a(x)=x+1 and b(x)=x/2 composed recursively; think functionally, trigonometrically)
depends_on:
  - THM-1875    # char_S(transitive) = ((x+1)^n+(x-1)^n)/2, the even form
  - THM-1555    # the half-dictionary = b∘a
  - THM-1810    # transitive = GIT nullcone char_A = x^n
related: [THM-1725, THM-293]
external:
  - "Chebyshev / Fibonacci–Pell polynomials; the staircase-skew cotangent spectrum; OEIS A000182 (tangent) via the succession GF W (THM-293)."
script: 04-computation/deep_archaeology_kps_S128c137.py (+ /tmp ab check)
---

# THM-1880 — the a/b functional frame

The owner's frame: build tournaments from `a(x) = x+1` and `b(x) = x/2`, composed recursively;
think functionally and trigonometrically. Here is the exact structure it names.

## The three polynomials are all a/b-composites

With `a(x)=x+1`, `ā(x)=x−1` (conjugate shift), `b(x)=x/2` (symmetriser), the transitive
tournament carries:

| object | polynomial | a/b form |
|---|---|---|
| adjacency `char_A` | `xⁿ` | the nullcone monomial (THM-1810) |
| skew `char_S` = `E_n` | `((x+1)ⁿ+(x−1)ⁿ)/2` | `b(aⁿ + āⁿ)` — even (THM-1875) |
| odd companion `O_n` | `((x+1)ⁿ−(x−1)ⁿ)/2` | `b(aⁿ − āⁿ)` |

## Chebyshev–Pell

Verified `n = 1..8`:

> **`E_n² − O_n² = (x²−1)ⁿ`** (Pell), and **`E_n = x·E_{n−1} + O_{n−1}`, `O_n = E_{n−1} + x·O_{n−1}`**
> (the cos/sin recursion), with `E_n + O_n = aⁿ`, `E_n − O_n = āⁿ`.

The first version transposed the `E/O` labels in this recurrence; it already
failed at `n=2`.  THM-2142 recorded the crossed formula, and MISTAKE-367 records
the live-canon repair.  The defining closed forms, Pell identity, spectra, and
coefficient statements were unaffected.

So `(E_n, O_n)` is the `(cosh, sinh)` pair of the quadratic metric `x²−1`: they are Chebyshev
polynomials of cotangent argument. The spectra are pure cotangent — `E_n` at
`i·cot((2k−1)π/2n)`, `O_n` at `i·cot(kπ/n)` — roots of unity of the odd/even circle.

## `b∘a` is the half-dictionary, and `b` is the ½

`b∘a = (x+1)/2` is exactly the `{−1,0,1} ↔ {0,½,1}` half-dictionary (THM-1555), inverse `2x−1`.
The `½` that recurs across the corpus — the tiling **fiber fraction** `(1/2)_{n−2}/(n−2)!`, the
**Legendre** generating exponent `(…)^{−1/2}`, the regular-tournament **`Re = −½`** line, the
**half-angle** of the cotangent spectrum — is one object: the generator `b` of this monoid. `a`
shifts, `b` halves; `b∘a` is the coordinate change between the sign world and the affine world.

## Triangular numbers are the coefficients

`E_n = Σⱼ C(n, 2j) x^{n−2j}`. Its **subleading coefficient is `C(n,2) = T_{n−1}` = the number of
arcs** of the tournament — the second elementary symmetric function of the skew spectrum
(`Σ_{i<j} S_{ij}² = C(n,2)`). `E_n` and `O_n` carry the even/odd binomial coefficients
`C(n,2j)` / `C(n,2j+1)`. So *"triangular numbers, and thus tournaments, composed of `a` and `b`
recursively"* is literal: the skew characteristic form is `a` and `b` folded `n` times, the
arc-count triangular number is one of its coefficients, and the spectrum is trigonometric.

## The unification

The odd/even (`sin/cos/tan`) axis — the `SL₂` Weyl involution `x ↦ −x` of the characteristic
binary form (THM-1810) — is, at the transitive vertex, a **two-generator functional monoid
`⟨a, b⟩`**. `E_n` (even, secant/cotangent side) is `char_S`; `O_n` (odd, tangent side) is its
companion; `b∘a` is the half-dictionary; the tangent numbers `A000182` sit on the odd side via the
succession GF `W` (THM-293/1875 named-next). Everything the owner has pulled — roots of unity, the
`½`, odd/even, triangular numbers — is this monoid acting on `xⁿ`.

## The two poles: spread vs concentrated (S128c139)

`char_S` is always the even form `∏ₖ (x² + λₖ²)` (skew spectrum `±iλₖ`). The `a/b` frame gives one
pole; the Gauss sum gives the other, and they are the two ends of GIT stability (THM-1810):

- **Transitive (GIT-unstable vertex).** `E_n = b(aⁿ+āⁿ)`, `λₖ = cot((2k−1)π/2n)` — the cotangent
  ladder, **maximally spread**: `var(λₖ²) = 74.7, 352, 2027` at `n = 7,11,19`.
- **Paley/regular (GIT-polystable).** `char_S(Paley_p) = x·(x²+p)^{(p−1)/2}` — **every** `λₖ² = p`
  (the imaginary Gauss sum `±i√p`, `|g(χ)|² = p`), so `var(λₖ²) = 0`. Verified `p = 7,11,19`
  (`x(x²+7)³, x(x²+11)⁵, x(x²+19)⁹`).

So the **spectral spread `var(λₖ²)` is a one-scalar GIT-instability measure**: transitive maximises
it (var `74.7` at `n=7`), Paley zeroes it (`0`), random tournaments sit between (`21–43`). The `a/b`
Chebyshev structure is the *spread* extreme; the concentrated extreme is not `b` of `aⁿ` but `b` of
a **constant** shift — `(x²+p)^m` — the Gauss-sum degeneration where the cotangent ladder collapses
to a single value. The `SL₂` Weyl axis (odd/even) carries GIT stability along it, from the monoid
`⟨a,b⟩` (spread, unstable) to the character sum (concentrated, stable).

## Named next

- **The tangent identity via `W` (THM-293).** `Σ O_n(0) t^n/n! = sinh t`, not `tan t`; the tangent
  numbers `A000182` come from the succession GF `W`, a *different* `a/b`-composite — recover it.
- **`O_n` as a tournament invariant.** `E_n = char_S`; does `O_n` (roots `i·cot(kπ/n)`) count a
  Pfaffian minor or a Seidel object?
- **The spread scalar across the census.** Is `var(λₖ²)` monotone with the `c₃`/intransitivity
  order, or does it separate a distinct stratum? (It is `0` exactly on the doubly-regular /
  Gauss-sum tournaments.)
- **The EGF / tangent identity.** `Σ O_n(0) t^n/n! = sinh t`, not `tan t`; the tangent numbers come
  from `W` (THM-293), a *different* a/b-composite — recover it and see which composite gives `tan`.
- **`O_n` as a tournament invariant.** `E_n = char_S`; does `O_n` (roots `i·cot(kπ/n)`) count a
  tournament quantity (a Pfaffian minor, a Seidel object)?
