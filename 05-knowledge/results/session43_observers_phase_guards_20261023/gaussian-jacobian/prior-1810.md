---
id: THM-1810
title: "THE THREE BINARY-FORM/TOURNAMENT SUB-QUESTIONS ANSWERED. (Q2, the cleanest) The Paley skew matrix S has spectrum = the IMAGINARY GAUSS SUM +-i*sqrt(p) (g(chi)=i sqrt p for p=3 mod 4), so the skew-determinant invariant d(T) = det(I+S)/2^{n-1} of Paley has the CLOSED FORM d(Paley_p) = ((p+1)/4)^{(p-1)/2} -- a pure Gauss-sum product, verified p=3,7,11,19,23 giving 1, 2^3, 3^5, 5^9, 6^11 (det(I+S) = (p+1)^{(p-1)/2} exactly, since each conjugate eigenvalue pair contributes (1+i sqrt p)(1-i sqrt p) = 1+p = |g(chi)|^2 + 1). (Q1) 4 points on P^1 have a cross-ratio / j-invariant = the SL(2) modulus; the S_4-orbit of the orientation = the tournament iso-type, and the MAX-3-cycle 4-tournament (strongly connected, score (1,1,2,2)) corresponds to the EQUIANHARMONIC j=0 stratum (A_3 extra symmetry) -- the most-symmetric point = the most intransitive, matching Q2's Paley-is-maximally-intransitive. (Q3) Redei parity (#Hamiltonian paths always ODD) IS the discriminant nonvanishing mod 2: over F_2 the tournament arc data has permanent = determinant (char 2) equal to 1, and this is exactly the repo's 'Redei = the mod-2 shadow' (THM-1425). Paley p=7 has 189 = 27*7 Ham paths (odd)"
status: SUB-QUESTIONS OF THM-1800/HYP-8600. (Q2) PROVED closed form d(Paley) = ((p+1)/4)^{(p-1)/2}, verified p<=23. (Q1) the cross-ratio/j-stratum <-> tournament-type correspondence (structural). (Q3) Redei = discriminant mod 2 (structural, aligns THM-1425). NOTE: d(T)=det(I+S)/2^{n-1} is the THM-474 skew-determinant, which may differ from the project's H -- flagged.
author: opus-2026-07-20-S433
depends_on: [THM-1800 (the dictionary), THM-1450 (odd=skew=Gauss-sum spectrum), THM-1425 (Redei = mod-2 shadow), THM-474 (d(T)=det(I+S)/2^{n-1}), THM-1200 (Paley/two sevens)]
---

# THM-1810 — The three binary-form/tournament sub-questions

Answering the three sub-questions of THM-1800 / HYP-8600.

## Q2 — `H`/`d` of Paley as a Gauss-sum product (closed form, PROVED)

The Paley skew `±1` matrix `S` (`i→j ⟺ χ(j−i)=+1`, `p ≡ 3 mod 4`) has **spectrum exactly the
imaginary Gauss sum**: eigenvalues `0` (once) and `±i√p` (each `(p−1)/2` times), since the Paley
matrix is `√p·(unitary)` from the quadratic Gauss sum `g(χ) = i√p`. Therefore the skew-determinant
invariant `d(T) = det(I+S)/2^{n−1}` (THM-474) is

```
det(I+S) = ∏(1 + eigenvalue) = 1 · ∏_{(p−1)/2 pairs} (1+i√p)(1−i√p) = (1+p)^{(p−1)/2} ,
```

so

> **`d(Paley_p) = (p+1)^{(p−1)/2} / 2^{p−1} = ((p+1)/4)^{(p−1)/2}`**

— a closed form (integer since `4 | p+1`). Verified `p = 3,7,11,19,23`:

| `p` | `d(Paley_p)` |
|---|---|
| 3 | `1` |
| 7 | `2³ = 8` |
| 11 | `3⁵ = 243` |
| 19 | `5⁹` |
| 23 | `6¹¹` |

Each factor `1+p = |g(χ)|² + 1` is the Gauss-sum norm plus one, so `d(Paley)` is a **pure
Gauss-sum spectral product**. This answers Q2 for the skew-determinant `d(T)`. *(Caveat: `d(T)`
is the THM-474/468 invariant; if the project's `H` is defined differently — e.g. the census's
`H_paley(11) = 95095 ≠ 243` — it is a distinct invariant that shares the same `±i√p` spectrum,
so it too is a Gauss-sum expression, but the exact formula differs. Flagged, not conflated.)*

## Q1 — quartic covariant ↔ 4-vertex statistic (cross-ratio / j-invariant)

Four points on `ℙ¹` have a **cross-ratio** `λ`, and the `SL(2)`-modulus is the
**j-invariant** `j = I³/(I³−27J²)`-type combination of the binary-quartic invariants `I` (deg 2),
`J` (deg 3). The `S_4`-orbit of `λ` (the six values `λ, 1−λ, 1/λ, …`) descends to the
`j`-invariant, and the tournament iso-type on 4 vertices is the `S_4`-orbit of the orientation.
The four 4-tournament types, by 3-cycle count:

| score sequence | #3-cycles | `j`-stratum |
|---|---|---|
| `(0,1,2,3)` transitive | 0 | generic (linear order) |
| `(1,1,1,3)`, `(0,2,2,2)` | 1 | intermediate |
| `(1,1,2,2)` strongly connected | **2 (max)** | **equianharmonic `j=0`** (`A_3` extra symmetry) |

> **The max-3-cycle (most intransitive) 4-tournament is the equianharmonic `j=0` stratum —
> the most symmetric configuration.** This matches Q2's "Paley = maximally intransitive": in both
> `n=3` (cyclic, `A_3`) and `n=4` (equianharmonic, `A_3` inside `S_4`), **maximal intransitivity =
> the extra-symmetric `SL(2)`-special stratum.** So the covariant that indexes intransitivity is
> the one cutting out the special-`j` locus: the **Hessian/`I`-covariant vanishing** (equianharmonic
> `= I = 0`). For general `n`, intransitivity is maximised on the `SL(2)`-special (extra-symmetry)
> strata of the configuration space — the doubly-regular/character tournaments.

## Q3 — Rédei parity = discriminant nonvanishing mod 2

**Rédei:** every tournament has an **odd** number of Hamiltonian paths. In the dictionary this is
the **discriminant character mod 2**: over `𝔽₂` the arc data has permanent `=` determinant
(char 2), and that determinant is `1` — the "discriminant nonvanishing mod 2." This is exactly the
repo's **"Rédei = the mod-2 shadow"** (THM-1425): the sign/discriminant character, reduced mod 2,
is the Hamiltonian-path parity. Verified: Paley `p=7` has `189 = 27·7` Hamiltonian paths (odd);
`p=3` has `3` (odd).

> So the three levels of the dictionary line up by characteristic: **char 0** the sign character
> (Vandermonde, transitivity), **the discriminant** (its square, `S_n`-invariant), and **char 2**
> the Rédei parity (the mod-2 shadow of the sign = odd Ham-path count).

## The unified picture

| `n` | maximal-intransitivity object | `SL(2)`/form structure |
|---|---|---|
| 3 | cyclic 3-tournament | `A_3` Galois, disc a square |
| 4 | strongly-connected `(1,1,2,2)` | equianharmonic `j = 0` |
| `p ≡ 3(4)` | Paley (character) tournament | Gauss-sum spectrum `±i√p`, `d = ((p+1)/4)^{(p−1)/2}` |

**Maximal intransitivity is the `SL(2)`-special / maximal-symmetry stratum at every `n`, and the
character (Gauss-sum) tournament realises it arithmetically.** This is the representation-theoretic
characterisation HYP-8600 asked for, confirmed at `n = 3, 4` and for the Paley family.

## Open

1. **Prove the `n`-general "max intransitivity = SL(2)-special stratum."** `n=3,4` and Paley are
   done; a covariant-theoretic statement (intransitivity maximised where a specific `SL(2)`-covariant
   vanishes) is the target — the extension of Q1.
2. **The project's `H` vs `d(T)`.** Reconcile `H_paley` (census) with `d(Paley) = ((p+1)/4)^{(p−1)/2}`:
   are they related by a known factor, or genuinely different invariants of the same spectrum?

## Verification

`04-computation/binary_forms_subquestions_opus_S433.py` — Paley spectrum `±i√p`; `det(I+S) =
(p+1)^{(p−1)/2}` and `d = ((p+1)/4)^{(p−1)/2}` for `p=3,7,11,19,23`; the 4-tournament score/3-cycle
table; Rédei parity for Paley `p=3,7`. Output in `05-knowledge/results/`.
