---
id: THM-1575
title: "THE CIRCULANT TOURNAMENT'S SKEW SPECTRUM IS A TANGENT, IN CLOSED FORM — and the eigenvalue equation IS the multiple-angle formula. For the rotational tournament R_n (i→j iff j−i mod n ∈ {1,…,(n−1)/2}), in the half-dictionary coordinate x = 2λ+1: (A) spec(S) = {i·tan(kπ/n) : k = 0,…,n−1} EXACTLY, and (B) char_S(x) = the ODD PART of (1+x)^n = ((1+x)^n − (1−x)^n)/2, so its coefficients are the odd-index binomials C(n,1), C(n,3), …, C(n,n). Both PROVED (closed-form sine sum) and verified n = 3..13. (C) THE TWO ARE THE SAME STATEMENT: substituting x = i·tan θ into the odd part of (1+x)^n yields i times the NUMERATOR of tan(nθ), whose zeros are θ = kπ/n — so 'what are the eigenvalues' and 'what is tan(nθ)' are one question. (D) COROLLARY, ARITHMETIC: since tan(kπ/n) depends only on d = n/gcd(k,n), the factorisation over Q is indexed by the DIVISORS of n with the factor for d of degree φ(d) — char_S(R_9) = x(x²+3)(x⁶+33x⁴+27x²+3), degrees 1+2+6 = 1+φ(3)+φ(9); at prime n the non-trivial factor is irreducible of degree n−1. The circulant skew spectrum SEES THE DIVISOR LATTICE OF n. (E) n = 7 RESOLVED CONCRETELY: there are exactly 3 regular tournaments on 7 vertices; R_7 gives x(x⁶+21x⁴+35x²+7), Paley-7 (connection set the residues {1,2,4}, doubly regular) gives x(x²+7)³, and the third gives x(x²+7)(x⁴+14x²+1). (F) HONEST NEGATIVE on the owner's x(x²+7)(x⁴+14x²+17): NO regular 7-tournament has x⁴+14x²+17 — all three are enumerated above and the closest has constant 1, not 17. The two are not variants of each other: x⁴+14x²+17 has x² = −7±4√2 whereas ours has x² = −7±4√3"
status: >
  (A) and (B) PROVED, in full, from the closed-form sine sum: the circulant symbol is
  μ_j = 2i·Σ_{k=1..m} sin(2πjk/n) with m = (n−1)/2 (THM-1440(D)), and the Lagrange
  identity Σ_{k≤m} sin(kφ) = sin(mφ/2)·sin((m+1)φ/2)/sin(φ/2) collapses it to
  μ_j = i·[cos(πj/n) − (−1)^j]/sin(πj/n), which is −i·tan(πj/2n) for even j and
  i·cot(πj/2n) for odd j; as j runs over 0..n−1 that multiset is exactly
  {i·tan(kπ/n)}.  VERIFIED n = 3, 5, 7, 9, 11, 13 (both the polynomial identity and the
  numerical spectrum).
  (C) PROVED — it is the standard tangent multiple-angle expansion, read backwards.
  (D) PROVED, immediate from (A) plus tan(kπ/n) = tan(k′π/d) for d = n/gcd(k,n);
  verified by explicit factorisation over Q at n = 3..13.
  (E) VERIFIED-EXACT by full enumeration of the regular tournaments on 7 vertices
  (3 iso classes, brute-force canonical form over all 5040 relabellings).
  (F) VERIFIED-EXACT — a negative, and stated as one.
  NOVELTY, HONESTLY: the closed form is elementary and is very likely classical for
  rotational tournaments; it has NOT been checked against the literature.  It is
  recorded because the repo previously had only the SINE-SUM form (THM-1440(D)) and
  because it closes three separate owner threads at once.
source: kind-pasteur-2026-07-20-S128c119 (owner: even/odd vs positive/negative; and the earlier thread on sin/cos and x(x²+7)(x⁴+14x²+17))
extends:
  - THM-1440    # (D) circulant symbol as a sine sum; Re = -1/2 critical line for circulants
depends_on:
  - THM-1555    # the half dictionary: x = 2λ+1 is the coordinate in which this is clean
related: [THM-1455, THM-1475]
script: 04-computation/circulant_tangent_spectrum_kps_S128c119.py, x_half_shift_regular7_kps_S128c119.py, even_n_perron_deficiency_kps_S128c119.py (+ .out)
---

# THM-1575 — the circulant skew spectrum is a tangent

THM-1440(D) established that circulant tournaments "make the sine literal":
`μ_j = 2i·Σ_{k∈C} sin(2πjk/n)`. That is a **sine sum**. For the standard connection set
it closes to a **closed form**, and the closed form is a tangent.

Let `R_n` be the rotational tournament: `i → j` iff `j − i mod n ∈ {1, …, (n−1)/2}`.

## A, B — the two statements

Work in the half-dictionary coordinate `x = 2λ + 1` (THM-1555: this is the inverse of
the owner's `x ↦ (1+x)/2`, and it is the coordinate in which `S = 2A − J + I` becomes
simply `S`). Then

> **(A)  `spec(S) = { i·tan(kπ/n) : k = 0, …, n−1 }`  exactly.**
>
> **(B)  `char_S(x) = ((1+x)ⁿ − (1−x)ⁿ)/2`**, the **odd part of `(1+x)ⁿ`**, whose
> coefficients are the odd-index binomials `C(n,1), C(n,3), …, C(n,n)`.

Verified at `n = 3, 5, 7, 9, 11, 13`, both as a polynomial identity and numerically:

| `n` | `char_S(x)` | odd binomials |
|---|---|---|
| 3 | `x(x² + 3)` | 3, 1 |
| 5 | `x(x⁴ + 10x² + 5)` | 5, 10, 1 |
| 7 | `x(x⁶ + 21x⁴ + 35x² + 7)` | 7, 35, 21, 1 |
| 9 | `x(x² + 3)(x⁶ + 33x⁴ + 27x² + 3)` | 9, 84, 126, 36, 1 |
| 11 | `x(x¹⁰ + 55x⁸ + 330x⁶ + 462x⁴ + 165x² + 11)` | 11, 165, 462, 330, 55, 1 |

*Proof.* The circulant symbol at `ω^j` is `μ_j = Σ_{k=1}^{m}(ω^{jk} − ω^{−jk})
= 2i·Σ_{k=1}^{m} sin(2πjk/n)` with `m = (n−1)/2`. The Lagrange identity
`Σ_{k≤m} sin(kφ) = sin(mφ/2)·sin((m+1)φ/2)/sin(φ/2)` at `φ = 2πj/n` gives, via
`sin a·sin b = [cos(a−b) − cos(a+b)]/2`,

> `μ_j = i·[cos(πj/n) − (−1)^j] / sin(πj/n)`.

For `j` even this is `−i·tan(πj/2n)`; for `j` odd it is `i·cot(πj/2n) =
i·tan(π(n−j)/2n)`. Writing `j = 2k` in the first case and `n − j = 2k′` in the second
(possible since `n` is odd), both are of the form `i·tan(kπ/n)`, and as `j` runs over
`0..n−1` the multiset obtained is exactly `{i·tan(kπ/n)}`. (B) follows since a monic
polynomial is determined by its roots. ∎

## C — the two statements are one, and it is the multiple-angle formula

Substitute `x = i·tan θ` into the odd part of `(1+x)ⁿ`. Since
`(1 ± i tan θ)ⁿ = (cos θ ± i sin θ)ⁿ / cosⁿθ = e^{±inθ}/cosⁿθ`, the odd part is
`i·sin(nθ)/cosⁿθ` — i.e. **`i` times the numerator of `tan(nθ)`** in its expansion in
powers of `tan θ`. Its zeros are `θ = kπ/n`. So

> *"what are the eigenvalues of `R_n`"* and *"what is the multiple-angle formula for the
> tangent"* are the same question.

This is the concrete form of the owner's sin/cos thread: THM-1440(A) says `char_S` is an
odd function at odd `n` "like sin"; here the odd function is *literally* `sin(nθ)`,
divided by `cosⁿθ`.

## D — corollary: the spectrum sees the divisor lattice of `n`

`tan(kπ/n)` depends only on `d = n/gcd(k,n)`. So the roots partition by divisors of `n`,
and the factorisation of `char_S` over `ℚ` is **indexed by the divisors**, the factor for
`d` having degree `φ(d)`:

> `char_S(R_n) = x · ∏_{d | n, d > 1} Ψ_d(x)`,  `deg Ψ_d = φ(d)`.

At `n = 9`: `1 + φ(3) + φ(9) = 1 + 2 + 6 = 9` ✓ — and `(x² + 3)` is exactly the `d = 3`
factor, appearing because `3 | 9`. At prime `n` there is a single non-trivial factor,
irreducible of degree `n − 1`. The constant term is `C(n,1) = n`, which is the classical
`∏_{k=1}^{n−1} tan(kπ/n) = ±n`.

This is a **tangent analogue of the cyclotomic factorisation**, and it means the
circulant tournament's skew spectrum carries arithmetic information about `n` — the same
kind of `d | n` structure the repo meets in the mod-`p` and Paley threads.

## E — the three regular tournaments on 7 vertices, resolved

There are exactly **3** regular tournaments on 7 vertices up to isomorphism (full
enumeration over all `2²¹` tournaments with canonical forms over all `5040`
relabellings):

| class | `char_S(x)` | which tournament |
|---|---|---|
| [1] | `x(x⁶ + 21x⁴ + 35x² + 7)` | `R_7`, rotational, connection set `{1,2,3}` |
| [2] | `x(x² + 7)(x⁴ + 14x² + 1)` | the third regular class |
| [3] | `x(x² + 7)³` | **Paley-7**, connection set the residues `{1,2,4}`, doubly regular |

Paley is the one with a *repeated* factor — doubly regular means exactly two distinct
non-Perron eigenvalues, here `±i√7`, matching THM-1440(D)'s `{0} ∪ {±i√q}`.

## F — honest negative on `x(x²+7)(x⁴+14x²+17)`

The owner's polynomial from the earlier sin/cos thread was `x(x²+7)(x⁴+14x²+17)`. The
`x` and the `(x²+7)` are exactly right and the `14` is exactly right — but **no regular
tournament on 7 vertices has `x⁴ + 14x² + 17` as a skew factor.** All three are
enumerated in E and the closest is class [2] with constant **1**, not 17.

They are not variants of one another: `x⁴+14x²+17` has `x² = −7 ± 4√2`, while
`x⁴+14x²+1` has `x² = −7 ± 4√3`. A `√2` where the tournament gives `√3` points at a
different structure — order-8 or a conference/Seidel object rather than a 7-vertex
tournament skew matrix. **Not chased here, and flagged as unresolved rather than
massaged into agreement.**

## Named next

- Identify the `√2` object. `−7 ± 4√2` and the shape `x⁴ + 14x² + 17` are specific
  enough to be searched for directly among Seidel matrices of graphs on 8 vertices and
  among conference matrices, which is where `√2` naturally appears.
- Prove `deg Ψ_d = φ(d)` implies `Ψ_d` is irreducible over `ℚ` (it is the minimal
  polynomial of `i·tan(π/d)` iff `[ℚ(i tan(π/d)) : ℚ] = φ(d)`; true at every `n` tested,
  and it should follow from the cyclotomic field degree).
- The even-`n` companion: no regular tournament exists, and by THM-1555 I.c the Perron
  deficiency is the spectral height above the line. The optima are algebraic of full
  degree — `n = 4` gives `λ⁴ − 2λ − 1` (irreducible, `ρ = 1.39534`), `n = 6` gives
  `λ⁶ − 8λ³ − 12λ² − 8λ − 2` (irreducible, `ρ = 2.43397`). So **the extremal even-`n`
  tournament has no rational eigenvalue at all**, while the extremal odd-`n` one has
  `ρ = (n−1)/2` rational. That parity split is worth its own statement.
