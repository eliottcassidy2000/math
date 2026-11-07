---
id: HYP-9500
title: "Strong equicentric tournaments exist at every even order n >= 6 (centroid = circumcentre of the signature simplex, i.e. M^2 u parallel to u). Their lambda-values (-M^2 u = lambda u) fill the whole residue class n-1 mod 4 between 3 or 1 and the largest product ab, a + b = n with a, b odd. This is proved for n = 2 (mod 4) and for every n = 0 (mod 4) that carries a skew-Hadamard matrix. Strong equicentric is much weaker than skew-Hadamard, so a construction at the remaining orders should not need the skew-Hadamard conjecture"
status: >
  OPEN (CONJECTURE), with large proved parts.
  PROVED:
  - (a) for n = 2 (mod 4), by THM-4911(3c): almost-regular equicentric tournaments with lambda = 1, which are strong;
  - (a) for n = 0 (mod 4) whenever a doubly regular tournament of order n - 1 exists, equivalently (Reid-Brown) a
    skew-Hadamard matrix of order n. Border it with a source, then switch one vertex of the doubly regular part. The
    result is skew-Hadamard, hence equicentric (THM-4911(4)), and strong (proof below). Checked by computation for
    Paley q = 7, 11, 19, 23, 31.
  OPEN: (a) at the orders n = 0 (mod 4) where skew-Hadamard matrices are not known; (b) in general.
  FINITE-EXACT: all isomorphism classes n <= 10 (THM-4911 table), with 2, 6, 128 strong equicentric classes at
  n = 6, 8, 10.
session: mac-mini-2026-11-07-sevensphere
source: 05-knowledge/results/sevensphere_anticomplementary_signature_simplex_20261107.md
scripts:
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/signature_simplex.py (+ .out)
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/switching_paley.py (+ .out)
related: [THM-4911, THM-1950, THM-447 (skew-Sylvester doubling = skew-Hadamard doubling), THM-4557, THM-211]
---

# HYP-9500 — strong equicentric tournaments at every even order

## Statement

Call a tournament *equicentric* if its signature simplex (THM-4911) has centroid equal to circumcentre, i.e.
`M²u = −λu`. By THM-4911(3):
- for odd `n` these are exactly the regular tournaments;
- the non-strong ones are exactly `R_a ⇒ R_b` (both parts regular, `λ = ab`);
- `λ ≡ n − 1 (mod 4)` for even `n`.

**Conjecture.**
- **(a)** For every even `n ≥ 6` there is a strong equicentric tournament on `n` vertices.
- **(b)** For even `n`, the set of values `λ` over equicentric tournaments of order `n` is the full progression

      Λ(n) = {λ ≡ n − 1 (mod 4) : λ_0(n) ≤ λ ≤ λ_1(n)},

  where:
  - `λ_0(n) = 1` if `n ≡ 2 (mod 4)` and `3` if `n ≡ 0 (mod 4)`;
  - `λ_1(n) = max_(a+b=n, a,b odd) ab`, which is `(n/2)²` or `(n/2)² − 1`. It is attained by `R_a ⇒ R_b`.

## What is proved

- **`n ≡ 2 (mod 4)`.** THM-4911(3c) constructs strong almost-regular equicentric tournaments (`λ = 1`): two regular
  halves of odd order `n/2` with a circulant cross pattern.
- **`n ≡ 0 (mod 4)` with a skew-Hadamard matrix of order `n`.**
  - Let `R` be doubly regular of order `m = n − 1` (Reid–Brown), and `T = R ⊕ s` with `s` a source.
  - `I + M_T` is skew-Hadamard. Switching at `w ∈ R` (reversing every arc at `w`) replaces `M` by `DMD`, with `D`
    diagonal and `D_ww = −1`, and keeps `(I + DMD)ᵀ(I + DMD) = nI`. So the new tournament `T'` is equilateral, hence
    equicentric.
  - **`T'` is strong.**
    - `s` beats every vertex of `R` except `w`, and `w → s`.
    - A vertex `y` that `w` used to beat now beats `w`.
    - A vertex `y` that used to beat `w` shares `(m−3)/4 ≥ 1` common out-neighbours `z` with `w` in `R`, and
      `y → z → w` survives the switch.
    - So every vertex reaches `w`, and `w → s → (R − w)`.
  - Computational check: Paley `q = 7, 11, 19, 23, 31` (`switching_paley.out`).

## Evidence for (b) (FINITE-EXACT, all isomorphism classes)

| `n` | equicentric classes | strong | `λ` values | predicted `Λ(n)` |
|---|---|---|---|---|
| 4 | 2 | 0 | 3 | {3} |
| 6 | 5 | 2 | 1, 5, 9 | {1, 5, 9} |
| 8 | 14 | 6 | 3, 7, 11, 15 | {3, 7, 11, 15} |
| 10 | 165 | 128 | 1, 5, 9, 13, 17, 21, 25 | {1, 5, …, 25} |

The strong classes realise:
- `{1}` at `n = 6`;
- `{3, 7, 11}` at `n = 8`;
- `{1, 5, 9, 13, 17}` at `n = 10`.

The non-strong `R_a ⇒ R_b` supply the top values `ab`.

## Why this matters (STRUCTURAL)

- Equicentric is the "first moment" shadow of skew-Hadamard. Skew-Hadamard asks for `I + MᵀM = nI`, i.e. all
  pairwise inner products of the signature vectors vanish. Equicentric asks only that the all-ones vector be an
  eigenvector of `I + MᵀM`.
- (a) at all `n ≡ 0 (mod 4)` is therefore implied by the skew-Hadamard conjecture, and is far weaker.
- A direct construction would be a small, honest step in the direction of that open problem: it would produce, at
  every such order, a strong tournament whose signature simplex has the symmetry `O_T = G_T`.

## Tests

- **Search.** Random search or hill-climbing for strong equicentric tournaments at `n = 12, 16, 20` with prescribed
  `λ`, to test (b). The orders are small enough for local search, though not for enumeration.
- **Construction.** Generalise THM-4911(3c) from two blocks with `x = ±1` to blocks with `x ∈ {±1, ±3}`. The strong
  `n = 8`, `λ = 3` class has `x = (−3, −1, −1, −1, 1, 1, 1, 3)`.
