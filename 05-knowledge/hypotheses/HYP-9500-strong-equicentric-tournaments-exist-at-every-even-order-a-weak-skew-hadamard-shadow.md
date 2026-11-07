---
id: HYP-9500
title: "The lambda-spectrum of equicentric tournaments. For even n, the values lambda with -M^2 u = lambda u over equicentric tournaments of order n (centroid = circumcentre of the signature simplex, THM-4911) fill the whole progression {lambda = n - 1 (mod 4): lambda_0(n) <= lambda <= lambda_1(n)}, with lambda_0 = 1 or 3 and lambda_1 = max ab (a + b = n, a, b odd). The companion existence statement, a strong equicentric tournament at every even n >= 6, was posed as part (a) and is PROVED (audit J1: substitution into a regular tournament, and skew-Walsh towers at powers of 2)"
status: >
  OPEN (CONJECTURE (b)).
  PROVED:
  - (a) at every even n >= 6 (audit J1, construction below). The first version mis-typed (a) as open and as a "shadow"
    of the skew-Hadamard conjecture (MISTAKE-732). It is elementary and carries no information about that conjecture.
  FINITE-EXACT: (b) for n = 4, 6, 8, 10, over all isomorphism classes (THM-4911 table).
session: mac-mini-2026-11-07-sevensphere
source: 05-knowledge/results/sevensphere_anticomplementary_signature_simplex_20261107.md
scripts:
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/signature_simplex.py (+ .out)
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/switching_paley.py (+ .out)
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/audit_J1/j1_hyp9500.py (+ .out): the substitution lemma and the construction at every even n <= 96
related: [THM-4911, THM-1950, THM-447 (skew-Sylvester doubling and the skew-Walsh tower), THM-472, Nozaki-Suda arXiv:1202.5374]
---

# HYP-9500 — the λ-spectrum of equicentric tournaments

## Setting

A tournament is *equicentric* if its signature simplex (THM-4911) has centroid equal to circumcentre, i.e.
`M²u = −λu`. Equivalently, every vertex is an endpoint of the same number `Δ` of directed 2-paths, with
`λ = (n−1)² − 2Δ`. By THM-4911(3):
- for odd `n`, equicentric means regular (`λ = 0`);
- the non-strong ones are exactly `R_a ⇒ R_b`, with `λ = ab`;
- for even `n`, `λ ≡ n − 1 (mod 4)`.

## Conjecture (b)

For even `n`, the set of `λ` over equicentric tournaments of order `n` is the full progression

    Λ(n) = {λ ≡ n − 1 (mod 4) : λ_0(n) ≤ λ ≤ λ_1(n)},

where:
- `λ_0(n) = 1` if `n ≡ 2 (mod 4)` and `3` if `n ≡ 0 (mod 4)`;
- `λ_1(n) = max_(a+b=n, a,b odd) ab`, which is `(n/2)²` or `(n/2)² − 1`, attained by `R_a ⇒ R_b`.

**Evidence** (FINITE-EXACT, all isomorphism classes).

| `n` | equicentric classes | strong | `λ` values | predicted `Λ(n)` |
|---|---|---|---|---|
| 4 | 2 | 0 | 3 | {3} |
| 6 | 5 | 2 | 1, 5, 9 | {1, 5, 9} |
| 8 | 14 | 6 | 3, 7, 11, 15 | {3, 7, 11, 15} |
| 10 | 165 | 128 | 1, 5, 9, 13, 17, 21, 25 | {1, 5, …, 25} |

**Composition rules** (audit J1). For a regular `R` of odd order `r`:
- `λ(R[T2]) = λ(T2)`;
- `λ(T1[R]) = r²·λ(T1)`.

Together with `λ(R_a ⇒ R_b) = ab`, these give part of `Λ(n)` at every even `n`. The middle of the progression
(e.g. `λ = 5, 13, 17` at `n = 10`) comes from strong classes found by enumeration only.

## Part (a), proved (audit J1)

**Substitution lemma.** Let `R` be regular of odd order `r ≥ 3`, and `T2` equicentric of order `m` with `−M2²u = λ2 u`.
Then `R[T2]` (each vertex of `R` replaced by a copy of `T2`) is strong and equicentric, with `λ = λ2`.
- Write `M = M_R ⊗ J + I ⊗ M2`. Since `M_R u = 0`, `Mu = u ⊗ x2`.
- Then `M²u = (M_R u) ⊗ (uᵀx2)u + u ⊗ M2 x2 = −λ2 u`.
- `R[T2]` is strong because `R` is.

**Construction at every even `n ≥ 6`.**
- **`n = 2^e·r` with `r` odd, `r ≥ 3`, `e ≥ 1`.** Take `R_r[R_1 ⇒ R_(2^e − 1)]`, with `λ = 2^e − 1`. For `e = 1` this
  is `R_r[single arc]`, with `λ = 1`.
- **`n = 2^e`, `e ≥ 3`.** THM-447(2)'s skew-Walsh tower gives a skew-Hadamard matrix of order `2^e`, hence a doubly
  regular `R` of order `2^e − 1`. Border `R` with a source, then switch one vertex of `R`.
  - The result stays skew-Hadamard, hence equicentric (THM-4911(4)).
  - It is strong. A vertex `y` with `y → w` shares `(m−3)/4 ≥ 1` common out-neighbours `z` with the switched vertex
    `w`, and `y → z → w` survives the switch.
  - Checked for Paley `q = 7, 11, 19, 23, 31` (`switching_paley.out`; audit J1 adds 43, 47 and the order-15 core of
    the skew-Walsh matrix of order 16).
- Verified for every even `n` from 6 to 96 (audit J1).

## Remarks

- **The skew-Hadamard conjecture.** The orders at which skew-Hadamard matrices are unknown all have an odd factor
  (e.g. 276 = 4·69), so (a) is no evidence for that conjecture. What is left of the original motivation is (b): which
  `λ` occur, and whether the strong part of the spectrum is an interval.
- **Prior work.** Nozaki–Suda (arXiv:1202.5374) study tournaments by the main angles of `√−1·M`. In their language,
  equicentric means at most two main eigenvalues `±√λ`. Their Theorem 2.7 treats `λ = 1` at order `4k+2` under an
  extra spectral condition, and shows those tournaments are equivalent to DRT(`4k+3`).

## Tests

- Local search for strong equicentric tournaments at `n = 12, 16, 20` with each prescribed `λ` in the conjectured
  progression.
- Find constructions for the middle values: blocks with `x ∈ {±1, ±3}`. The strong `n = 8`, `λ = 3` class has
  `x = (−3, −1, −1, −1, 1, 1, 1, 3)`.
