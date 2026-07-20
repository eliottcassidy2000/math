---
id: THM-1525
title: "THE GMC(4) COUNTEREXAMPLE: verified, improved from 6 terms to 4, and explained — and NO, we did not have it. (0) PRIORITY: the repo's backlog states the Zhao vanishing / image / Mathieu witness is 'OPEN, unclaimed' and that nobody in repo or literature has one; my own THM-1430 produced a symmetric-case Keller counterexample on ℂ⁶ that was explicitly NOT a VC witness, and our planned route was N = 158 with a 1660-monomial quartic. The outside result reaches the same target at n = 4 with a cubic. It is strictly better than anything we had. (I) VERIFIED exactly, with the ambiguous '4 real Gaussians' reading pinned: W_j = conj(Z_j) is forced, because independent W's make every Q-moment vanish identically by charge and could never give m!. (II) IMPROVED: the published 6 terms are not minimal — P₄ = (1+Z)(conj(Z) − |Z₁|²) is a FOUR-term cubic with the same E[P^m] = 0 and E[Z P^m] = m!, verified to m = 9. (III) EXPLAINED: E[P^m] = Σᵢ(−1)ⁱC(m,i)²(m−i)!μᵢ with μᵢ = E[Yⁱ], and μᵢ = i! makes C(m,i)²(m−i)!i! = m!·C(m,i), collapsing the sum to m!·Σ(−1)ⁱC(m,i) = 0. The vanishing is the binomial theorem. (IV) AND n = 4 IS FORCED, not chosen: the triangular system determines μᵢ = i! UNIQUELY, so Y must be Exp(1); a single real Gaussian squared has moments (2i−1)!! ≠ i! for any scaling, so Y costs two real Gaussians and Z two more. (V) A ONE-PARAMETER FAMILY moves the pole: P_a = (1+aZ₂)(conj(Z₁)(1−aZ₁)+conj(Z₂)) gives E[Q P_a^m] = a^{m−1}m!, i.e. E[Q exp(tP_a)] = t/(1−at) against E[exp(tP_a)] ≡ 1 — an mgf identically 1 that acquires a simple pole at ANY prescribed t = 1/a"
status: >
  (I) VERIFIED-EXACT to m = 7 by termwise complex-Gaussian moment arithmetic (no
  sampling); the interpretation question is settled by computing both readings.
  (II) VERIFIED-EXACT to m = 9, found by exhaustive search over all subsets of the
  published six terms.  A 5-term witness also exists.
  (III),(IV) PROVED — a binomial identity and a triangular linear system, both
  elementary and both checked symbolically.
  (V) VERIFIED-EXACT to m = 8 at a = 1, 2, 1/2, −1, 3, −3/2.
  WITHDRAWN in the same session: my own "quadratic counterexamples", produced by a
  sieve that only imposed E[P^m] = 0 for m ≤ 3.  Tested properly they satisfy the
  MATHIEU PROPERTY — every Q tried has E[QP^m] = 0 by m ≈ 4 — so they are not
  counterexamples at all.  Recorded because the death depth is the useful number.
  HONEST: the counterexample is the outside source's.  Items II–V are this session's.
source: kind-pasteur-2026-07-20-S128c117 (owner: verify an outside GMC(4) counterexample, say whether we had it or anything stronger, and improve it)
depends_on: []
related: [THM-1430]
script: 04-computation/gmc4_verify_kps_S128c117.py, gmc4_extend_kps_S128c117.py, gmc4_family_kps_S128c117.py, gmc4_quadratic_decide_kps_S128c117.py, gmc4_four_term_mechanism_kps_S128c117.py (+ .out)
---

# THM-1525 — GMC(4): confirmed, sharpened to four terms, and accounted for

## 0. Did we have it? No.

`00-navigation/INVESTIGATION-BACKLOG.md:10293` is explicit: transporting `F` to a Zhao
vanishing-conjecture / image-conjecture / Mathieu-subspace witness is *"OPEN, unclaimed"*
and *"Nobody (repo or literature) has an explicit witness."* My own THM-1430 got as far as
an explicit symmetric-Jacobian Keller counterexample on `ℂ⁶` and said plainly it was **not**
a VC witness, because skipping BCW leaves `P` inhomogeneous and `Hess P` non-nilpotent. The
route we had planned would have landed at `N = 158` with a quartic of 1660 monomials.

The outside construction reaches the same target at **`n = 4` with a cubic**. It is
strictly better than anything we had, and nothing in canon anticipates it.

## I. Verified — and the reading was ambiguous

With `Z₁, Z₂` standard complex Gaussians and `E[∏ Z^a conj(Z)^b] = ∏ δ_{ab} a!`:

| m | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| `E[P^m]` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `E[QP^m]` | 1 | 2 | 6 | 24 | 120 | 720 | 5040 |

"4 real Gaussians via complex `Z_j, W_j`" admits two readings and only one survives: if
`Z₁,Z₂,W₁,W₂` were independent, `P` would carry `W`'s with no `W̄`'s, every monomial of
`QP^m` would have unbalanced `W`-charge, and **both** expectations would vanish
identically — no `m!` possible. So `W_j = conj(Z_j)` is forced, which is also what makes
the count 4 real.

## II. Four terms, not six

Exhaustive search over subsets of the published six terms:

> **`P₄ = (1 + Z)(conj(Z) − |Z₁|²)`** — four terms, still cubic, still `n = 4`,
> with `E[P₄^m] = 0` and `E[Z P₄^m] = m!` verified to `m = 9`.

(A 5-term witness `(1+Z₂)(conj(Z₂) − |Z₁|²) + conj(Z₁)` also exists.) So `Z₁` enters only
through `|Z₁|²`, and the object is really *one complex Gaussian times an independent
exponential*.

## III. Why it works — the binomial theorem

Write `P = (1+Z)(Z̄ − Y)` with `Y` independent, `μ_i = E[Y^i]`. Since
`E[Z^j Z̄^k] = δ_{jk} j!`,

> `E[P^m] = Σ_{i=0}^{m} (−1)^i C(m,i)² (m−i)! μ_i`.  (∗)

If `μ_i = i!` then `C(m,i)² (m−i)! i! = m!·C(m,i)` exactly, so (∗) collapses to

> `m! · Σ_i (−1)^i C(m,i) = 0`.

**The vanishing is nothing but the alternating binomial sum.** Both the identity and the
collapse are checked symbolically.

## IV. n = 4 is forced by the construction

Read (∗) `= 0` as a triangular system in `μ`. It determines the moments one at a time:
`μ₁ = 1, μ₂ = 2, μ₃ = 6, μ₄ = 24, …` — **`μ_i = i!` is the unique solution**, so `Y` must
be `Exp(1)`.

Could `Y` come from one real Gaussian, giving `n = 3`? No: `Y = cX²` has
`E[Y^i] = c^i(2i−1)!!`, and `(2i−1)!! ≠ i!` beyond `i = 1` for every `c` — already `3` vs
`2` at `i = 2`. So `Y` needs a genuine `Exp(1)` (two real Gaussians) and `Z` two more.

> **`n = 4` is not a choice in this construction; it is forced. Beating it requires a
> different mechanism, not a cheaper `Y`.**

## V. The pole is tunable

`E[P^m] = 0 ∀m` says exactly `E[exp(tP)] ≡ 1`; `E[QP^m] = m!` says
`E[Q exp(tP)] = Σ_{m≥1} t^m = t/(1−t)`. So the counterexample is *an mgf identically 1
that acquires a simple pole at `t = 1` once tilted by `Q`* — the failure sits at a finite
`t`, not at infinity.

And it moves. Along the diagonal of `P_a = (1+aZ₂)(conj(Z₁)(1−aZ₁)+conj(Z₂))`:

> `E[P_a^m] = 0` and `E[Q P_a^m] = a^{m−1} m!`, hence `E[Q exp(tP_a)] = t/(1−at)`,

verified at `a = 1, 2, ½, −1, 3, −3/2` to `m = 8`. **The pole can be placed at any
non-zero real `1/a`.**

## Honest negatives from the same session

- The naive three-variable chain `(1+Z₃)((1+Z₂)(…)+conj(Z₃))` **does not** extend the
  construction: `E[P₃^m] = m! ≠ 0`, so it is not even in the vanishing set.
- **My own "quadratic counterexamples" are withdrawn.** They came from a sieve imposing
  `E[P^m] = 0` only for `m ≤ 3` — three conditions on a large coefficient space, so hits
  are expected by chance. Tested properly, they do satisfy `E[P^m] = 0` to `m = 11`, but
  **every** `Q` tried has `E[QP^m] = 0` by `m ≈ 4`, so the Mathieu property *holds* for
  them. Not counterexamples. The useful number is the depth: a sieve on `E[P^m]` alone,
  at `m ≤ 3`, means nothing here — and the second half of the definition (some `Q` with
  non-vanishing moments) is where the content is.

## Named next

- Is 4 terms minimal? Only subsets of the published `P` were searched; a different cubic
  might do it in 3.
- Is `n = 3` reachable at all? §IV kills it for *this* mechanism only. The general
  question — the least `n` with GMC(`n`) false — is open and now sharply posed.
- The `Exp(1)` in §IV is the same `i!` that drives `E[Z^m conj(Z)^m] = m!`. Whether that
  coincidence is structural or bookkeeping is worth one pass.
