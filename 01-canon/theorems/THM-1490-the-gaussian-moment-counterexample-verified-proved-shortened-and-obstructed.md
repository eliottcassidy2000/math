---
id: THM-1490
title: "THE OUTSIDE GAUSSIAN-MOMENT COUNTEREXAMPLE: INDEPENDENTLY VERIFIED, THEN PROVED IN CLOSED FORM, SHORTENED, AND EXPLAINED. (0) NO, WE DID NOT ALREADY HAVE IT, and we had strictly LESS in the direction that matters — THM-1435 brackets the Zhao-VC witness dimension at 5 ≤ vcwd ≤ ~20 and explicitly does not produce a witness. (1) VERIFIED EXACTLY by Wick calculus, not sampling — but the claim AS STATED IS AMBIGUOUS: 'via complex Z_j, W_j' in 4 real Gaussians admits several conjugate pairings and only W_j = conj(Z_j) works; under Z₂ = conj(Z₁) it FAILS outright (E[P²] = −2 ≠ 0). (2) PROVED FOR ALL m, not checked to m = 10: E[P^m] = m!·(1−1)^m and E[QP^m] = m!·(1 − (1−1)^m) come from ONE binomial computation — the counterexample IS the identity Σ(−1)^k C(m,k) = 0, and the Q-insertion merely deletes the k = −1 term that made the sum vanish. (3) SHORTENED 6 TERMS → 4: the constant in (1−Z₁) is INERT, because E_A[A^k(c−Ā)^k] = (−1)^k k! for every c; so P′ = (1+B̄)(B − |A|²) works, still cubic. (4) A ONE-PARAMETER FAMILY: P = (1+μB̄)(A(1−λĀ)+B) gives E[P^m] = m!(μ−λ)^m, so the whole diagonal λ = μ ≠ 0 refutes GMC, with E[QP^m] = m!·μ^{m−1}. (5) A PROVED OBSTRUCTION explaining why nothing smaller was found: the Gaussian torus action grades monomials by weight, E ≠ 0 forces weight 0, hence single-weight P is useless, 0 must lie in conv(weights of P), and real-valued P is trivial. Bounded exhaustive searches (2 real Gaussians; ≤3 monomials; degree ≤2) come back empty, with a positive control that fires."
status: >
  (0) ASSESSMENT, checked against the repo's own files.
  (1) VERIFIED-EXACT.  Integer Wick arithmetic, m = 1..10, no sampling.  Five conjugate
      pairings tested; exactly two (mutually conjugate) reproduce the claim.
  (2) PROVED for all m >= 1.  Elementary and complete; machine-checked against brute
      expansion at m = 1..8.
  (3) PROVED (the c-independence is one line) + VERIFIED m = 1..10.
  (4) PROVED by the same computation + VERIFIED at mu = 1,2,3,-1.
  (5) The three structural statements are PROVED (one line each) and machine-checked
      against two deliberate violators.  The MINIMALITY claims are BOUNDED SEARCHES, NOT
      PROOFS, and the bounds are stated explicitly below.
  NO PRIORITY CLAIM.  The counterexample is not ours; see attribution.
source: klein-2026-07-20-S341 (owner: take the outside result disproving the Gaussian moment conjecture at n>=4, determine if we already had it or anything stronger, then explore and improve)
attribution: >
  THE COUNTEREXAMPLE IS NOT OURS.  P = (1+Z2)(W1(1-Z1)+W2), Q = Z2 was supplied by the
  owner from an outside source.  Everything claimed here is verification, proof,
  simplification and obstruction ON TOP of someone else's construction.  This sits in the
  same field-wide scramble recorded in THM-1300's attribution blocks (Alpoge/Mathew,
  19 July 2026); no priority is available or claimed on any of it.
depends_on: []
related:
  - THM-1435  # klein-S337: the Zhao-VC witness machinery and the 5 <= vcwd <= ~20 bracket
  - THM-1430  # kind-pasteur: the explicit symmetric Keller counterexample on C^6
  - THM-1300  # Alpoge's Keller map + the attribution record
script: 04-computation/gmc_counterexample_klein_S341.py, gmc_minimality_klein_S341.py (+ .outs)
---

# THM-1490 — the Gaussian moment counterexample: verified, proved, shortened, obstructed

**GMC(n).** For the standard Gaussian on `ℝ^n` and `P, Q ∈ ℂ[x_1,…,x_n]`: if `E[P^m] = 0`
for all `m ≥ 1`, then `E[Q·P^m] = 0` for all `m ≫ 0`. (The Mathieu–Zhao subspace property of
`ker E`.) The outside claim refutes it at `n = 4`.

## 0. Did we already have this? No — and we had less

Asked directly, the answer is no, and the comparison is worth stating precisely because it
lands on the repo's own witness-effective/truth-only axis.

| what the repo holds | what it is |
|---|---|
| THM-1300 | Alpöge's Keller map, correctly attributed; notes Zhao's VC, the image conjecture and Mathieu subspaces as corollaries "being chased by the whole field" |
| THM-1430 (kind-pasteur) | an explicit symmetric Keller counterexample on `ℂ⁶` |
| THM-1435 (klein-S337) | Zhao-VC transport machinery, control-validated; **`5 ≤ vcwd ≤ ~20`**; *the witness is not produced* |

So we held an *inference* plus *machinery* plus a *bracket* — and no witness. The outside
result is an explicit witness, in dimension **4**, with 4 monomials and degree 3. On the
axis this repo keeps score on, that is strictly better than anything here.

**One caveat that must not be blurred: GMC is not Zhao's Vanishing Conjecture.** They are
different members of the Mathieu–Zhao family. This witness does **not** close `vcwd`, and
THM-1435's bracket stands untouched. Anyone tempted to write "so vcwd = 4" should not.

## 1. Verified — and the statement needs repairing first

"in 4 real Gaussians via complex `Z_j, W_j`" is ambiguous: four complex symbols in four
*real* Gaussians forces conjugate identifications, and the statement does not say which.
Testing all five natural pairings by exact Wick arithmetic:

| pairing | `E[P^m]`, m=1..5 | verdict |
|---|---|---|
| `Z₂ = conj(Z₁)`, `W₂ = conj(W₁)` | `0, −2, 0, 60, 0` | **fails** |
| `W_j = conj(Z_j)` | `0,0,0,0,0` | **works**, `E[QP^m] = m!` |
| `Z₁=conj(A), Z₂=A, …` | `0, −2, 0, 60, 0` | fails |
| `Z₁=A, Z₂=conj(B), W₁=B, W₂=conj(A)` | `1,−2,−30,−120,1320` | fails |
| `W_j = conj(Z_j)` (global conjugate) | `0,0,0,0,0` | works |

**Only `W_j = conj(Z_j)` realises the claim.** Writing `A = Z₁`, `B = Z₂` and `X̄` for
conjugation, the object is

```text
P = (1 + B̄)·( A(1 − Ā) + B ),      Q = B̄
```

with `A, B` independent standard complex Gaussians (`E[XX̄] = 1`), i.e. 4 real Gaussians.
Exactly: `E[P^m] = 0` and `E[QP^m] = m!` for `m = 1,…,10`.

## 2. The mechanism — a closed form, hence a proof for all `m`

`A ⊥ B`, so expanding `U = A(1−Ā) + B` binomially makes the expectation factor:

```text
E[P^m] = Σ_k C(m,k) · E_A[A^k(1−Ā)^k] · E_B[(1+B̄)^m B^{m−k}]
```

Both factors are one-line Wick evaluations, since `E[A^p Ā^q] = δ_{pq} p!`:

```text
E_A[A^k(1−Ā)^k]     = (−1)^k k!            (only the Ā^k term survives)
E_B[(1+B̄)^m B^{m−k}] = C(m,k)(m−k)!
```

Using the collapse `C(m,k)² k!(m−k)! = m!·C(m,k)`:

```text
E[P^m]   = m! · Σ_k (−1)^k C(m,k)     = m!·(1−1)^m         = 0    (m ≥ 1)
E[QP^m]  = m! · Σ_k (−1)^k C(m,k+1)   = m!·(1 − (1−1)^m)   = m!   (m ≥ 1)
```

> **The counterexample is the binomial theorem.** One computation gives both lines; the `Q`
> insertion shifts the binomial index by one, deleting exactly the `k = −1` term whose
> presence made the alternating sum vanish. That is the whole content.

Note this is a *strong* refutation: `E[QP^m] = m!` is never zero and grows, so no
"eventually" clause can rescue GMC. *(Machine-checked against brute expansion, `m = 1..8`.)*

## 3. Shortened: 6 terms → 4

`E_A[A^k(c − Ā)^k] = (−1)^k k!` for **every** constant `c`, since only the `Ā^k` term
survives. So the `1` in `(1 − Z₁)` is inert, and `c = 0` is legal:

```text
P′ = (1 + B̄)(B − |A|²),      Q = B̄        — 4 monomials, still cubic
```

Verified: `E[P′^m] = 0` and `E[QP′^m] = m!` for `m = 1..10`. The `A`-block enters only
through `|A|² ~ Exp(1)`, whose moments `k!` are the entire role it plays.

## 4. It is a line, not a point

```text
P_{λ,μ} = (1 + μB̄)( A(1 − λĀ) + B )   ⟹   E[P^m] = m!·(μ − λ)^m
```

so **every** `λ = μ` kills every moment, and then `E[QP^m] = m!·μ^{m−1}`, nonzero for
`μ ≠ 0`. Verified at `μ = 1, 2, 3, −1`. The published example is the point `μ = 1` of a
one-parameter family.

## 5. Why nothing smaller was found — a proved obstruction

The Gaussian on `ℂ^K` is invariant under the torus `A_j ↦ e^{iθ_j}A_j`. Grade the monomial
`∏ A_j^{a_j} Ā_j^{b_j}` by its **weight** `w = (a_j − b_j)_j`. Weights add under
multiplication, and Wick gives `E[monomial] ≠ 0 ⟹ w = 0`. Three consequences, one line each:

1. **Single-weight `P` is useless.** If all monomials of `P` share a weight `v ≠ 0`, then
   `QP^m` has weights `w_Q + mv`; a fixed `Q` has bounded weights, so no weight is `0` for
   large `m`. (This kills the cheap "take `P` holomorphic" idea outright.)
2. **`0 ∈ conv(weights of P)` is necessary** — otherwise all weights lie in an open
   half-space `⟨w,u⟩ > 0`, so `P^m` has `⟨w,u⟩ ≥ mc` with `c > 0`, and `Q` only shifts
   boundedly.
3. **Real-valued `P` is trivial**: `E[P²] = 0` forces `P ≡ 0`. GMC has content only for
   *complex* `P` — which is precisely why the witness is built from complex combinations.

The witness has weights `{(0,1), (0,0), (0,−1)}`: it straddles `0` in the `B` direction and
is identically `0` in the `A` direction. Two deliberate violators behave as predicted —
`A + B` (holomorphic) has `E[P^m] = 0` but **no** `Q` works, and `|A|² − |B|²` (weight 0)
fails at `E[P²] = 2 ≠ 0`.

**Bounded minimality searches** — exhaustive in the stated boxes, with the known 4-term
witness as a firing positive control:

| search | box | result |
|---|---|---|
| 2 real Gaussians (one complex) | ≤4 monomials, degree ≤4, coeffs `±1,±2` | **none** |
| 4 real Gaussians, fewer terms | 2 and 3 monomials, degree ≤3, coeffs `±1,±2` | **none** |
| 4 real Gaussians, lower degree | degree ≤2, ≤4 monomials, coeffs `±1,±2` | **none** |

These are **searches, not proofs** — they do not establish that GMC(2) or GMC(3) is true, and
they do not prove 4 terms minimal. What §5(1)–(3) *does* prove is that any smaller witness
must still straddle `0` in weight, be complex, and use at least two opposite weights; and the
two-term case dies by hand, since `P = c₁A + c₂Ā` gives `E[P^{2k}] = C(2k,k)(c₁c₂)^k k! ≠ 0`.

## 6. Open, in priority order

- **Is GMC(2) or GMC(3) true?** The obstruction constrains but does not decide. A proof at
  `n ≤ 3` would make the outside witness dimension-minimal.
- **Is 4 monomials minimal at `n = 4`?** Same status.
- **Does a GMC witness transport to a Zhao-VC witness?** If the two Mathieu–Zhao members are
  linked constructively, this would bear on THM-1435's `5 ≤ vcwd ≤ ~20`. *Not* assumed here.

*Files: `04-computation/gmc_counterexample_klein_S341.py`, `gmc_minimality_klein_S341.py` (+ `.out`s).*
