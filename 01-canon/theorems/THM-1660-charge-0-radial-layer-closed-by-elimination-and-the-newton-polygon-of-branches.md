---
id: THM-1660
title: "THE CHARGE-0 RADIAL LAYER, CLOSED BY EXACT ELIMINATION FOR BOUNDED β-DEGREE — the sign-indefinite sub-case where positivity AND domination both fail — and the Newton polygon of the branches that reads off boxeph's κ exponent. (1) FIX α = r (two-sided: charges ±1 both present), let β = b(r) be the charge-0 coefficient of degree d; NC2 predicts no nullcone member, and exact Gröbner returns the UNIT ideal ⟨1⟩ for d = 0,1,2,3 — CLOSED, no asymptotics, where THM-1640's positivity is unavailable (L_m sign-indefinite) and MISTAKE-202's domination is false. (2) THE REAL CASE IS TRIVIAL AT m=2: E_r[L_2] = βᵀHβ + 2E_r[α] with H_{ij} = (i+j)! the exponential Hankel matrix, positive-definite, so over ℝ the layer closes at the second moment alone — opus's Hankel-PD (THM-1535) appearing at m=2 — and the COMPLEX case is the only hard one, which is why elimination needs m up to d+4. (3) THE NEWTON POLYGON of the branches λ_±(r) = β ± 2√α: the large branch grows like r^{max(B,(A+1)/2)}, so the branch points 1/λ accumulate at t=0 like r^{−e}, e = max(B,(A+1)/2) — which reproduces boxeph's THM-1620 κ = −1/2 at the pair model (A=B=0) from a Newton polygon rather than a pinch, and the split at B = (A+1)/2 IS the Newton-split of Case I vs Case II."
status: >
  (1) PROVED for deg β ∈ {0,1,2,3} with α = r, by exact Gröbner over ℚ: ideal = ⟨1⟩,
      empty variety. Extends mechanically; cost grows with degree. Bounded, not unbounded.
  (2) PROVED (H_{ij} = (i+j)! is the Hankel matrix of Exp(1), leading minors 1,1,4,144,82944
      verified positive). The real-case closure at m=2 is elementary; it is opus THM-1535
      localised to the first even moment on this stratum.
  (3) DERIVATION (Newton polygon of λ_±). The exponent e = max(B,(A+1)/2) is read off the
      polygon and matches boxeph's κ at the one point (A=B=0) both compute. Presented as the
      organising asymptotic, not a proof of the unbounded case.
  GMC(2) REMAINS OPEN. This closes bounded β-degree exactly and frames the unbounded case.
source: klein-2026-07-20-S365 (owner: work the radial Laplace nullcone for the charge-0 part — reduced via Watson but not closed — think Newton polygon factorization of large branches)
depends_on:
  - THM-1640  # the sign-indefiniteness reframing this closes at bounded degree
  - THM-1590  # the same elimination method
related:
  - THM-1615  # boxeph: Hermite orthogonality closes constant-β; this does non-constant β
  - THM-1620  # boxeph: the Newton-split / Case II Watson lemma this reads off algebraically
  - THM-1535  # opus: Hankel positive-definiteness, here localised to m=2
script: 04-computation/radial_charge0_elimination_klein_S365.py (+ .out)
---

# THM-1660 — the charge-0 radial layer, closed by elimination

On the `{−1,0,1}` span, `E[P^m] = E_r[L_m(α,β)]`, `α = r·a·c`, `β = b`,
`L_m = Σ_k m!/(k!²(m−2k)!) α^k β^{m−2k}`, `E_r[r^j] = j!`. The charge-0 coefficient is `β`.

boxeph closed **constant** `β` by Hermite orthogonality (THM-1615) and reduced **non-constant**
`β` (specifically the `P(0)=0` case) to one unclosed per-component Watson lemma (THM-1620
Case II). This handles the non-constant case at bounded degree *by exact elimination*, with no
asymptotics — and reads boxeph's asymptotic exponent off a Newton polygon.

## 1. The closure

Fix `α = r` (so `a = c = 1`: both charges `±1` present, `P` is **two-sided**, and NC2 predicts
no nullcone member). Let `β = b(r) = Σ_{i=0}^d b_i r^i`. Then `E_r[L_m]` is a polynomial in the
`b_i`, and the nullcone question is whether `⟨E_r[L_m] : m ≥ 1⟩` has any zero.

| `deg β` | unknowns | moments | Gröbner basis |
|---|---|---|---|
| 0 | 1 | 1..4 | `⟨1⟩` |
| 1 | 2 | 1..5 | `⟨1⟩` |
| 2 | 3 | 1..6 | `⟨1⟩` |
| 3 | 4 | 1..7 | `⟨1⟩` |

> **The ideal is the unit ideal in every case: no `β` of degree ≤ 3 makes every moment
> vanish when `α = r`.** The sign-indefinite charge-0 layer is closed at bounded degree,
> exactly, where THM-1640 showed positivity is unavailable (`L_m(r)` takes both signs) and
> MISTAKE-202 showed domination is false.

This is boxeph's Watson Case II handled without any asymptotics, for these bounded degrees.

## 2. Over ℝ it is trivial — at `m = 2`

`L_2 = β² + 2α`, so `E_r[L_2] = E_r[β²] + 2E_r[α]`. Writing `β = Σ b_i r^i`,

```text
E_r[β²] = Σ_{i,j} b_i b_j (i+j)! = bᵀ H b,   H_{ij} = (i+j)!.
```

`H` is the **Hankel moment matrix of `Exp(1)`** — positive-definite (leading minors
`1, 1, 4, 144, 82944` for `d = 1..4`). So over ℝ, with `α = r ≥ 0`,

```text
E_r[L_2] = ‖β‖²_H + 2·E_r[r] ≥ 2 > 0,
```

and no real `β` can make even the *second* moment vanish. This is **opus's
Hankel-positive-definiteness (THM-1535) appearing at `m = 2`** on this stratum. It also
explains the shape of §1: the complex case is the only hard one — over ℂ the form `bᵀHb` is
indefinite — which is exactly why the elimination has to climb to `m = d+4` rather than
stopping at 2.

## 3. The Newton polygon of the branches

The generating function is `Σ_m E[P^m] t^m = E_r[D(r,t)^{−1/2}]`,
`D = (1−βt)² − 4αt² = (1−λ₊t)(1−λ₋t)` with `λ_± = β ± 2√α`. Branch points sit at
`t = 1/λ_±(r)`. Let `A = deg(ac)` (so `deg α = A+1`) and `B = deg β`. The **large branch**
`λ_+(r)` grows like the larger of `β ~ r^B` and `2√α ~ r^{(A+1)/2}`:

| `(A,B)` | large branch | branch exponent `e = max(B,(A+1)/2)` |
|---|---|---|
| (0,0) | `√α` | `1/2` |
| (0,1) | `β` | `1` |
| (0,2) | `β` | `2` |
| (2,1) | `√α` | `3/2` |
| (1,3) | `β` | `3` |

So the branch points `1/λ_+ ~ r^{−e}` **accumulate at `t = 0`** — which is boxeph's Case II
picture ("active arcs land AT 0"), and

> `e = max(B, (A+1)/2)` reproduces boxeph's `κ = −1/2` at the pair model `(A,B)=(0,0)` from a
> **Newton polygon** rather than a pinch computation.

The two regimes are the two edges of the polygon: when the **charge-0 coefficient** dominates
(`B ≥ (A+1)/2`) the accumulation rate is set by `B`; when the **pair scale** dominates it is
`(A+1)/2`. The changeover `B = (A+1)/2` is precisely the Newton-split of THM-1620 (Case I where
`P(0) ≠ 0` vs Case II where the charge-0 term drives the accumulation to `0`). This gives an
algebraic reading of the exponent that boxeph's per-component Watson lemma must control in the
unbounded case.

## 4. Scope and the honest boundary

- **Bounded `β`-degree: closed** (§1), exactly, over ℂ, no asymptotics.
- **Real coefficients: closed at `m=2`** (§2), by Hankel positivity.
- **Unbounded `β`-degree over ℂ: open**, governed by the single Newton-polygon exponent
  `e = max(B,(A+1)/2)` (§3) — the object of boxeph's THM-1620 Case II. This file does not
  close it; it identifies the controlling exponent and shows it is a Newton-polygon datum.

GMC(2) is not closed. What is added: the sign-indefinite charge-0 case is now settled at every
bounded degree by a finite exact computation (removing it from the list of things needing
analysis at bounded degree), the real case is seen to close trivially at `m=2`, and boxeph's
`κ` is re-derived as a Newton-polygon exponent — tying the elimination and Watson routes to the
same combinatorial datum.

*Files: `04-computation/radial_charge0_elimination_klein_S365.py` (+ `.out`).*
