---
id: THM-1590
title: "THE ONE-SIDED CONJECTURE, PROVED BY EXACT ELIMINATION ON BOUNDED CHARGE SPAN — and the Laplace/GMC(1) layer settled the same way. With w = e^{iθ}, E[P^m] = ∫₀^∞ CT_w(P(r,w)^m)e^{−r}dr, and on the {−1,0,1} span CT_w(P^m) = L_m(α,β) = Σ_k m!/(k!²(m−2k)!)α^kβ^{m−2k} with α = r·a·c, β = b. Since L_m sees (a,b,c) only through (α,β), bounding deg α ≤ A and deg β ≤ B makes the whole conjecture a FINITE polynomial system — E_r[L_m] = 0 with E_r[r^j] = j! — and a Gröbner basis decides it outright. PROVED for all seven strata with A,B ≤ 2: the variety is exactly {α = β = 0}, i.e. P is one-sided. No sampling, no tolerance, no asymptotics. Separately the pure Laplace layer (β ≡ 0, where E[P^{2j}] = C(2j,j)E_r[α^j]) is settled by the same elimination for deg ≤ 3 — EMP re-proved with no tail estimate."
status: >
  PROVED on the stated bounded strata, by exact Gröbner elimination over Q.
  {-1,0,1} span: (A,B) in {(0,0),(1,0),(0,1),(1,1),(2,1),(1,2),(2,2)} -- variety = origin only.
  Laplace/GMC(1) layer: deg h <= 3 -- variety = origin only.
  Stability checked: the (2,1), (1,2), (2,2) verdicts are unchanged at m = 7, 11, 15
  (resp. 8, 12, 16) moments, so they are not truncation artefacts.
  NOT a proof of GMC(2): the strata are bounded in BOTH the charge span and the coefficient
  degrees.  What it removes is the reliance on asymptotics within those bounds.
  A first run reported three of these strata as "nontrivial"; that was an artefact of MY
  radical test, not of the mathematics.  See §4.
source: klein-2026-07-20-S357 (owner: prove the one-sided conjecture for bounded charge span by exact elimination and settle the Laplace-GMC(1) layer; with w = e^{iθ}, show ∫₀^∞ CT_w(P^m)e^{−s}ds ≠ 0 for some m)
depends_on:
  - THM-1510  # EMP and the two-weight theorem -- this re-proves EMP by elimination
  - THM-1550  # the toral criterion; CT_w is the same projection used there
related:
  - THM-1530  # the toral nullcone framing
  - "boxeph: the One-Sided Nullcone Conjecture; opus THM-1535/1540: the charge lattice"
script: 04-computation/exact_elimination_klein_S357.py (+ .out)
---

# THM-1590 — one-sided by exact elimination

## 1. The integral, and why it is finite-dimensional

With `Z = √r·e^{iθ}`, `r ~ Exp(1)` independent of `θ ~ Unif`, and `w = e^{iθ}`, the polynomial
`P` is a Laurent polynomial in `w` with `r`-dependent coefficients and

```text
E[P^m] = ∫₀^∞ CT_w( P(r,w)^m ) e^{−r} dr
```

— the owner's integral. On the charge span `{−1,0,1}`, writing `P = Z̄a(r) + b(r) + Zc(r)`,
the `CT_w` step is classical:

```text
CT_w(P^m) = L_m(α, β) := Σ_k m!/(k!²(m−2k)!) α^k β^{m−2k},    α = r·a(r)c(r),  β = b(r).
```

**`L_m` depends on `(a,b,c)` only through the pair `(α, β)`.** So bounding `deg α ≤ A` and
`deg β ≤ B` leaves `(A+1)+(B+1)` unknowns, and the conjecture becomes a finite polynomial
system — `E_r[L_m] = 0` with `E_r[r^j] = j!` — decidable by elimination.

## 2. The elimination

| `(A,B)` | unknowns | moments used | verdict |
|---|---|---|---|
| (0,0) | 2 | 4 | **variety = {0}** |
| (1,0) | 3 | 5 | **variety = {0}** |
| (0,1) | 3 | 5 | **variety = {0}** |
| (1,1) | 4 | 6 | **variety = {0}** |
| (2,1) | 5 | 7, 11, 15 | **variety = {0}** |
| (1,2) | 5 | 7, 11, 15 | **variety = {0}** |
| (2,2) | 6 | 8, 12, 16 | **variety = {0}** |

> **Theorem.** On the `{−1,0,1}` charge span with `deg α ≤ 2` and `deg β ≤ 2`, if
> `∫₀^∞ CT_w(P^m)e^{−r}dr = 0` for every `m ≥ 1` then `α = β = 0` — i.e. `a·c ≡ 0` and
> `b ≡ 0`, so `P` is one-sided. Equivalently: for any `P` on that span that is not one-sided,
> `∫₀^∞ CT_w(P^m)e^{−r}dr ≠ 0` for some `m`.

Exact Gröbner bases over `ℚ`; no floating point anywhere. Sample bases: `(1,0)` gives
`[a₀+a₁, a₁², b₀]` and `(0,1)` gives `[2a₀+b₁², b₀+b₁, b₁³]` — visibly cutting out the origin
with multiplicity.

## 3. The Laplace / GMC(1) layer

With `β ≡ 0` the odd moments vanish by parity and `E[P^{2j}] = C(2j,j)·E_r[α^j]`, so the
`{−1,+1}` span **is** the pure radial layer: `E_r[h^j] = 0 ∀j ⟹ h = 0` for `r ~ Exp(1)`.
Elimination settles it for `deg h ≤ 3`:

| `deg h` | Gröbner basis | only `h = 0`? |
|---|---|---|
| 0 | `[h₀]` | yes |
| 1 | `[h₀+h₁, h₁²]` | yes |
| 2 | `[h₀+h₁+2h₂, h₁²+8h₁h₂ …]` | yes |
| 3 | `[h₀+h₁+2h₂+6h₃, h₁² + …]` | yes |

This is **EMP (THM-1510) re-proved by elimination instead of Laplace asymptotics** — exact and
finite, with no tail estimate to defend. Within these degree bounds it fully replaces the
saddle-point argument.

## 4. Correction — my first run, not the mathematics

The first run reported `(2,1)`, `(1,2)`, `(2,2)` as having a **nontrivial** variety. A
nontrivial solution there would be a GMC(2) *counterexample*, so it had to be settled rather
than reported. Two errors of mine, both in the test and not in the algebra:

- **The radical test was too shallow.** I checked whether `v^k` lies in the ideal only for
  `k ≤ 6`; these ideals need higher powers. At `k ≤ 12` all three flip to `variety = {0}`, and
  the verdict is stable at `m = 7, 11, 15`.
- **"Is the basis `{1}`?" was the wrong question.** `{1}` would mean *no* solutions, but
  `α = β = 0` is always a solution — so the basis is never `{1}` and that column was
  meaningless. The correct test is whether the variety is *only* the origin.

Recorded because the failure mode is specific: *when the trivial solution is always present,
inconsistency is the wrong thing to test for.*

## 5. Scope

Bounded in **both** the charge span (`{−1,0,1}`) and the coefficient degrees (`≤ 2`; `≤ 3` for
the radial layer). It is not a proof of GMC(2). What it does is discharge the asymptotic
machinery inside those bounds: on these strata the one-sided conclusion is now a finite exact
computation, and the method extends mechanically as compute allows.

## 6. Not delivered this session

The owner also asked for the relation to Erdős problems #1016, 506, 742, 19, 580, 547, 460,
556, 475, 848. The lookup task **failed with an API error** and returned nothing usable. No
statements about those problems are recorded here rather than guessed — they remain open for a
future session.

*Files: `04-computation/exact_elimination_klein_S357.py` (+ `.out`).*
