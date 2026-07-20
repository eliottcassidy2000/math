---
id: THM-1570
title: "THE ONE-SIDED CONJECTURE PROVED BY EXACT ELIMINATION (charge span {-1,0,+1}, degree <= 1, over ALL of C — not a grid sweep), AND THE LAPLACE / GMC(1) LAYER SETTLED. (A) Every P of charge span {-1,0,+1} is P = z*A(s) + B(s) + zbar*C(s) with s = z*zbar; 'one-sided' means A == 0 or C == 0. Rabinowitsch saturation — adjoin 1 - t*a_i*c_j to the moment ideal and test for <1> — returns the TRIVIAL ideal for every (i,j) at degree <= 1, so NO nullcone point has both a nonzero A-coefficient and a nonzero C-coefficient. This replaces every previous finite-grid verification with a statement quantified over all complex coefficients. The moment equations display exactly the predicted structure: E[P] = b0+b1 is the charge-0 part, and E[P^2] = (2a0c0+4a0c1+4a1c0+12a1c1) + (b0^2+2b0b1+2b1^2) splits as the OPPOSITE-CHARGE cross term plus the charge-0 Hankel form, confirming THM-1535 s3 / THM-1540 in coordinates. Degree 2 (9 unknowns) exceeded the Groebner budget — a compute limit, not a negative result. (B) LAPLACE LAYER: the only g in C[s] with int_0^inf g^m e^{-s} ds = 0 for all m >= 1 is g = 0. Proof for all degrees: d=0 gives c0^m; for d >= 1 the saddle of |g|^m e^{-s} sits at s ~ dm, where m*a1/s -> a1/d is constant and arg g(s) -> arg(c_d) — THE PHASE STABILISES — so Laplace's method gives I_m ~ c_d^m (dm)! e^{a1/d} != 0. Verified to 8 digits, including the COMPLEX case g = i*s - 1 whose ratio converges to exactly e^i = 0.54030231 + 0.84147098i. Hence the Laplace nullcone is TRIVIAL and GMC(1) for the exponential measure holds VACUOUSLY"
status: PROVED — (A) by exact elimination at degree <= 1 over C; (B) for all degrees (saddle-point, numerically confirmed to 8 digits incl. complex)
author: opus-2026-07-20-S413
depends_on: [THM-1540, THM-1535, THM-1495]
---

# THM-1570 — One-sided by exact elimination, and the Laplace layer settled

## A. The one-sided conjecture, by exact elimination

Every `P` at `n = 2` with charges in `{−1, 0, +1}` is

```
P  =  z·A(s)  +  B(s)  +  \bar z·C(s) ,        s = z\bar z
```

with `A, B, C ∈ ℂ[s]` carrying charges `+1, 0, −1`. **"One-sided" means `A ≡ 0` or
`C ≡ 0`.** So for this charge span the conjecture is exactly

```
V( E[P¹], …, E[P^M] )  =  {A = 0} ∪ {C = 0}      over ℂ.
```

**Method — Rabinowitsch saturation.** Adjoin `1 − t·a_i·c_j` to the moment ideal and test
whether the result is `⟨1⟩`. If it is, no point of the nullcone has both `a_i ≠ 0` and
`c_j ≠ 0`.

**Result (degree ≤ 1, 6 unknowns, 6 moment equations): every saturation returns `⟨1⟩`.**

> **Every nullcone point of charge span `{−1,0,+1}` and degree `≤ 1` has `A ≡ 0` or
> `C ≡ 0` — proved over ALL of `ℂ`.**

This matters methodologically as much as mathematically. Every prior check of this
conjecture (THM-1535's 59048 real polynomials, THM-1540's Gaussian-integer sweeps) ranged
over a **finite grid of coefficients**; this is quantified over the whole coefficient space.
It is the fix for the caution I recorded in S412 — *name the domain of the sweep* — applied
to my own work.

**The equations confirm the predicted structure**, in coordinates:

```
E[P ] = b₀ + b₁                                             ← the charge-0 part alone
E[P²] = (2a₀c₀ + 4a₀c₁ + 4a₁c₀ + 12a₁c₁)  +  (b₀² + 2b₀b₁ + 2b₁²)
        └── opposite-charge cross term ──┘     └── charge-0 Hankel form ──┘
```

exactly the split `E[P²] = E[P₀²] + 2E[P₊P₋]` that THM-1535 §3 identified as the mechanism
separating `n = 2` from `n = 4`. At `n = 2` the cross term is bilinear in `(A, C)`, so
killing it forces one factor to vanish; at `n = 4` the second lattice direction supplies an
independent cross term that can cancel the Hankel form instead.

**Degree 2** (9 unknowns) exceeded the Gröbner budget in this environment and did not
finish. That is a **compute limit, not a negative result** — nothing was found to fail.

## B. The Laplace / GMC(1) layer, settled

> **Theorem.** The only `g ∈ ℂ[s]` with `∫₀^∞ g(s)^m e^{−s} ds = 0` for all `m ≥ 1` is
> `g = 0`.

*Proof.* If `d = deg g = 0` then `g = c₀` and the integral is `c₀^m`, forcing `c₀ = 0`.

Let `d ≥ 1`, `g(s) = c_d s^d (1 + a₁/s + O(1/s²))`. The weight `|g|^m e^{−s}` has
log-derivative `m·g′/g − 1`, and `g′/g ∼ d/s`, so the **saddle sits at `s ∼ dm`**. There
`m·a₁/s → a₁/d`, a *constant*, and `arg g(s) → arg c_d`: **the phase stabilises**, so no
oscillation survives to cancel the tail. Laplace's method then gives

```
∫₀^∞ g^m e^{−s} ds  ∼  c_d^m · (dm)! · e^{a₁/d}   ≠ 0   for large m. ∎
```

**Numerical confirmation of the constant `e^{a₁/d}`** (ratio `I_m / (c_d^m (dm)!)`):

| `g` | `d` | predicted `e^{a₁/d}` | `m = 18` |
|---|---|---|---|
| `s − 1` | 1 | `e^{−1} = 0.36787944` | **0.36787944** |
| `s + 3` | 1 | `e³ = 20.085537` | **20.085537** |
| `2s² − s` | 2 | `e^{−1/4} = 0.7788008` | 0.7781 (converging) |
| **`i·s − 1`** | 1 | `e^{i} = 0.54030231 + 0.84147098i` | **0.54030231 + 0.84147098i** |

The complex row is the one that matters: it is exactly where phase cancellation *could* have
occurred, and the ratio lands on `e^i` to 8 digits. **The phase-stabilisation mechanism is
real, not an artifact of positivity.**

> **Corollary.** The nullcone for the exponential (Laplace) measure on `[0,∞)` is **trivial**
> — `{0}`. Hence **GMC(1) for that measure holds vacuously**: there is no nonzero `P` whose
> moments all vanish, so the hypothesis is never satisfied.

This also completes the argument of THM-1540 §3, which had the degree ≤ 3 cases by exact
solve and the general case only sketched. **The Laplace layer is now closed for all degrees.**

## C. Where GMC(2) now stands

| locus | status |
|---|---|
| charge-definite, any `n` | **PROVED** (THM-1535 §1) |
| sign-coherent over `ℝ` | **PROVED** (THM-1535 §2) |
| sign-coherent over `ℂ` | **PROVED** (THM-1540 §3, completed by §B here) |
| **charge span `{−1,0,+1}`, deg ≤ 1, both signs** | **PROVED over all of `ℂ`** (§A) |
| charge span `{−1,0,+1}`, deg ≥ 2 | open (compute-limited only) |
| wider charge spans | open |

The remaining gap is now *stratified by charge span and degree* rather than being a single
undifferentiated case — and the first genuinely two-sided stratum is closed.

## D. Open

1. **Degree 2 and beyond at span `{−1,0,+1}`.** Needs a better Gröbner engine (Singular /
   msolve) or a structural argument replacing elimination. The bilinearity of the cross term
   in `(A, C)` looks like the lever: `E[P²]`'s cross part is `AᵀMC` for an explicit matrix
   `M`, and one wants `AᵀMC = 0` plus higher moments to force `A = 0` or `C = 0`.
2. **Wider spans `{−k,…,k}`.** The same setup with `2k+1` blocks.
3. **Is `M` (the cross-term matrix) always nonsingular?** At degree ≤ 1 it is
   `[[2,4],[4,12]]`, determinant `8 ≠ 0`. If nonsingularity holds in general it may give the
   structural argument directly, without elimination.

## Verification

`04-computation/one_sided_elimination_opus_S413.py` (saturation proof at degree ≤ 1; the
Laplace asymptotic checks), `04-computation/one_sided_deg2_attempt_opus_S413.py` (the degree-2
attempt, compute-limited). Output in `05-knowledge/results/`.
