---
id: THM-1185
title: THE DELSARTE LP CANNOT PROVE LRC(14) — EVEN IN PRINCIPLE — the positive-definite test-function bound reduces to a symmetric TOEPLITZ eigenvalue problem whose matrix is M[k,l] = Σ_{v ∈ V, v | (l−k)} ĥ((l−k)/v), with diagonal 13·2λ = 13/7, and the criterion is λ_min(M) < 1. But the Toeplitz SYMBOL is f(θ) = Σ_d c_d e(dθ) = Σ_v h(vθ) — exactly the COVERING MULTIPLICITY function — so λ_min converges to its ESSENTIAL infimum, and the LP criterion is precisely "the uncovered set has POSITIVE MEASURE". That is strictly stronger than LRC(14), which permits a measure-zero lonely set, and it FAILS on all three tight families: {1,…,13}, {1,…,11,13,24} and 2·{1,…,13} have uncovered measure exactly 0 (so λ_min = 1.00000) while LRC(14) holds there with gap exactly 1/14 at the six points p/14
status: PROVED structurally (the symbol computation is exact) and verified numerically — λ_min at K = 160 reads 1.00000, 1.00000, 1.00000 on the tight families against exact uncovered measure 0, and 0.0265 / 0.1499 on non-tight families whose uncovered measure is positive. The conclusion is a genuine impossibility for this method, not a failed attempt at it
source: opus-2026-07-17-S386 (owner: work a new angle on the open LRC 14 math)
depends_on: [THM-1120/1125 (the tight families that defeat the method), THM-1075 (the Fourier frame ĥ(n) = sin(πn/7)/(πn))]
scripts: 04-computation/delsarte_lp_opus_S386.py -> 05-knowledge/results/delsarte_lp_opus_S386.out
---

# THM-1185 — the standard weapon, and why it misses

The Delsarte / positive-definite LP is the usual tool for covering and
packing problems, and it had not been tried here. It directly confronts the
13/7 obstruction with a *free optimisation* instead of a fixed bound, which
is exactly what every failed approach in this programme lacked. It still
cannot work, and the reason is structural.

## The setup

For P ≥ 0 a trig polynomial, P = |Σ a_k e(kt)|²,

> ∫P·1_U = ∫P − ∫P·1_{⋃D_v} ≥ ∫P − Σ_v ∫P·1_{D_v}.

Since 1_{D_v}(t) = h(vt), its Fourier mass sits on multiples of v, so
∫P·1_{D_v} = Σ_n P̂(nv)ĥ(n). With P̂(m) = Σ_k a_k a_{k+m} this is a quadratic
form aᵀMa in the symmetric **Toeplitz** matrix

> **M[k,l] = c_{l−k},  c_d = Σ_{v ∈ V, v | d} ĥ(d/v),  c₀ = 13·2λ = 13/7.**

So **λ_min(M) < 1 ⟹ the uncovered set has positive measure ⟹ LRC(14) for V.**

## Why it degenerates

The arcs are not arbitrary sets — they are *dilates of one function*,
D_v = h(v·). That forces the symbol to collapse:

> f(θ) = Σ_d c_d e(dθ) = Σ_v Σ_n ĥ(n)e(nvθ) = **Σ_v h(vθ)**

which is exactly the **covering multiplicity**. As K → ∞, λ_min(M) tends to
the *essential* infimum of f. So the LP criterion is not a relaxation of the
problem — it **is** the problem, in the form "some point is covered zero
times on a set of positive measure".

Measured convergence at K = 30, 80, 160:

| family | λ_min → | exact uncovered measure |
|---|---|---|
| {1,…,13} | 1.0794 → 1.0001 → **1.00000** | **0** |
| {1,…,11,13,24} | 1.0981 → 1.0008 → **1.00000** | **0** |
| random spread | 1.2628 → 0.3999 → 0.0265 | > 0 |
| AP d = 8 | 1.3691 → 0.7958 → 0.1499 | > 0 |

## The impossibility

λ_min sees only the **essential** infimum, so the LP certifies exactly the
statement *"the uncovered set has positive measure"* — strictly stronger
than LRC(14), which allows the lonely set to be null. And that stronger
statement is **false**: all three tight families have uncovered measure
exactly 0 while LRC(14) holds for them, with gap exactly 1/14 attained at
the six points p/14.

So the Delsarte LP cannot prove LRC(14) **even in principle**. It is not
that the bound is weak or the polynomial degree too small — the quantity it
computes is blind to measure-zero witnesses, and the extremal families are
precisely those whose only witnesses are measure-zero.

## What this adds

A clean impossibility for the strongest standard technique available, so
nobody spends a session on it. It also sharpens the picture from THM-1170:
the difficulty of LRC(14) is concentrated entirely in the null-set regime.
Every method in this programme that works with measures — Bonferroni
(THM-1095), the density bounds (THM-1155/1165), and now the LP — is
structurally unable to see the extremal families, because those families
have nothing to see at the level of measure. Methods that survive contact
with the tight families have all been *pointwise*: the classical sieve, the
essential-region criterion, and the beat-frequency structure.
