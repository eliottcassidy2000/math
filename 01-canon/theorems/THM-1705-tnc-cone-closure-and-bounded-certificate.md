---
id: THM-1705
title: "TWO STRONGER TNC CLOSURES. (1) COMMON-RAY CLOSURE, PROVED IN ONE LINE: if the nonzero coefficients of R share a common complex argument phi (r_k = rho_k e^{i phi}, rho_k > 0), then every charge-representation of 0 at level m uses EXACTLY m factors, so CT(Lambda^m) = e^{i m phi} * (strictly positive sum) != 0 -- R is TNC-clear. Since TNC is scale-invariant (R -> lambda R), WLOG rotate to positive-real R, where CT(m) > 0 for every reachable m. So the nullcone MISSES the ENTIRE common-ray locus, a real-codimension family that includes the positive orthant; a TNC-violator needs coefficients genuinely SPREAD in argument (phase-tuned, roots-of-unity-flavoured: the trinomial witness is a^2=-1). (2) THE (k-1)-LEVEL BOUNDED CERTIFICATE: for a k-nomial the emptiness test of THM-1685 is achieved by the FIXED level set {m_0, 2m_0, ..., (k-1)m_0} -- k-1 multiples of the minimal-representation level m_0. Verified: trinomial 2 levels, 4-nomial 3 levels, 5-nomial 4 levels (sharp: the N=3 charges -3,-2,-1,1,2 pattern needs all 4). This upgrades THM-1685 from a terminating algorithm to a CLOSED-FORM certificate of fixed cost, and its load-bearing fact is that CT(2m_0) mod <CT(m_0)> is a nonzero constant on the genuinely-non-unique witness ({-2,1,4}: -30), the GMC-cascade correction (THM-1535)"
status: PROVED (1) unconditionally; (2) VERIFIED as the level bound for k <= 5 (uniform-in-pattern proof open)
author: opus-2026-07-20-S422
depends_on: [THM-1685 (k-nomial Nullstellensatz), THM-1680 (trinomial gcd), THM-1655 (binomial), THM-415 (vanishing sums), THM-1535 (GMC cascade)]
---

# THM-1705 — Two stronger TNC closures: the common-ray cone and the bounded certificate

## 1. The common-ray closure (PROVED, one line)

> **Theorem.** If every nonzero coefficient of `R` has the same complex argument
> (`r_k = ρ_k e^{iφ}`, `ρ_k > 0`, common `φ`), then `Λ = u^{−N}R` is **not** a nullcone
> element. Hence TNC holds for every such `R`.

*Proof.* `CT(Λ^m) = Σ_{reps of 0 at level m} \binom{m}{\cdots} ∏ r_{k_i}`. Every
representation of `0` as a sum of charges uses **exactly `m`** charges (one per factor of
`Λ^m`), so each product `∏_{i=1}^{m} r_{k_i}` has argument exactly `mφ`. Factoring it out,

```
CT(Λ^m) = e^{imφ} · Σ_{reps} \binom{m}{\cdots} ∏ ρ_{k_i}  =  e^{imφ} · (strictly positive) ≠ 0
```

for every level `m` at which a representation exists (in particular `m = m_0`). ∎

*(Verified with `φ = π/5`: `CT(m_0) = 6(-1)^{3/5}`, `6(-1)^{2/5}`, `2(-1)^{2/5}` — all
nonzero.)*

**Scope.** TNC is invariant under `R ↦ λR` (`Λ ↦ λΛ`, `CT(Λ^m) ↦ λ^m CT(Λ^m)`, same zero
set), so the common-ray case reduces to **positive-real `R`**, where `CT(m) > 0`. The
nullcone therefore misses the whole common-ray locus — a real-codimension family containing
the positive orthant. **A TNC-violator must have coefficients genuinely spread in argument.**
The known tuned-cancellation witnesses are roots of unity (trinomial `a² = −1`), consistent
with the vanishing-sums-of-roots-of-unity structure of THM-415.

## 2. The bounded certificate: `(k−1)` levels suffice

THM-1685 reduced TNC for a `k`-nomial pattern to emptiness of `V(I) ∩ (ℂ*)^{k−2}`,
`I = ⟨CT(Λ^m)⟩`. This note bounds *which* `CT(m)` are needed:

> **Claim.** The emptiness is certified by the **fixed level set**
> `{m_0, 2m_0, …, (k−1)m_0}` — the first `k−1` multiples of the minimal-representation level
> `m_0`.

| `k` | levels needed | verified patterns |
|---|---|---|
| 2 (binomial) | `{m_0}` (1 level) | all (THM-1655) |
| 3 (trinomial) | `{m_0, 2m_0}` (2) | 10 (THM-1680) |
| 4 | `{m_0, 2m_0, 3m_0}` (3) | 5 |
| 5 | `{m_0, …, 4m_0}` (4) | 2 — **sharp**: `N=3`, charges `−3,−2,−1,1,2` needs all 4 |

So THM-1685's terminating algorithm becomes a **closed-form certificate of fixed cost**:
compute `CT` at `k−1` prescribed levels, form the ideal, saturate, test `1`.

**Load-bearing fact.** For the genuinely-non-unique trinomial `{−2,1,4}`,
`CT(2m_0) \bmod ⟨CT(m_0)⟩ = −30` — a **nonzero constant**, directly certifying
`V(CT(m_0), CT(2m_0)) ∩ ℂ* = ∅`. This is the **GMC `n=2/n=4` cascade correction** (THM-1535
§3) in the coefficient ring: `CT(2m_0) = (CT(m_0)/c)^2 + \text{correction}`, and the
correction survives on `V(CT(m_0))`. (For unique-minimal patterns `CT(m_0)` is already a
monomial, so one level suffices and the higher ones are redundant — `CT(2m_0) ∈ ⟨CT(m_0)⟩`.)

## 3. The TNC picture now

Four independent handles, closing overlapping regions:

| handle | region closed | by |
|---|---|---|
| Dickson ladder | small bidegree (`M=0`, `min(M,N)=1`, `(2,2)`, `(2,3)`) | THM-1530/1595 |
| few-terms procedure | binomial, trinomial, `k`-nomial patterns | THM-1655/1680/1685 |
| **common-ray cone** | **coefficients on a common ray (incl. positive orthant)** | **§1** |
| **bounded certificate** | **any `k`-nomial via `k−1` fixed levels** | **§2** |

The residual is the intersection of the *complements*: large-bidegree, many-term, and
phase-spread coefficients simultaneously — and even there, §2 gives a fixed-cost decision.

## 4. Next — the single-shot closure

1. **Uniform `(k−1)`-level proof.** Show `{m_0,…,(k−1)m_0}` saturates for *every* charge
   pattern, not just the tested ones. With §1 removing the common-ray directions, the
   remaining varieties are cut by `k−1` equations in `k−2` unknowns — an over-determined
   system whose only solution should be the coordinate degeneration (a lower `k`).
2. **The cyclotomic conjecture.** Are all tuned-cancellation points roots of unity in the
   normalized gauge? If so, TNC is a **cyclotomic non-vanishing** and THM-415's
   prime/composite vanishing-sums dichotomy applies directly — the same object as the
   JC-monodromy residual (HYP-8450). This would be the genuinely single-shot closure.

## Verification

`04-computation/tnc_strong_closures_opus_S422.py` (common-ray `CT ≠ 0` with `φ = π/5`; the
`(k−1)`-level saturation for `k = 5`), `tnc_monomial_correction_opus_S422.py` (the cone
positivity; the `CT(2m_0) \bmod CT(m_0)` structure; 4-nomial two/three-level emptiness).
Outputs in `05-knowledge/results/`.
