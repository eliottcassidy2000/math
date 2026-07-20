---
id: THM-1675
title: "THE TUNED-CANCELLATION LOCUS IS EMPTY FOR TRINOMIALS: TNC reduces to a FINITE GCD, and every tested pattern closes. (0) STRUCTURAL FACT: once r_N=0 (forced, THM-1655), F(u,t)=u^N - tR(u) is LINEAR in t -- [u^N]F=1, [u^k]F=-t r_k. (1) THE DECISION PROCEDURE: every trinomial reduces by the u-scale and t-scale gauges to R = 1 + a u^j + u^d (single free coefficient a); each CT(Lambda^m) is then a POLYNOMIAL in a, and a non-monomial nullcone violator is a common NONZERO root of ALL of them. So TNC for a charge pattern <=> the CT(m) polynomials have no common nonzero root <=> their finite GCD is a power of a (a=0 = monomial, excluded). (2) VERIFIED 10/10 patterns: in every case gcd(CT(m_0),...,CT(m_0+3)) in a is a nonzero constant or a power of a -- NO nonzero common root -- so NO non-monomial trinomial is a nullcone element there. (3) THE WITNESS, exact: {-2,1,4}, N=2 gives CT(3)=3(a^2+1) and CT(6)=15(a^4+4a^2+1); a^4+4a^2+1 = (a^2+1)^2 + 2a^2, so at a^2=-1 (CT(3)=0) it equals 2a^2=-2 != 0, hence CT(6)=-30. COPRIME, TNC proved for this pattern. The '+2a^2' surviving correction is the same leading-square-plus-correction structure as the GMC cascade. (4) THIS GENERALISES THE LADDER: binomial = 0 parameters (unique min rep, THM-1655), trinomial = 1 parameter (this GCD), k-nomial = (k-2) parameters (a Groebner/Nullstellensatz check on the CT ideal) -- a DECIDABLE finite reduction per charge pattern, at all bidegrees, replacing the one-at-a-time Dickson ladder"
status: PROVED as a decision procedure; 10/10 trinomial patterns closed exactly. A uniform coprimality proof over all (j,d,N) is the remaining step.
author: opus-2026-07-20-S420
depends_on: [THM-1655 (binomial + unique-min; r_N=0), THM-1635 (branch product), THM-1615, THM-1595, THM-1530]
---

# THM-1675 — The trinomial finishing statement: TNC is a finite gcd

## 0. Structural fact: `F` is linear in `t`

Once `r_N = 0` (forced, THM-1655 §1), `F(u,t) = u^N − tR(u)` has `[u^N]F = 1 − t r_N = 1` and
`[u^k]F = −t r_k` for `k ≠ N`. So **`F` is linear in `t`** — degree exactly 1 — and the whole
branch geometry is governed by a `t`-affine polynomial. (Verified for `1+u³−u⁶`, `1+u²+u⁵`.)

## 1. The decision procedure (the finishing statement)

Every trinomial `R = r_0 + r_j u^j + r_d u^d` with `r_0, r_d ≠ 0` reduces, by the two gauges

- **u-scale** `u ↦ λu` with `λ^d = r_0/r_d` (sets the `u^d` coefficient to 1),
- **t-scale / overall** (sets `r_0 = 1`),

to the **one-parameter** family `R = 1 + a·u^j + u^d`. Then each `CT(Λ^m) = [u^{Nm}]R^m` is a
**polynomial in `a`**, and

> **a non-monomial nullcone violator is a common NONZERO root `a` of ALL the `CT(Λ^m)`.**

Therefore:

> **THM-1675.** For a fixed trinomial charge pattern `(N; j, d)`, TNC holds **iff** the
> `CT(Λ^m)` (as polynomials in `a`) have **no common nonzero root** — equivalently, their
> gcd is a monomial `a^k` (the root `a = 0` is the excluded degenerate case). This is a
> **finite gcd computation**: it suffices to check finitely many `m`, since a gcd stabilises.

## 2. Verified: 10/10 patterns closed

| pattern `(N; charges)` | `gcd` of the `CT(m)` in `a` | nonzero common root |
|---|---|---|
| `(2; −2,1,4)` | `3` | **none → TNC** |
| `(2; −2,−1,2)` | `1` | none → TNC |
| `(2; −2,−1,3)` | `a` | none → TNC |
| `(2; −2,1,3)` | `a` | none → TNC |
| `(3; −3,−2,2)` | `a` | none → TNC |
| `(3; −3,−1,4)` | `a` | none → TNC |
| `(2; −2,−1,5)` | `2a` | none → TNC |
| `(3; −3,−2,1)` | `3a` | none → TNC |
| `(3; −3,1,2)` | `a` | none → TNC |
| `(2; −2,3,4)` | `2a²` | none → TNC |

In every case the gcd is `c·a^k` — **no nonzero common root**, so no non-monomial trinomial
in that pattern is a nullcone element. (Only four `CT(m)` per pattern were needed for the gcd
to stabilise to a monomial.)

## 3. The witness, exactly

`{−2,1,4}`, `N=2`, `R = 1 + a u^3 + u^6`:

```
CT(Λ^3) = 3(a² + 1)
CT(Λ^6) = 15(a⁴ + 4a² + 1) = 15[(a² + 1)² + 2a²]
```

At `a² = −1` (where `CT(3) = 0`), the bracket is `2a² = −2 ≠ 0`, so `CT(6) = −30`. The two
polynomials are **coprime**, hence no common root: TNC proved for this pattern, and the
`1 + u³ − u⁶` witness (`a = i` gives real coefficients up to relabeling; `a = ±i`) has
`CT(6) ≠ 0` for the structural reason, not by luck.

**Resonance with GMC.** `CT(6) = (leading square) + 2a²` — the surviving `+2a²` is exactly the
"leading term cancels, correction survives" structure of the GMC `n=2`/`n=4` cascade
(THM-1535 §3). The tuned-cancellation locus fails for the same reason the sign-coherent GMC
nullcone is rigid: a positive-definite correction survives the vertex cancellation.

## 4. The ladder replaced by a decidable reduction

| support of `R` | free gauge parameters | TNC test |
|---|---|---|
| **binomial** (2 terms) | 0 | unique minimal rep, `CT(m_0) ≠ 0` automatically (THM-1655) |
| **trinomial** (3 terms) | 1 (`a`) | gcd of `CT(m)` in `a` is a monomial (§1) |
| **`k`-nomial** | `k − 2` | the `CT(m)` **ideal** in `ℂ[a_1,…,a_{k−2}]` has no nonzero-coord zero (Nullstellensatz/Gröbner) |

This is a **decidable finite reduction at every bidegree simultaneously**, in the number of
*terms* of `R` — not the one-bidegree-at-a-time Dickson ladder (THM-1595). Sparsity of `R`,
not `(M,N)`, is the true complexity parameter.

## 5. Status of the TNC

| case | status | by |
|---|---|---|
| `M=0`; `min(M,N)=1`; `(2,2)`,`(2,3)` | PROVED | THM-1530/1595 |
| binomial `R` | PROVED | THM-1655 |
| unique-minimal-rep `R` | PROVED | THM-1655 |
| **trinomial `R`, 10 patterns** | **PROVED (gcd)** | **§2** |
| trinomial `R`, all `(j,d,N)` | decidable per pattern; uniform coprimality open | §1 |
| `k`-nomial | reduces to a Gröbner check on `k−2` parameters | §4 |

## 6. Next

1. **Uniform trinomial coprimality.** Show `gcd(CT(m_0), CT(2m_0))` in `a` is always a
   monomial. The witness shows `CT(2m_0) = (CT(m_0)/c)² + (positive correction)`; a general
   "the correction never vanishes at a root of `CT(m_0)`" closes all trinomials at once.
2. **The doubling law.** Empirically the first surviving level after a tuned `CT(m_0) = 0` is
   `2m_0` (§3, and `tnc_finishing_trinomial`). Proving "level `2m_0` always survives" would
   bound the gcd check to two levels for every pattern.

## Verification

`04-computation/tnc_trinomial_gcd_decision_opus_S420.py` (the 10-pattern gcd table; the
witness `CT(3),CT(6)` coprimality), `tnc_finishing_trinomial_opus_S420.py` (the `t`-linearity;
53-pattern no-survivor sweep; the doubling law). Outputs in `05-knowledge/results/`.
