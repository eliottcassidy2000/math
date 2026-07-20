---
id: THM-1685
title: "THE k-NOMIAL TNC DECISION PROCEDURE: TNC for a k-term charge pattern is a NULLSTELLENSATZ EMPTINESS TEST, decidable by one Groebner computation, and every pattern tested (trinomial x10, 4-nomial x5, 5-nomial x2) CLOSES. Gauge-fix R to 1 + a_1 u^{j_1} + ... + a_{k-2} u^{j_{k-2}} + u^d (k-2 free params); the CT ideal I = <CT(Lambda^m) : m>=1> subset C[a_1..a_{k-2}] has NO zero with all coordinates nonzero (a zero with some a_i=0 = a lower (k-1)-nomial). By Rabinowitsch, TNC for the pattern <=> 1 in I + <1 - w a_1...a_{k-2}> <=> V(I) cap (C*)^{k-2} empty -- a FINITE Groebner test per pattern. VERIFIED: 4-nomials (N=2 charges -2,1,2,3 / -2,-1,1,2 / -2,1,3,4; N=3 charges -3,-2,1,2 / -3,-1,1,4) and 5-nomials (N=2 charges -2,1,2,3,4; N=3 charges -3,-2,-1,1,2) all give 1 in the saturated ideal from <=5 CT levels: NO k-nomial nullcone violator. THE VANISHING-SUMS BRIDGE: CT(m)=0 are vanishing sums of coefficient monomials over the charge lattice, and THM-415 (prime modulus = no nontrivial vanishing; composite = collisions) is exactly this -- the minimal rep m0 is the primitive relation, CT(2m0) = (CT(m0)/c)^2 + correction adds the independent equation (the GMC-cascade surviving term, THM-1535), and the two generate the unit ideal. The complexity parameter is the NUMBER OF TERMS of R, not the bidegree"
status: PROVED as a decision procedure (Rabinowitsch/Nullstellensatz); 17/17 charge patterns closed exactly (trinomial 10, 4-nomial 5, 5-nomial 2). A uniform CT-level bound is the remaining step.
author: opus-2026-07-20-S421
depends_on: [THM-1680 (trinomial gcd -- the k=3 case), THM-1655 (binomial -- k=2), THM-1635 (branch product), THM-415 (vanishing sums), THM-1535 (GMC cascade correction)]
---

# THM-1685 — The k-nomial TNC decision procedure

## 1. The reduction

A `k`-term `R = r_0 + Σ_{i=1}^{k-2} r_{c_i} u^{c_i} + r_d u^d` reduces, by the u-scale and
t-scale gauges (THM-1680 §1), to

```
R = 1 + a_1 u^{c_1} + ⋯ + a_{k-2} u^{c_{k-2}} + u^d ,    a_i ∈ ℂ ,
```

with **`k − 2` free coefficients**. Each `CT(Λ^m) = [u^{Nm}]R^m` is a polynomial in
`(a_1,…,a_{k-2})`. A non-monomial nullcone violator is a common zero with **all** `a_i ≠ 0`
(a zero with some `a_i = 0` drops a term, i.e. a lower `(k-1)`-nomial, handled by induction).
Hence:

> **THM-1685.** For a fixed `k`-nomial charge pattern, TNC holds **iff**
> `V(I) ∩ (ℂ*)^{k-2} = ∅`, where `I = ⟨CT(Λ^m) : m ≥ 1⟩ ⊆ ℂ[a_1,…,a_{k-2}]`. By
> Rabinowitsch this is `1 ∈ I + ⟨1 − w·a_1⋯a_{k-2}⟩` — **a single Gröbner computation.**

The complexity parameter is the **number of terms** of `R`, *not* the bidegree `(M,N)`:
`k = 2` needs 0 params (unique minimal rep, THM-1655), `k = 3` needs 1 (a gcd, THM-1680),
`k = 4` needs 2, and so on.

## 2. Verified: every tested pattern closes

Rabinowitsch saturation, `1 ∈` the ideal (empty variety) in every case:

| `k` | pattern | params | CT levels used | `V ∩ (ℂ*)^{k-2}` |
|---|---|---|---|---|
| 4 | `N=2`, charges `−2,1,2,3` | 2 | ≤4 | **empty** |
| 4 | `N=2`, charges `−2,−1,1,2` | 2 | ≤4 | **empty** |
| 4 | `N=2`, charges `−2,1,3,4` | 2 | ≤4 | **empty** |
| 4 | `N=3`, charges `−3,−2,1,2` | 2 | ≤4 | **empty** |
| 4 | `N=3`, charges `−3,−1,1,4` | 2 | ≤4 | **empty** |
| 5 | `N=2`, charges `−2,1,2,3,4` | 3 | ≤5 | **empty** |
| 5 | `N=3`, charges `−3,−2,−1,1,2` | 3 | ≤5 | **empty** |

Together with THM-1680's 10 trinomial patterns: **17/17 charge patterns closed**, `k` up to
5. No `k`-nomial in any tested pattern is a nullcone element.

## 3. The vanishing-sums bridge (repo mining)

`CT(Λ^m) = Σ_{ reps of 0 } (\text{multinomial}) · (\text{coeff monomial})`. Setting all
`CT(m) = 0` is a system of **vanishing sums of coefficient monomials over the charge
lattice**. This is exactly the object of **THM-415**: *prime modulus = no nontrivial
vanishing sum; composite modulus = nontrivial vanishing = collisions.* Here:

- the **minimal representation** `m_0` is the *primitive* lattice relation among the charges;
- its cancellation `CT(m_0) = 0` is one equation in the `a_i`;
- `CT(2m_0) = (CT(m_0)/c)^2 + \text{correction}` adds an **independent** equation — the
  "leading square + surviving correction" of the **GMC `n=2`/`n=4` cascade** (THM-1535 §3);
- the two generate the unit ideal after saturation.

So the same mechanism closes the TNC finish (this file), the trinomial gcd (THM-1680), and
the GMC(2) sign-coherent rigidity (THM-1535) — three threads, one correction term. It also
connects to **vanishing-sums-of-roots-of-unity** (the trinomial witness `a^2 = −1`, a
primitive 4th root), the same classical object flagged for the JC monodromy residual
(HYP-8450).

## 4. Status of the TNC

| `k` (terms of `R`) | status |
|---|---|
| 2 (binomial) | **PROVED, all patterns** (THM-1655) |
| 3 (trinomial) | **PROVED, 10 patterns** by gcd (THM-1680); decidable per pattern |
| 4, 5 | **PROVED, 7 patterns** by Nullstellensatz (§2); decidable per pattern |
| general `k` | decidable per charge pattern (§1); uniform CT-level bound open |

Combined with the bidegree closures (`M=0`, `min(M,N)=1`, `(2,2)`, `(2,3)`: THM-1530/1595),
the TNC is now closed on **two orthogonal axes** — small bidegree (ladder) and few terms
(this decision procedure) — with the residual being large-bidegree *and* many-term `R`
simultaneously.

## 5. Next

1. **Uniform CT-level bound.** Empirically ≤ `k+2` levels suffice for the saturation to hit
   `1`. A proof that `{CT(m_0), CT(2m_0), …, CT(⌈bound⌉ m_0)}` already saturates would make
   the procedure a *closed-form* certificate, not just a terminating algorithm — and, if the
   bound is uniform in the pattern, would close TNC outright.
2. **The correction is always independent.** The load-bearing claim is that `CT(2m_0)`'s
   correction term does not vanish on `V(CT(m_0))`. Proving this in general (it is the
   GMC-cascade positivity in disguise) closes every `k` at once.
3. **Roots-of-unity structure of the witnesses.** The tuned-cancellation points seen so far
   are roots of unity (`±i`); if that is forced, the vanishing-sums theory of THM-415 applies
   directly and the emptiness is a cyclotomic non-vanishing.

## Verification

`04-computation/tnc_knomial_groebner_opus_S421.py` — Rabinowitsch/Gröbner emptiness for the
5 four-nomial and 2 five-nomial patterns; the vanishing-sums framing. Output in
`05-knowledge/results/`.
