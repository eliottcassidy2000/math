---
id: THM-499
title: The spectral-reframe boundary — cycle counts are power sums, but H is strictly finer than the eigenvalue spectrum (cospectral tournaments with distinct H from n=6)
status: PROVED (exhaustive, exact integer spectral signatures + exact H; explicit cospectral/distinct-H witnesses)
source: kind-pasteur-2026-06-13-S5
depends_on:
  - THM-118   # c_k = tr(A^k)/k for k<=5 (cycle counts ARE power sums)
  - THM-498   # cycle-spectrum gaps = power-sum exclusions
related:
  - THM-029   # H-forbidden 7
  - THM-079   # H-forbidden 21
  - THM-133   # spectral-OCF chain
  - HYP-2492  # cycle gaps = skew-spectral exclusions (this REFINES its overreach to H)
---

# THM-499 — where the spectral reframe stops

This session's reframe (THM-498/HYP-2492) showed the cycle-count gaps are
**skew-spectral exclusions**: `c_k = tr(A^k)/k` for `k ≤ 5` (THM-118), so a
forbidden cycle-count is a non-realizable power-sum value. HYP-2492 conjectured
the OCF forbidden values `H ∈ {7,21}` are *also* power-sum exclusions. This file
draws the exact boundary: **they are not** — `H` is a strictly finer invariant
than the eigenvalue spectrum.

## Statement

Let `sig(T) = (tr A, tr A², …, tr Aⁿ)` (the exact integer power-sum signature ⟺
the characteristic polynomial of `A` ⟺ the eigenvalue spectrum).

> **THM-499.**
> 1. `H` (the Hamiltonian-path count = the OCF `I(Ω,2)`) **is** determined by
>    `sig(T)` at `n ≤ 5`, but **is NOT** at `n ≥ 6`: there exist **cospectral
>    tournaments with distinct `H`**. Witnesses (n=6, exhaustive): the spectral
>    class `sig = (0,0,18,28,30,120)` carries `H ∈ {25,29,33}`; `(0,0,12,12,10,48)`
>    carries `H ∈ {13,17}`; `(0,0,21,36,35,159)` carries `H ∈ {33,37}`.
> 2. Even the odd-cycle-count pair `(c3,c5)` does not determine `H` (6 split
>    classes at n=6).
> 3. Consequently the forbidden OCF values `{7,21}` are **NOT power-sum
>    exclusions** — unlike the cycle gaps (`c5=10`), they live in the
>    conflict-graph **disjointness** layer (`α₂` = # vertex-disjoint odd-cycle
>    pairs in `H = 1 + 2α₁ + 4α₂ + …`), a strictly higher invariant than the
>    spectrum.

*Status:* exhaustive over all `2^{C(n,2)}` labeled tournaments at `n=5,6`; `sig`
computed as exact integer traces, `H` exact (subset-DP Hamiltonian-path count);
witnesses explicit. (`04-computation/spectral_reframe_boundary_kps5.py` + `.out`.)

## The two kinds of forbidden value (the refined picture)

| forbidden value | layer | reason it's forbidden |
|---|---|---|
| `c5 = 10` (n=6) | **spectral** (power sum `tr A⁵`) | no tournament spectrum realizes `tr A⁵ = 50` (THM-498) |
| `H ∈ {7,21}` | **conflict-graph / disjointness** (`α₂`) | `α₁=3 ⟹ α₂≥1` etc. (THM-029/079) — NOT a spectral statement |

So "every forbidden value is a power-sum exclusion" is FALSE. The OCF invariant
hierarchy stratifies by *which* data the value depends on: the low cycle counts
are pure spectrum (`tr A^k`), but `H` injects the disjointness structure of the
odd-cycle conflict graph `Ω`, which cospectral tournaments can differ in.

## Why this matters (the reframe map gains an edge)

- `H` is **strictly finer** than the eigenvalue spectrum: cospectral tournaments
  (eigenvalue twins) are separated by `H`. This sharpens the repo's "H as a
  universal tournament fingerprint" (engineering domain 12): the fingerprint
  succeeds *because* it sees past the spectrum into `Ω`.
- It tells the spectral-exclusion program (HYP-2492, the route to proving `c5=10`
  and beyond) exactly how far it reaches: it can prove the **cycle-count** gaps
  spectrally, but the **H** gaps need a conflict-graph (independence-polynomial)
  argument — the two are genuinely different theorems, and conflating them was the
  HYP-2492 overreach this file corrects.
- The cospectral/distinct-H pairs are concrete new objects: minimal "eigenvalue
  twins, OCF-distinct" tournaments, the n=6 frontier of spectral non-determination
  for the Rédei–Berge / Hamiltonian-path statistic.

**Artifacts:** `04-computation/spectral_reframe_boundary_kps5.py` (+ `.out`).
Reflection: extends `efficiency-becomes-proof-kps4` (the spectral reframe and its
boundary).
