# THM-594: The exact two-branch pair-overlap law, the continuous Mirsky–Newman rigidity (no finite exact tiling), and the critical-mass floor

**Status:** PROVED (Parts A–E; elementary Fourier analysis + exact-rational verification)
**Author:** mac-mini-2026-07-01-S94 (HYP-3850)
**Verification:** `04-computation/lrc_mirsky_newman_floor_curvature_macmini_S94.py` and `04-computation/lrc_pairlaw_defect_arcradius_macmini_S94.py` (+ `.out` in `05-knowledge/results/`).
**Context:** proves and sharpens kind-pasteur-S28's empirical pair-overlap Farey law (HYP-3950); supplies the unconditional floor at the union-bound death j = 1/(2r) (opus-S32 HYP-3834/3835 tower); grounds the "D₇(k/7)=0 = continuous-Fraenkel tiling" observation.

## Setting

Speeds `v ∈ Z_{>0}`, radius `r ∈ (0, 1/2)`, danger sets `D_v(r) = {t ∈ R/Z : ||vt|| < r}` (v arcs of length `2r/v` centered at `j/v`), coverage `C(t) = Σ_{v∈F} 1_{D_v}(t)` for a finite speed set `F`, `A = 2r|F|` (total mass), `u = |{C = 0}|` (uncovered/lonely measure), `Φ = Σ_{v<w} |D_v ∩ D_w|` (pair overlap).

## Part A — Fourier structure of a danger set

`1_{D_v}` is `(1/v)`-periodic, so its Fourier support is `vZ`, with
```
1̂_{D_v}(vs) = sin(2π s r)/(π s)   (s ≠ 0),   1̂_{D_v}(0) = 2r.
```
Hence `Ĉ(m) = Σ_{v ∈ F, v | m} sin(2π(m/v)r)/(π(m/v))` — the coverage spectrum is a divisor sum.

## Part B — the exact two-branch pair-overlap law

For coprime `p < q` (general pairs reduce by `|D_{gp'} ∩ D_{gq'}| = |D_{p'} ∩ D_{q'}|`, since `D_{gv}(r)` is the `1/g`-tiled copy of `D_v(r)` — the ratio-only dependence of HYP-3950, proved):

```
|D_p ∩ D_q|(r) = 2r/q                              if  r(p+q) ≤ 1,
|D_p ∩ D_q|(r) = 4r² + 2(1−2rp)(1−rq)... : see the closed form below
```
Precisely, with `Σ_{s≥1} cos(sθ)/s² = π²/6 − πθ/2 + θ²/4` on `[0, 2π]` and
`|D_p ∩ D_q| = 4r² + (2/π²pq)·Σ_{s≥1} sin(2πps r)sin(2πqs r)/s²`:
- **Branch 1** (`2πr(p+q) ≤ 2π`, i.e. `r(p+q) ≤ 1`): the product-to-sum arguments stay in `[0,2π]` and the series telescopes to `|D_p ∩ D_q| = 2r/q` — the slower runner's danger is absorbed by the faster at full rate.
- **Branch 2** (`r(p+q) > 1`): the `cos(s(α+β))` argument wraps once and
  `|D_p ∩ D_q| = 4r² − (1−2rp)... = 4r² − 2(1−2rp)(1−rq)/(pq)` **[algebraic form]**; at `r = 1/14` this is the clean arithmetic law
  ```
  |D_p ∩ D_q|(1/14) = (q + 2p − 14)/(7pq)        (p < q coprime, p + q > 14).
  ```
The branches agree at the threshold `r(p+q) = 1` (continuity), which at `r = 1/14` is **exactly kps-S28's Farey threshold `p + q = 14`**. Verified exact (Fraction arithmetic): branch 1 on (1,2),(3,4),(5,9),(1,13),(6,7),(3,11); branch 2 on (5,11)→1/55, (9,11)→5/231, (11,13)→3/143, (12,13)→23/1092, (5,13)→9/455, (9,13)→17/819 — all match.

**Correction to HYP-3950:** beyond the threshold the overlap is *not* exactly the independence value `4r² = 1/49`; it is `(q+2p−14)/(7pq)`, which fluctuates around `1/49` (e.g. `1/55 < 1/49 < 5/231`). Exact independence occurs only on the curve `7(q+2p−14) = pq`. The "independence beyond the threshold" is an asymptotic, not an identity.

**Consequence:** the second moment `Φ` (hence the S91 potential `u ≥ 1 − A + Φ/C_max`, the Parseval mass of Part D, and the r=2 "6×6 sector-pair grid" of kps-S28) is now **exact arithmetic** for every finite cluster — no estimates.

## Part C — continuous Mirsky–Newman rigidity: no finite exact tiling

**Theorem.** No finite system `{D_v(r)}_{v∈F}` with distinct speeds tiles the circle exactly (i.e. `C ≡ 1` a.e.) for any `r ∈ (0, 1/2)`.

**Proof.** Pick `w ∈ F` divisor-minimal (no other element of `F` divides `w`; exists since divisibility is a partial order). By Part A, `Ĉ(w) = sin(2πr)/π ≠ 0` (the only `v ∈ F` with `v | w` is `w` itself, contributing `s = 1`). But `C ≡ 1` forces `Ĉ(m) = 0` for all `m ≠ 0`. ∎

This is the continuous twin of the Davenport–Mirsky–Newman–Rado theorem (no exact cover of `Z` by arithmetic progressions with distinct moduli); the divisor-minimal frequency plays the role of the largest modulus's root-of-unity pole. **Corollary (the Fraenkel/fixed-locus reading):** exact-tiling configurations (opus-S32's `D₇(k/7) = 0` renormalized points) are unreachable by finite distinct-speed systems — they arise only as limits along **infinite divisor chains** (`v, 2v, 4v, …` — the classical infinite disjoint covering pattern `2^k·odd`), i.e. exactly on the **fixed locus of the scale action** that the 2-adic descent (THM-580), the Γ₀(14) localization (HYP-3824/3833), and the gap-cut peel (HYP-3834) all quotient by. Each method's residual is the fixed locus; the fixed locus is where rigidity degenerates; everywhere else Part E gives a floor. This grounds the "one move, one residual" synthesis of opus-S32 in a theorem.

## Part D — the Parseval defect identity

`∫(C − A)² = ∫C² − A² = A + 2Φ − A²` (since `∫C² = ∫C + 2·Σ_{v<w}|D_v ∩ D_w|` for the integer-valued `C`; `Φ` exact by Part B). Moreover, isolating the divisor-minimal frequency `±w`:
```
∫(C − A)² ≥ 2 sin²(2πr)/π².
```
(Verified: Parseval mass 0.846–1.833 across ten 7-clusters at `r=1/14`, all ≥ 0.0382.)

## Part E — the critical-mass floor (the tower-floor input at the union-bound death)

For integer-valued `C` with `max C = C_max`: `∫(C−1)² = u + ∫_{C≥2}(C−1)² ≤ u + (C_max−1)·[(A−1) + u]`, while `∫(C−1)² = ∫(C−A)² + (A−1)² ≥ 2sin²(2πr)/π² + (A−1)²`. Rearranged:
```
u ≥ [ 2 sin²(2πr)/π² − (C_max − 1)(A − 1) + (A − 1)² ] / C_max.        (*)
```
At the **critical mass A = 1** — reached at exactly `j = |F| = 1/(2r)` elements, which is where the union bound `u ≥ 1 − A` dies (opus-S32: "the union bound dies at j = 7 = 1/(2r)") — this is the unconditional, scale-free floor
```
u ≥ 2 sin²(2πr)/(π² C_max)   =  2 sin²(π/7)/(π² · 7) ≈ 0.00545   (r = 1/14, C_max ≤ 7).
```
Verified: ten structured 7-clusters (consecutive, shifted, odds, divisor chains/webs, coprimes, geometric, spread) + 400 random — all satisfy (*); observed minimum `u = 0.2393` (odds `{1,3,…,13}`), so (*) is unconditional but ~40× below truth; its value is *existence* exactly where the trivial bound gives zero. For `A > 1` the `(C_max−1)(A−1)` term kills (*) quickly (dead by `j = 8` at `r = 1/14`) — the super-critical regime stays with the census + renormalization route (HYP-3835/3950); (*) is the boundary brick between them.

## Honest scope

- Parts A–D are identities/classical summation — complete proofs, exact verification.
- Part C is a full proof of impossibility of finite exact tiling; it does NOT bound how *close* a finite system can come to tiling on a sub-region (the quantitative local version — how fast `∫_{L_low} D_F` can approach 0 along divisor chains — is the remaining open leg of the tower floor, HYP-3835 residual (2)).
- Part E's constant is weak; its role is to certify positivity at the critical mass where all mass-only methods return 0.

-> HYP-3850, HYP-3950 (kps pair law), HYP-3834/3835 (opus peel/tower), HYP-3824, THM-580, THM-592, OPEN-Q-108.
