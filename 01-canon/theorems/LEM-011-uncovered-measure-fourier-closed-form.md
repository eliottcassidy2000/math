---
id: LEM-011
title: The exact Fourier closed form for the uncovered-measure function 𝒲 on T^{k-1} — 𝒲̂(n) = (6/7)^z · [Π_{i: n_i≠0} b(n_i)] · Q(N), where b(m)=(1−e(m/7))/(2πim), Q(0)=6/7, Q(N)=(e(−N/7)−1)/(2πiN), z=#{zeros}, N=Σn_i. This is the a-priori "𝒲̂-decay constant" shared by opus-S157's density-floor tail rate and THM-664's Weyl grid-resonance sum: |𝒲̂(n)| ≤ (6/7)^z Π|sin(πn_i/7)|/(π|n_i|)·|Q(N)|, vanishing whenever 7|n_i (nonzero) or 7|N — replacing both numerical certifications with a closed form
status: PROVED (elementary Fourier computation, below) + VERIFIED (direct T^2/T^3 Fourier integrals match to ≤2e-5 = integration error, k=3,4; 𝒲̂(0)=(6/7)^k exact; the 7-commensurate zeros 𝒲̂(7,0)=𝒲̂(1,6)=0 exact). Consequence: opus-S157's mixed-variation constants V_j (|Ĝ^j(a,b)|≤V_j/|ab|) and THM-664's grid resonance sum R=Σ_{n≠0,Vmax|n·e}𝒲̂(n) are now A-PRIORI (explicit closed form, no grid certification). Instantiates LEM-007's master law Ŵ(m)=Σ_{Σn=0,Σn·e=m}Πĝ as the T^{k-1} coefficient
source: mac-mini-2026-07-08-S59
depends_on:
  - THM-657   # W = uncovered measure = 1 − coverage of the k arcs of length 1/7
  - LEM-007   # the overlap Fourier master law (ĝ_0=6/7, ĝ_n=−ĥ_n); this is its T^{k-1} form
related:
  - THM-664   # Weyl first-moment: E_grid[W]=(6/7)^k + Σ_{Vmax|n·e}𝒲̂(n) — bounded a-priori by this
  - LEM-009   # the density-floor tail (opus-S157 rate |D3−D3_∞|≤C/(pd)); C now a-priori
  - THM-661   # the covering-moment density floor (the tail leg this certifies)
---

# LEM-011 — The uncovered-measure Fourier closed form

## Setup

For phases `0 = φ_0, φ_1, …, φ_{k−1}` on the circle, let `A_i = [φ_i, φ_i + ℓ)`, `ℓ = 1/7`, and
`𝒲(φ_1,…,φ_{k−1}) := 1 − meas(⋃_i A_i) = ` the **uncovered measure** (THM-657's `W`, with the
observer phase `φ_0` pinned at `0`). `𝒲 : T^{k−1} → [0, 6/7]` is continuous. Its Fourier
coefficients `𝒲̂(n) = ∫_{T^{k−1}} 𝒲(φ) e(−n·φ) dφ` (`n·φ = Σ_{i≥1} n_i φ_i`) have a closed form.

## The closed form (PROVED)

Write `𝒲(φ) = ∫_0^1 ∏_{i=0}^{k−1}(1 − 1_{A_i}(x)) dx`, and note `1 − 1_{A_i}(x) = 1 − 1[φ_i ∈
(x−ℓ, x]]`. Then

`𝒲̂(n) = ∫_0^1 (1 − 1[0 ∈ (x−ℓ,x]]) ∏_{i≥1} F_i(x) dx`,  `F_i(x) = ∫_0^1 (1 − 1[φ∈(x−ℓ,x]]) e(−n_i φ) dφ`.

For `n_i = 0`: `F_i = 1 − ℓ = 6/7`. For `n_i ≠ 0`: `∫_0^1 e(−n_iφ)dφ = 0`, so
`F_i(x) = −∫_{x−ℓ}^{x} e(−n_iφ)dφ = e(−n_i x)·b(n_i)`, `b(m) := (1 − e(mℓ))/(2πim)`. Factoring the
constant `(6/7)` from each zero coordinate and collecting the `e(−n_i x)`:

`𝒲̂(n) = (6/7)^z · [∏_{i: n_i≠0} b(n_i)] · ∫_0^1 (1 − 1[x∈[0,ℓ)]) e(−N x) dx`,  `N := Σ_{i≥1} n_i`, `z := #{i : n_i = 0}`.

The last integral is `Q(N) := 6/7` if `N=0`, else `(e(−Nℓ) − 1)/(2πiN)`. Hence

> **`𝒲̂(n) = (6/7)^z · [∏_{i: n_i≠0} b(n_i)] · Q(N)`,  `b(m) = (1−e(m/7))/(2πim)`, `Q(N) = ĝ_N`.** ∎

(`Q(N) = b(−N)` for `N≠0`; this is LEM-007's master law `Ŵ(m) = Σ_{Σn=0, Σn·e=m} ∏ ĝ` read as the
`T^{k−1}` coefficient, the `Q(N)` factor being the pinned `φ_0`-coordinate's `ĝ_{−N}`.)

## The a-priori decay bound

`|b(m)| = |1−e(m/7)|/(2π|m|) = |sin(πm/7)|/(π|m|) ≤ 1/(π|m|)`, and `= 0` iff `7 | m`. Likewise
`|Q(N)| = |sin(πN/7)|/(π|N|) ≤ 1/(π|N|)` (`N≠0`), `Q(0)=6/7`. Therefore

> **`|𝒲̂(n)| ≤ (6/7)^z · ∏_{i: n_i≠0} \frac{|sin(πn_i/7)|}{π|n_i|} · |Q(N)|`,**  and
> **`𝒲̂(n) = 0` unless every nonzero `n_i` and (if `N≠0`) `N` is coprime to `7`** (7-commensurate support).

This is fully explicit — no grid, no numerical certification. Two structural corollaries:

- **Single-frequency `1/m²` decay.** For support-1 `n` (only `n_i = m ≠ 0`), `N = m`, so
  `|𝒲̂(n)| = (6/7)^{k−2} |sin(πm/7)|² / (π²m²) = O(1/m²)`. This `1/m²` is what powers opus-S157's
  `O(1/(pd))` tail rate.
- **Balanced (doubly-balanced) resonances.** `N = 0 ⟹ Q = 6/7` (no extra decay), so the
  `Σn_i = 0` coefficients `|𝒲̂(n)| ≤ (6/7)^{z+1}∏|b(n_i)|` are the slow directions — exactly
  LEM-007's doubly-balanced resonances and the tent `T(V')` of THM-664.

## Closes both certifications

- **opus-S157 / LEM-009 (density-floor tail).** The tail rate `|D3(E_{d,p}) − D3_∞| ≤ C/(pd)` used
  `|Ĝ^j(a,b)| ≤ V_j/|ab|` with `V_j` *numerically* certified (0.28/0.16/0.10). Now `Ĝ^1(a,b) =
  Σ_{n: Σ i n_i = a, n_L = b} 𝒲̂(n)` is an explicit sum of the closed form (and `Ĝ^{j}` its
  `j`-fold convolution), so each `V_j` — hence `C` — is **a-priori** (bounded by the `∏|sin|/(π|n_i|)`
  product summed over the pullback fibre; the `1/m²` single-frequency decay gives the `Σ_κ 1/(pd)`
  convergence directly). The density-floor tail loses its last numerical certification.
- **THM-664 (Weyl grid resonance).** `E_grid[W] = (6/7)^k + Σ_{n≠0, Vmax|n·e} 𝒲̂(n)`; each summand is
  now a-priori, and the `7`-commensurate vanishing + `1/(|n_i||N|)` decay bound the resonance sum
  explicitly (small once the small-`n` resonances are counted — the same additive structure as the
  density-floor tail).

The **decay constant is one object** — the closed form `𝒲̂(n)` — shared by the continuous tail
(`n·e = 0`) and the grid sum (`Vmax | n·e`); proving it a-priori discharges both certifications at
once. File: `04-computation/lrc14_What_closedform_macmini_S59.py` (+ `.out`).
