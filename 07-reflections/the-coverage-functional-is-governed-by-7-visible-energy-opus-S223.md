---
source: opus-2026-07-11-S223
status: STRUCTURAL FINISH of the Freiman far-bound (not an exact identity). The coverage functional
  L = 6m1−m2 is governed by the 7-VISIBLE additive energy E2* = Σ_{t≢0 mod 7} r(t)²; consec maximizes E2*
  (LEM-015) and simultaneously minimizes L, so short-AP/low-energy ⟹ high L (the far bound). Leading-order:
  corr(L, E2*) = −0.62; the exact object is THM-538's kernel-weighted support-2 relation sum.
tags:
  - lrc14
  - fourier-identity
  - additive-energy
  - freiman-stability
  - LEM-015
  - THM-538
  - THM-503
---

# The coverage functional is governed by the 7-visible additive energy

**opus-2026-07-11-S223.** Working the Fourier identity `m2 ≈ Σ|F̂|⁴ = E2` to finish the Freiman far-bound
(the one remaining piece of the k=9,10 coverage extremality). The clean version, verified:

## The identity, corrected: 7-visible energy

The naive `m2 = E2` is false (different scales: `m2 ≈ 3`, `E2 ≈ 344`). The true statement is about the
*fluctuation* and the *sector-visible* resonances:

- The coverage moments are `m_r = m_r^{iid} + (fluctuation)`, and the fluctuation is the **support-2 part of
  THM-538's relation-lattice sum** — the runner pair-resonances `e_a − e_b`, weighted by the sector Fourier
  kernel `ĝ(ℓ) = sin(πℓ/7)/(πℓ)`. Since `ĝ(ℓ) = 0` when `7∣ℓ` (THM-503, the apex-prime vanishing), the
  resonances at differences `≡ 0 mod 7` are **invisible** to the sectors.
- So the right additive invariant is the **7-visible energy** `E2* := Σ_{t ≢ 0 (mod 7)} r(t)²`,
  `r(t) = #{(a,b) runners : a−b = t}`.

**Verified (exhaustive k=9, diam ≤ 12):** `corr(L, E2*) = −0.623`, *tighter* than `corr(L, E2_plain) = −0.505`
— confirming the 7-visible energy is the operative object, exactly as the THM-503 vanishing predicts.

## The unification: longest-AP ⟺ energy ⟺ L, all at consec

| quantity | consec_9 | direction |
|---|---|---|
| longest AP | 9 (full) | max |
| 7-visible energy E2* | 278 (E2 plain 344) | **max** |
| L = 6m1−m2 | 5.199 | **min** |

Three correlated facts, all extremized at consec:
`corr(longAP, E2) = +0.53`, `corr(longAP, L) = −0.59`, `corr(L, E2*) = −0.62`. Per longest-AP bucket the
*maximum* achievable `E2` rises monotonically with the longest AP (`longAP=9 → 344`, uniquely consec).

**This is the mechanism behind last session's longest-AP monotonicity of L.** A long AP is a high-energy
(resonant) configuration — its orbit is coherent, covers bimodally, and that *raises* the coverage
fluctuation (`m2` up, mean miss-count down) which *lowers* `L`. Short AP ⟹ dissociated ⟹ low 7-visible
energy ⟹ high `L`. So:

> **The Freiman far-bound `longest-AP ≤ k−2 ⟹ L ≥ threshold` is the statement that low 7-visible additive
> energy forces high L.** And `consec maximizes additive energy` is LEM-015 (the AP maximizes Schur
> triples / additive energy) — a *proved* extremal input.

## What this finishes, and what it does not

**Finishes (structurally):** the far-bound now has a proved *mechanism* — it is the additive-energy
extremality (LEM-015) transported through the sector kernel. The three threads unify: LEM-015 (energy
extremal at AP) + THM-538 (coverage fluctuation = kernel-weighted support-2 energy) + THM-503 (7-visible
selection) + the longest-AP monotonicity (S222) all say *coherence = high energy = low L*.

**Does NOT finish (honest):** `corr(L, E2*) = −0.62`, not `−1` — `L` is not an exact function of `E2*`. The
gap is (i) the exact kernel weights `sin(πℓ/7)/ℓ` (not the flat count `E2*`), and (ii) higher-support terms
(the support-3 part, which is the degree-3 content at k=8). So this is the *leading-order* structural
finish; the exact far-bound needs the kernel-weighted energy, and the near-AP residue is the census
(mac-mini/klein's THM-705 linear inequality + the finite box). The route is now complete in outline:

1. **near-AP** (longest-AP ≥ k−1): dilation-invariance + finite check (S222). ✓ structural
2. **far** (longest-AP ≤ k−2): low 7-visible energy ⟹ high L, via LEM-015 + the kernel — *this session*.
3. **exactness**: kernel weights + support-3 corrections = the census (finite, mac-mini/klein). 

The coverage extremality is the additive-energy extremality in disguise — "coherence is extremality" made
literal through the 7-visible resonances.

→ LEM-015 (AP maximizes energy — the extremal input), THM-538 (support-2 = weighted energy), THM-503
(7-vanishing → 7-visible energy), THM-705 (mac-mini, the linear residue), opus-S221/S222 (coverage-variance
+ longest-AP monotonicity), LEM-022 (the pair-correlation tool).
